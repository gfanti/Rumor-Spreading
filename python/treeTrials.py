import scipy.io as io
import numpy as np
from scipy import stats
import runExperiments
import utilities
import sys

if __name__ == "__main__":

    # Parse the arguments
    trials, write_results, alt = utilities.parse_args(sys.argv)
    
    # Irregular infinite graph
    # xks = [np.arange(3,5), np.arange(3,6), np.arange(3,20,14)]
    # pks = [(0.5, 0.5), (0.5, 0.25, 0.25), (0.9, 0.1)]
    # max_times = [5, 4, 4]
    xks = [np.arange(3,5) for i in range(9)]
    pks = [(0.1*i, 1-0.1*i) for i in range(1,10)]
    max_times = [5 for i in range(9)]
    # xks = [np.arange(2,3)]
    # pks = [(1.0)]
    # max_times = [10]
    # additional_time = 30
    max_infection = 2
    additional_time = 0
        
    for (xk, pk, max_time) in zip(xks, pks, max_times):
        print('Checking xks = ',xk)
        degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
        
        # Check if the tree is regular
        if isinstance(pk, list) == 1:
            num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 3)
            pd_rc, pd_jc = results
            
            if write_results:
                filename = 'results/regular_tree_results_d_' + str(xk) + '.mat'
                io.savemat(filename, {'pd_rc':np.array(pd_rc), 'pd_jc':np.array(pd_jc), 'num_infected':np.array(num_infected)})
        else:
            # Check for weighted spreading
            if alt:
                num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 1)[:2]
            # Use regular spreading
            else:
                if additional_time:
                    num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 0,
                                                                      additional_time = additional_time)[:2]
                else:
                    num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 0)[:2]
                pd_ml, additional_pd = results
            print("pd ML: ", pd_ml)
            print("additional pd: ", additional_pd)
            print('num_infected: ', num_infected)

            if write_results:
                xk_str = [str(i) for i in xk]
                pk_str = [str(round(i,1)) for i in pk]
                if alt:
                    filename = 'results/irregular_tree_alt_results_' + "_".join(xk_str) + "_" + "_".join(pk_str) + '.mat'
                else:
                    filename = 'results/irregular_tree_results_' + "_".join(xk_str) + "_" + "_".join(pk_str) + '.mat'
                io.savemat(filename,{'pd_ml':np.array(pd_ml), 'num_infected':np.array(num_infected)})
