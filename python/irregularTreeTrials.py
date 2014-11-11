import scipy.io as io
import numpy as np
from scipy import stats
import runExperiments
import utilities

'''Run the irregular tree algorithm'''
if __name__ == "__main__":

    trials = 1
    write_results = False
    alt = True # Use the alternative method of spreading virtual sources?
    
    # Irregular infinite graph
    xks = [np.arange(2,4), np.arange(3,5), np.arange(2,5)]
    pks = [(0.5, 0.5), (0.5, 0.5), (0.5, 0.25, 0.5)]
    max_times = [13, 7, 10]
    max_infection = 2
    
    for (xk, pk, max_time) in zip(xks, pks, max_times):
        print('Checking xks = ',xk)
        degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
        if alt:
            num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 1)[:2]
        else:
            # Use regular spreading
            num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 0)[:2]
        print(pd_ml)

        if write_results:
            xk_str = [str(i) for i in xk]
            filename = 'results/irregular_tree_alt_results_' + "_".join(xk_str) + '.mat'
            io.savemat(filename,{'pd_ml':np.array(pd_ml), 'num_infected':np.array(num_infected)})
