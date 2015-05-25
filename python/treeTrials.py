import scipy.io as io
import numpy as np
from scipy import stats
import runExperiments
import utilities
import sys

if __name__ == "__main__":

    # Parse the arguments
    args = utilities.parse_args(sys.argv)
    trials = args.get('trials', 1)
    write_results = args.get('write_results', False)
    alt = args.get('alt', False)
    run = args.get('run', 0)
    diffusion = args.get('diffusion', False)
    spy_probability = args.get('spy_probability')
    
    # Irregular infinite graph
    # xks = [np.arange(3,5), np.arange(3,6), np.arange(3,20,14)]
    # pks = [(0.5, 0.5), (0.5, 0.25, 0.25), (0.9, 0.1)]
    # max_times = [5, 4, 4]
    
    # Regular tree sizes:
    # d=3, T=12 => E[N] = 155
    # d=4, T=8 => E[N] = 120
    # d=5, T=7 => E[N] = 157
    xks = [np.array([3]) for i in range(1)]
    # pks = [(0.5, 0.5) for i in range(4)]
    pks = [(1.0) for i in range(1)]
    max_times = [16 for i in range(1)]
    # est_times = [6,8,10,12,14,16] # d=3
    est_times = [8,10,11,12] # d=3
    # est_times = [6,8,9,10,11] # d=4
    # est_times = [6,7,8,9] # d=5
    # est_times = [50,100,150,200] # d=2
    
    # xks = [np.arange(2,3)]
    # pks = [(1.0)]
    # max_times = [10]
    # additional_time = 30
    max_infection = 0 #min(xks) - 1
    additional_time = 0
        
    for (xk, pk, max_time) in zip(xks, pks, max_times):
        print('Checking xks = ',xk)
        degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
        max_infection = max_infection + 1
        
        print('Diffusion: ', diffusion)
        
        # Check if the tree is regular
        if isinstance(pk, list) == 1:
            if diffusion:
                print('Diffusion code')
                num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 4)[:2]
            else:
                num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 3)
            pd_rc, pd_jc = results
            
            if write_results:
                filename = 'results/regular_tree_results_d_' + str(xk) + '.mat'
                io.savemat(filename, {'pd_rc':np.array(pd_rc), 'pd_jc':np.array(pd_jc), 'num_infected':np.array(num_infected)})
        else:
            if diffusion:
                num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection,
                                                                    degrees_rv, method=4, p=0.5, 
                                                                    spy_probability = spy_probability,
                                                                    est_times = est_times)[:2]
                pd_ml, hop_distances, pd_spy, spy_hop_distances = results
                print('hop distances',hop_distances)
                if write_results:
                    if isinstance(pk, float):
                        xk_str = str(xk[0])
                        pk_str = str(pk)
                        filename = 'results/spies/regular_trees/results_' + xk_str + "_" + pk_str + '_spies_'+str(spy_probability) + '_run_' + str(run) + '.mat'
                    else:
                        xk_str = [str(i) for i in xk]
                        pk_str = [str(round(i,1)) for i in pk]
                        filename = 'results/spies/regular_trees/results_' + "_".join(xk_str) + "_" + "_".join(pk_str) + 'spies_'+str(spy_probability) + '_run_' + str(run) + '.mat'
                    io.savemat(filename,{'pd_ml':np.array(pd_ml), 'hop_distances':np.array(hop_distances),
                                         'pd_spy':np.array(pd_spy), 'spy_hop_distances':np.array(spy_hop_distances),
                                         'num_infected':np.array(num_infected)})
                    continue
            # Check for weighted spreading
            elif alt:
                num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 1)[:2]
            # Use regular spreading
            else:
                if additional_time:
                    num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 0,
                                                                      additional_time = additional_time)[:2]
                    print("additional pd: ", additional_pd)
                else:
                    print(max_infection, degrees_rv)
                    num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 0)[:2]
                pd_ml, additional_pd = results
            print("pd ML: ", pd_ml)
            print('num_infected: ', num_infected)

            if write_results:
                xk_str = [str(i) for i in xk]
                pk_str = [str(round(i,1)) for i in pk]
                if alt:
                    filename = 'results/irregular_tree_alt_results_' + "_".join(xk_str) + "_" + "_".join(pk_str) + '.mat'
                else:
                    filename = 'results/irregular_tree_results_' + "_".join(xk_str) + "_" + "_".join(pk_str) + '.mat'
                io.savemat(filename,{'pd_ml':np.array(pd_ml), 'num_infected':np.array(num_infected)})
