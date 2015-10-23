import scipy.io as io
import numpy as np
from scipy import stats
import runExperiments
import utilities
import sys
from multinomial import *
from decimal import Decimal

if __name__ == "__main__":

    # Parse the arguments
    args = utilities.parse_args(sys.argv)
    trials = args.get('trials', 1)
    write_results = args.get('write_results', False)
    alt = args.get('alt', False)
    run = args.get('run', 0)
    diffusion = args.get('diffusion', False)
    spies = args.get('spies', False)
    spy_probability = args.get('spy_probability')
    
    if diffusion:
        delay_parameter = args.get('delay_parameter')
    
    # ----Irregular tree graph-----
    # xks = [np.arange(3,5), np.arange(3,6), np.arange(3,20,14)]
    # pks = [(0.5, 0.5), (0.5, 0.25, 0.25), (0.9, 0.1)]
    # max_times = [5, 4, 4]
    
    # ----Regular tree graph, diffusion-------
    # Regular tree sizes:
    # d=3, T=12 => E[N] = 155
    # d=4, T=8 => E[N] = 120
    # d=5, T=7 => E[N] = 157
    xks = [np.array([2,3]) for i in range(1)]
    pks = [(0.3,0.7) for i in range(1)]
    # pks = [(1.0) for i in range(1)] 
    
    # est_times: the timestamps at which to estimate the source
    # est_times = [6,8,10,12,14,16] # d=3
    # est_times = [16] # Adaptive diffusion irregular tree trials
    # est_times = [6,8,9,10,11] # d=4
    # est_times = [6,7,8,9] # d=5
    # est_times = [50,100,150,200] # d=2
    est_times = [9]
    
    max_times = [max(est_times) for i in range(1)] # the maximum time we run the algorithm
    
    max_infection = 100 #min(xks) - 1
    additional_time = 0  # collect additional_time more estimates after the first snapshot
        
    for (xk, pk, max_time) in zip(xks, pks, max_times):
        print('Checking xks = ',xk)
        # degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))

        # Make a random variable out of the support and probabilities listed
        degrees_rv2 = Multinomial(xk, pk)
        # max_infection = max_infection + 1
        # max_infection = max(xk) - 1
        
        print('Diffusion: ', diffusion)
       
        # Run regular diffusion
        if diffusion:
            num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection,
                                                                degrees_rv2, method=4, q=delay_parameter, 
                                                                spy_probability = spy_probability,
                                                                est_times = est_times)[:2]
            pd_ml, hop_distances, pd_spy, spy_hop_distances, pd_lei, lei_hop_distances = results
            print('hop distances',hop_distances)
            if write_results:
                prob = str(Decimal(spy_probability).quantize(Decimal('.01')))
                if isinstance(pk, float):
                    xk_str = str(xk[0])
                    pk_str = str(pk)
                    filename = 'results/spies/regular_trees/results_' + xk_str + "_" + pk_str + '_spies_'+prob + '_q_' + str(delay_parameter) + '_run_' + str(run) + '.mat'
                else:
                    xk_str = [str(i) for i in xk]
                    pk_str = [str(round(i,1)) for i in pk]
                    filename = 'results/spies/regular_trees/results_' + "_".join(xk_str) + "_" + "_".join(pk_str) + 'spies_'+ prob + '_run_' + str(run) + '.mat'
                io.savemat(filename,{'pd_ml':np.array(pd_ml), 'hop_distances':np.array(hop_distances),
                                     'pd_spy':np.array(pd_spy), 'spy_hop_distances':np.array(spy_hop_distances),
                                     'pd_lei':np.array(pd_lei), 'lei_hop_distances':np.array(lei_hop_distances),
                                     'num_infected':np.array(num_infected), 'est_times':np.array(est_times)})
                continue
        # Run adaptive diffusion (weighted spreading by neighbor degree)
        elif alt:
            num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv2, method = 1)[:2]
            hop_distances = []
        # Run adaptive diffusion (normal)
        else:
            if additional_time:
                num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv2, method = 0,
                                                                  additional_time = additional_time)[:2]
                print("additional pd: ", additional_pd)
            else:
                num_infected, results = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv2, method = 0,
                                                                    spy_probability=spy_probability,
                                                                    diffusion=diffusion,
                                                                    spies=spies,
                                                                    est_times = est_times)[:2]
            pd_ml, additional_pd, hop_distances = results
        print("pd ML: ", pd_ml)
        print('num_infected: ', num_infected)

        # write results, if desired
        if write_results:
            if spies:
                filename = 'results/spies/results_'
            else:
                filename = 'results/snapshot/results_'
            xk_str = [str(i) for i in xk]
            if isinstance(pk, (int, float, complex)):
                pk_str = str(round(pk,1))
                filename += xk_str[0] + "_" + pk_str
            else:
                pk_str = [str(round(i,1)) for i in pk]
                filename += "_".join(xk_str) + "_" + "_".join(pk_str)
            if spy_probability > 0.0:
                filename += 'spies' + str(Decimal(spy_probability).quantize(Decimal('.01')))
            if alt:
                filename += '_alt'
            filename += '_run_' + str(run) + '.mat'
            
            io.savemat(filename,{'pd_ml':np.array(pd_ml), 'hop_distances':np.array(hop_distances),
                                 'num_infected':np.array(num_infected), 'est_times':np.array(est_times)})
