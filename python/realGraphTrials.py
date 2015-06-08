import buildGraph
import random
import scipy.io as io
import numpy as np
import utilities, runExperiments, infectionModels
import sys, time
from sys import platform as _platform
  
if __name__ == "__main__":

    # Parse the arguments
    # trials, write_results, database, run = utilities.parse_args(sys.argv)
    args = utilities.parse_args(sys.argv)
    trials = args.get('trials', 1)
    write_results = args.get('write_results', False)
    database = args.get('database', 'fb')
    spy_probability = args.get('spy_probability', 0.0)
    diffusion = args.get('diffusion', False)
    spies = args.get('spies', False)
    q = args.get('delay_parameter', 1.0)
    run = args.get('run', 1)

    
    # Dataset options:
    # 1: Facebook
    # 2: Power Grid
    if _platform == "linux" or _platform == "linux2":
        # linux
        prefix = '../data/'
    elif _platform == "win32":
        prefix = '..\\data\\'

        
    ##---------- Real Graphs -------------------#
    
    # facebook
    if database == 'fb':
        filename = prefix + 'out.facebook-wosn-links'
    else:
        filename = prefix + 'out.opsahl-powergrid'
    
    min_degree = 3;
    max_time = 10
    max_infection = 20
    max_graph_size = 10000
    
    start = time.clock()
    p_fb, num_infected, results = runExperiments.run_dataset(filename, min_degree, trials, max_time, max_infection, max_graph_size,
                                                             spies=spies, spy_probability=spy_probability,
                                                             diffusion=diffusion, q=q)
    print('Experiment took ', time.clock() - start, ' seconds.')
    if spies:
        if diffusion:
            pd_ml_leaf, pd_spy, ml_leaf_dists = results
            print('Accuracy using first-spy: ',pd_spy)
        else:
            pd_ml_leaf, ml_leaf_dists = results
            
    else:
        pd_jordan, pd_rumor, pd_ml_leaf, ml_leaf_dists = results
        print('Accuracy using Jordan centrality: ',pd_jordan)
        print('Accuracy using rumor centrality: ',pd_rumor)
    
    print('Accuracy using ML leaf: ',pd_ml_leaf)
        
    print('Facebook result: ',p_fb,', Avg. num infected: ', num_infected)
    
    # print('Average distances of the ML estimate: ', ml_leaf_dists)
    
    
    if write_results:
        if _platform == "linux" or _platform == "linux2":
            # linux
            prefix = 'results/'
        elif _platform == "win32":
            # windows
            prefix = 'results\\'
        if database == 'fb':
            filename = prefix + 'pd_facebook'
        else:
            filename = prefix + 'pd_power_grid'
        filename += '_' + str(run)
        io.savemat(filename,{'pd_jordan':np.array(pd_jordan),'pd_rumor':np.array(pd_rumor), 'pd_ml_leaf':np.array(pd_ml_leaf), 'time':np.array([i for i in range(max_time)]), 'num_infected':np.array(num_infected), 'ml_leaf_dists':np.array(ml_leaf_dists)})