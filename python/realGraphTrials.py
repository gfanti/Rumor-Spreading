import buildGraph
import random
import scipy.io as io
import numpy as np
import utilities, runExperiments, infectionModels
import sys, time
from sys import platform as _platform
  
if __name__ == "__main__":

    # Parse the arguments
    trials, write_results, database, run = utilities.parse_args(sys.argv)

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
    max_infection = 3
    
    start = time.clock()
    p_fb, num_infected, results = runExperiments.runDataset(filename, min_degree, trials, max_time, max_infection, 10000)
    print('Experiment took ', time.clock() - start, ' seconds.')
    pd_jordan, pd_rumor, pd_ml_leaf, ml_leaf_dists = results
    print('Facebook result: ',p_fb,num_infected)
    print('Accuracy using Jordan centrality: ',pd_jordan)
    print('Accuracy using rumor centrality: ',pd_rumor)
    print('Accuracy using ML leaf: ',pd_ml_leaf)
    print('Average distances of the ML estimate: ', ml_leaf_dists)
    
    
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
        filename += '_do3_' + str(run)
        io.savemat(filename,{'pd_jordan':np.array(pd_jordan),'pd_rumor':np.array(pd_rumor), 'pd_ml_leaf':np.array(pd_ml_leaf), 'time':np.array([i for i in range(max_time)]), 'num_infected':np.array(num_infected), 'ml_leaf_dists':np.array(ml_leaf_dists)})