import buildGraph
import infectionModels
import random
import scipy.io as io
import numpy as np
import utilities
import runExperiments
import sys
from sys import platform as _platform
  
if __name__ == "__main__":

    # Parse the arguments
    trials, write_results, database = utilities.parse_args(sys.argv)

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
    p_fb, num_infected, pd_jordan, pd_rumor, pd_ml = runExperiments.runDataset(filename, min_degree, trials, max_time, max_infection, 10000)
    print('Facebook result: ',p_fb,num_infected)
    print('Accuracy using Jordan centrality: ',pd_jordan)
    print('Accuracy using rumor centrality: ',pd_rumor)
    print('Accuracy using ML: ',pd_ml)
    
    if write_results:
        if dataset == 'fb':
            filename = 'pd_facebook'
        else:
            filename = 'pd_power_grid'
        io.savemat('pd',{'pd_jordan':np.array(pd_jordan),'pd_rumor':np.array(pd_rumor), 'time':np.array([i for i in range(max_time)])})