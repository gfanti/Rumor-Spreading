import buildGraph
import infectionModels
import random
import scipy.io as io
import numpy as np
import utilities
  
if __name__ == "__main__":

    trials = 2  
    prefix = '..\\data\\'

        
    ##---------- Real Graphs -------------------#
    
    # facebook
    filename = prefix + 'out.facebook-wosn-links'
    min_degree = 3;
    max_time = 10
    max_infection = 3
    p_fb, num_infected, pd_jordan, pd_rumor, pd_ml = utilities.runDataset(filename, min_degree, trials, max_time, max_infection)
    print('Facebook result: ',p_fb,num_infected)
    print('Accuracy using Jordan centrality: ',pd_jordan)
    print('Accuracy using rumor centrality: ',pd_rumor)
    print('Accuracy using ML: ',pd_ml)
    
    io.savemat('pd',{'pd_jordan':np.array(pd_jordan),'pd_rumor':np.array(pd_rumor), 'time':np.array([i for i in range(max_time)])})
    exit()
    
    # power grid
    filename = prefix + 'out.opsahl-powergrid'
    min_degree = 5;
    p_powergrid, num_infected = runDataset(filename, min_degree, trials)
    print('power grid result: ',p_powergrid,num_infected)
    exit() 
    
    
    
    
    
    # ## 3-regular tree
    # adjacency = [[]]
    # max_degree = 3
    # max_time = 5
    # alpha = 0.5;
    # beta = 1;
    # adjacency, num_infected = infectionModels.infect_nodes_adaptive_diff_tree(0, adjacency, max_degree, max_time, alpha, beta)
    # print('num infected ',num_infected)
    # exit()