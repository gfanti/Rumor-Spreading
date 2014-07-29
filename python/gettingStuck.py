import buildGraph
import infectionModels
import random
import scipy.io as io
import numpy as np


'''Runs a spreading algorithm over a real dataset. Either pramod's algo (deterministic) or message passing (ours)'''    
def runDataset(filename, min_degree, trials, max_time=100, max_infection = -1):
    # Run dataset runs a spreading algorithm over a dataset. 
    # Inputs:
    #
    #       filename:               name of the file containing the data
    #       min_degree:             remove all nodes with degree lower than min_degree
    #       trials:                 number of trials to run
    #       max_time(opt):          the maximum number of timesteps to use
    #       max_infection(opt):     the maximum number of nodes a node can infect in any given timestep
    #
    # Outputs:
    #
    #       p:                      average proportion of nodes reached by the algorithm 
    #       num_infected:           total number of nodes infected by the algorithm
    #
    # NB: If max_infection is not set, then we'll run the deterministic algorihtm. Otherwise, it's the message passing one.
    
    adjacency = buildGraph.buildDatasetGraph(filename, min_degree)
    # print(adjacency)
    num_nodes = len(adjacency)
    num_true_nodes = sum([len(item)>0 for item in adjacency])
    print('nonzero nodes:',num_true_nodes)
    num_infected = 0
    p = 0
    
    for trial in range(trials):
        if trial % 20 == 0:
            print('Trial ',trial, ' / ',trials)
        while True:
            source = random.randint(0,num_nodes-1)
            if len(adjacency[source]) > 0:
                break
        if max_infection == -1:      # i.e. we're running the deterministic version
            num_infected, infection_pattern = infectionModels.infect_nodes_deterministic(source,adjacency)
        else:
            num_infected, infection_pattern = infectionModels.infect_nodes_adaptive_diff(source,adjacency,max_time,max_infection)
        p += num_infected / (1.0 * num_true_nodes)
    p = p / trials
    return p, num_infected
  
if __name__ == "__main__":

    trials = 100
    prefix = '..\\data\\'

        
    ##---------- Real Graphs -------------------#
    
    # facebook
    filename = prefix + 'out.facebook-wosn-links'
    min_degree = 3;
    max_time = 10
    alpha = 0.5;
    beta = 1;
    max_infection = 3
    p_fb, num_infected = runDataset(filename, min_degree, trials, max_time, max_infection)
    print('Facebook result: ',p_fb,num_infected)
   
    exit()
    
    # power grid
    filename = prefix + 'out.opsahl-powergrid'
    min_degree = 5;
    p_powergrid, num_infected = runDataset(filename, min_degree, trials)
    print('power grid result: ',p_powergrid,num_infected)
    exit() 
    
    
    
    ##---------- Synthetic Graphs ---------------#
    num_nodes_range = range(10,661,50)

    p_reached_smallworld = []
    p_reached_scalefree = []
    
    for num_nodes in num_nodes_range:
        p_smallworld = 0
        p_scalefree = 0
        print('Nodes: ',num_nodes)
        for trial in range(trials):
            #------ Build graphs---------#
            # build a barabasi-albert (scale-free) graph
            adjacency_scalefree = buildGraph.buildBarabasiAlbertGraph(num_nodes);
            
            # build a watts-strogatz graph
            desired_degree = int(sum([len(neighbors) for neighbors in adjacency_scalefree]) / 2 / num_nodes)
            desired_degree = min(int((num_nodes-1)/2),desired_degree) # make sure we don't have too many connections
            # print('degree: ',desired_degree)
            adjacency_smallworld = buildGraph.buildSmallWorldGraph(num_nodes, desired_degree);
           
            #--------- Spread the message on both graphs ------------#
            source = random.randint(0,num_nodes-1)
            num_infected, infection_pattern = infectionModels.infect_nodes_deterministic(source,adjacency_scalefree)
            p_scalefree += num_infected / (1.0 * num_nodes)
            
            num_infected, infection_pattern = infectionModels.infect_nodes_deterministic(source,adjacency_smallworld)
            p_smallworld += num_infected / (1.0 * num_nodes)
        
        p_scalefree = p_scalefree / trials
        p_smallworld = p_smallworld / trials
        
        p_reached_scalefree = p_reached_scalefree + [p_scalefree]
        p_reached_smallworld = p_reached_smallworld + [p_smallworld]
    print("Proportion of nodes reached (small-world): ", p_reached_smallworld)
    print("Proportion of nodes reached (scale-free): ", p_reached_scalefree)
    
    n = [i for i in num_nodes_range]
    io.savemat('synthetic_results',{'p_scalefree':np.array(p_reached_scalefree), 'p_smallworld':np.array(p_reached_smallworld), 'n':n})
    
    
    # ## 3-regular tree
    # adjacency = [[]]
    # max_degree = 3
    # max_time = 5
    # alpha = 0.5;
    # beta = 1;
    # adjacency, num_infected = infectionModels.infect_nodes_adaptive_diff_tree(0, adjacency, max_degree, max_time, alpha, beta)
    # print('num infected ',num_infected)
    # exit()