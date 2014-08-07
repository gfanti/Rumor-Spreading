import buildGraph
import infectionModels
import random
import scipy.io as io
import numpy as np
    
  
if __name__ == "__main__":

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
