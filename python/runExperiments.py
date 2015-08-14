# runExperiments
import utilities
import buildGraph
import infectionModels
import spyInfectionModels
import exportGraph
import adversaries
import random
import numpy as np

'''Runs a spreading algorithm over a real dataset. Either pramod's algo (deterministic) or message passing (ours)'''    
def run_dataset(filename, min_degree, trials, max_time=100, max_infection = -1, max_num_nodes = 4941,
                spies = False, spy_probability = 0.0, diffusion = False, q = 1.0):
    '''Run dataset runs a spreading algorithm over a dataset. 
    Inputs:
    
          filename:               name of the file containing the data
          min_degree:             remove all nodes with degree lower than min_degree
          trials:                 number of trials to run
          max_time(opt):          the maximum number of timesteps to use
          max_infection(opt):     the maximum number of nodes a node can infect in any given timestep
    
    Outputs:
    
          p:                      average proportion of nodes reached by the algorithm 
          num_infected:           total number of nodes infected by the algorithm
    
    NB: If max_infection is not set, then we'll run the deterministic algorihtm. Otherwise, it's the message passing one. '''
    
    print('spy probability',spy_probability)
    adjacency = buildGraph.buildDatasetGraph(filename, min_degree, max_num_nodes)
    num_nodes = len(adjacency)
    num_true_nodes = sum([len(item)>0 for item in adjacency])
    
    print('----Graph statistics:-----')
    degrees = [len(item) for item in adjacency]
    mean_degree = np.mean(degrees)
    print('the mean degree is',mean_degree)
    print('num true nodes',num_true_nodes)
    print('---------------------------\n')
    
    num_infected = 0
    p = 0
    pd_jordan = [0 for i in range(max_time)]
    pd_rumor = [0 for i in range(max_time)]
    pd_ml_leaf = [0 for i in range(max_time)]
    pd_spy = [0 for i in range(max_time)]
    avg_num_infected = [0 for i in range(max_time)]
    # avg_leaf_dists = [[0,0] for i in range(max_time)]
    avg_leaf_dists = [0 for i in range(max_time)]

    for trial in range(trials):
        # if trial % 20 == 0:
        print('Trial ',trial, ' / ',trials)
        while True:
            source = random.randint(0,num_nodes-1)
            if len(adjacency[source]) > 0:
                break
        if max_infection == -1:      # i.e. we're running the deterministic version
            num_infected, infection_pattern = infectionModels.infect_nodes_deterministic_reinfect(source,adjacency)
            print('num nodes',sum([1 for r in adjacency if len(r)>0]))
        else:
            if spies:
                # WIth spies
                
                if diffusion:  
                    # Run diffusion + first-spy estimator
                    diff_infector = spyInfectionModels.DatasetDiffusionInfector(adjacency, spy_probability, max_infection,q=q)
                    infection_details = diff_infector.infect(source, max_time)
                    who_infected, num_infected = infection_details
                    
                    print('infection done!')
                    # Estimate the source
                    adversary = adversaries.DiffusionSpiesAdversary(source, diff_infector.spies_info, who_infected)
                    results = adversary.get_estimates(max_time)
                    ml_leaf_correct, spy_correct, hop_distances = results
                    
                else:
                    # Run adaptive diffusion + spies ML estimator
                
                    # Infect nodes
                    up_down_infector = spyInfectionModels.DatasetUpDownInfector(adjacency, spy_probability, max_infection)
                    infection_details = up_down_infector.infect(source, max_time)
                    who_infected, num_infected = infection_details
                    
                    # Estimate the source
                    adversary = adversaries.DatasetUpDownAdversary(source, up_down_infector.spies_info, who_infected, adjacency, max_infection)
                    results = adversary.get_estimates(max_time)
                    ml_leaf_correct, hop_distances = results
                    

            else:   # snapshot adversary
                num_infected, infection_pattern, who_infected, results = infectionModels.infect_nodes_adaptive_diff(source,adjacency,max_time,max_infection)
                # unpack the results
                jordan_correct, rumor_correct, ml_leaf_correct, ml_leaf_dists = results
                print('dists', ml_leaf_dists)
            
                pd_jordan = [i+j for (i,j) in zip(pd_jordan, jordan_correct)]
                pd_rumor = [i+j for (i,j) in zip(pd_rumor, rumor_correct)]
                # avg_leaf_dists = [[k+m for (k,m) in zip(i,j)] for (i,j) in zip(ml_leaf_dists, avg_leaf_dists)]
                avg_leaf_dists = [i+j for (i,j) in zip(ml_leaf_dists, avg_leaf_dists)]
            
            pd_ml_leaf = [i+j for (i,j) in zip(pd_ml_leaf, ml_leaf_correct)]
            if spies and diffusion:
                pd_spy = [i+j for (i,j) in zip(pd_spy, spy_correct)]
            # write the infected subgraph to file
            # filename = 'infected_subgraph_'+str(trial)
            # exportGraph.export_gexf(filename,who_infected,source,infection_pattern,adjacency)
        print("infected", num_infected)
        p += float(num_infected[-1]) / num_true_nodes
        avg_num_infected = [i+j for (i,j) in zip(avg_num_infected, num_infected)]
    p = p / trials
    pd_jordan = [float(i)/trials for i in pd_jordan]
    pd_rumor = [float(i)/trials for i in pd_rumor]
    pd_ml_leaf = [float(i) / trials for i in pd_ml_leaf]
    if spies and diffusion:
        pd_spy = [float(i) / trials for i in pd_spy]
    avg_num_infected = [ float(i) / trials for i in avg_num_infected]
    # avg_leaf_dists = [[float(k) / trials for k in i] for i in avg_leaf_dists]
    avg_leaf_dists = [float(i) / trials for i in avg_leaf_dists]
    if spies:
        if diffusion:
            results = (pd_ml_leaf, pd_spy, avg_leaf_dists)
        else:
            results = (pd_ml_leaf, avg_leaf_dists)
    else:
        results = (pd_jordan, pd_rumor, pd_ml_leaf, avg_leaf_dists)
    return p, avg_num_infected, results
    
'''Run a random tree'''    
def run_randtree(trials, max_time, max_infection, degrees_rv, method=0, known_degrees=[], additional_time = 0,
                 q = 0.5, spies=False, spy_probability = 0.0, diffusion=False, est_times=None):
    ''' Run dataset runs a spreading algorithm over a dataset. 
    
    Arguments:    
        trials:                 number of trials to run
        max_time:               the maximum number of timesteps to use
        max_infection:          max number of nodes an infected node can infect in 1 timestep
        degrees_rv:             random tree degree random variable
        method:                 which approach to use: 0 = regular spreading to all neighbors, 
                                                       1 = alternative spreading (VS chosen proportional to degree of candidates) 
                                                       2 = pre-planned spreading
                                                       3 = Old adaptive diffusion (i.e. line algorithm)
                                                       4 = regular diffusion (symmetric in all directions)
        known_degrees           the degrees of nodes in the tree for likelihood computation
        additional_time         the number of additional estimates to make after the initial snapshot
        p                       probability with which diffusion passes the message in a timestep
        spy_probability         probability of a node being a spy
        est_times               times at which we estimate the source (must be <= max_time)
    
    Outputs:
    
          p:                      delay rate 
          avg_num_infected:           total number of nodes infected by the algorithm
    
    NB: If max_infection is not set, then we'll run the deterministic algorihtm. Otherwise, it's the message passing one.'''
    
    pd_ml = [0 for i in range(max_time)] # Maximum likelihood
    pd_spy = [0 for i in range(max_time)] # Maximum likelihood
    pd_lei = [0 for i in range(max_time)] # Maximum likelihood
    additional_pd_mean = [0 for i in range(additional_time)] # Pd from additional measurements
    pd_rc = [0 for i in range(max_time)] # Rumor centrality
    pd_jc = [0 for i in range(max_time)] # Jordan centrality
    pd_rand_leaf = [0 for i in range(max_time)]
    avg_num_infected = [0 for i in range(max_time)]
    avg_hop_distance = [0 for i in range(max_time)]
    avg_spy_hop_distance = [0 for i in range(max_time)]
    avg_lei_hop_distance = [0 for i in range(max_time)]
    
    for trial in range(trials):
        if trial % 10 == 0:
            print('\nTrial ',trial, ' / ',trials-1)
        source = 0
        if method == 0:      # Infect nodes with adaptive diffusion over an irregular tree (possibly with multiple snapshots)
            if spies:
                # WIth spies
                
                if not diffusion:   # adaptive diffusion
                    # Infect nodes
                    up_down_infector = spyInfectionModels.UpDownInfector(spy_probability,degrees_rv)
                    infection_details = up_down_infector.infect(source, max_time)
                    who_infected, num_infected = infection_details
                    
                    # Estimate the source
                    adversary = adversaries.UpDownAdversary(source, up_down_infector.spies_info, who_infected, degrees_rv)
                    results = adversary.get_estimates(max_time, est_times)
                    ml_correct, hop_distances = results
                    additional_pd = [0 for i in range(additional_time)]
            else:
                # snapshot, adaptive diffusion
                infection_details, ml_correct, additional_pd = infectionModels.infect_nodes_adaptive_irregular_tree(source, max_time, max_infection,
                                                                                                 degrees_rv, additional_time = additional_time, alt=False)
                num_infected, infection_pattern, who_infected, additional_hops = infection_details
                hop_distances = []
            # print(additional_hops)
            additional_pd_mean = [i+j for (i,j) in zip(additional_pd, additional_pd_mean)]
        elif method == 1:    # Infect nodes with adaptive diffusion over an irregular tree, alternative spreading
            infection_details, ml_correct = infectionModels.infect_nodes_adaptive_irregular_tree(source, max_time, max_infection, 
                                                                                                 degrees_rv, True, spy_probability=spy_probability)
            num_infected, infection_pattern, who_infected = infection_details
        elif method == 2:    # Infect nodes with adaptive diffusion over a pre-determined irregular tree
            infection_details, ml_correct, rand_leaf_correct, known_degrees = infectionModels.infect_nodes_adaptive_planned_irregular_tree(source, max_time, max_infection, degrees_rv, known_degrees)
            num_infected, infection_pattern, who_infected = infection_details
            pd_rand_leaf = [i+j for (i,j) in zip(pd_rand_leaf, rand_leaf_correct)]
        elif method == 3:   # Infect nodes with adaptive diffusion over a line
            infection_details, results = infectionModels.infect_nodes_line_adaptive(source, max_time, max_infection, degrees_rv)
            num_infected, who_infected = infection_details
            pd_jc = [i+j for (i,j) in zip(pd_jc, results[0])]
            pd_rc = [i+j for (i,j) in zip(pd_rc, results[1])]
            # We don't actually compute the ML estimate here because it's computationally challenging
            ml_correct = pd_ml 
        elif method == 4:   # infect nodes with regular diffusion
            infection_details, results = infectionModels.infect_nodes_diffusion_irregular_tree(source, max_time, degrees_rv, q, spy_probability,
                                                                                               est_times = est_times, diffusion=True)
            ml_correct, spy_correct, lei_correct = results
            num_infected, who_infected, hop_distances, spy_hop_distances, lei_hop_distances = infection_details
            avg_hop_distance  = [i+j for (i,j) in zip(avg_hop_distance, hop_distances)]
            avg_spy_hop_distance  = [i+j for (i,j) in zip(avg_spy_hop_distance, spy_hop_distances)]
            avg_lei_hop_distance  = [i+j for (i,j) in zip(avg_lei_hop_distance, lei_hop_distances)]
            pd_spy = [i+j for (i,j) in zip(pd_spy, spy_correct)]
            pd_lei = [i+j for (i,j) in zip(pd_lei, lei_correct)]
        # Infect nodes with adaptive diffusion over an irregular tree, alternative spreading
        # unpack the results
        pd_ml = [i+j for (i,j) in zip(pd_ml, ml_correct)]
        avg_num_infected = [i+j for (i,j) in zip(avg_num_infected, num_infected)]
        
    pd_ml = [float(i) / trials for i in pd_ml]
    pd_rand_leaf = [float(i) / trials for i in pd_rand_leaf]
    avg_num_infected = [ float(i) / trials for i in avg_num_infected]
    
    if method == 0:
        additional_pd_mean = [i/trials for i in additional_pd_mean]
        results = (pd_ml, additional_pd_mean, hop_distances)
    elif method == 1:
        results = (pd_ml)
    elif method == 2:
        results = (pd_ml, pd_rand_leaf)
    elif method == 3:
        pd_rc = [float(i) / trials for i in pd_rc]
        pd_jc = [float(i) / trials for i in pd_jc]
        results = (pd_rc, pd_jc)
    elif method == 4:
        pd_spy = [float(i) / trials for i in pd_spy]
        avg_spy_hop_distance = [float(i) / trials for i in avg_spy_hop_distance]
        pd_lei = [float(i) / trials for i in pd_lei]
        avg_lei_hop_distance = [float(i) / trials for i in avg_lei_hop_distance]
        avg_hop_distance = [float(i) / trials for i in avg_hop_distance]
        print('pd_ml: ', pd_ml, 'avg_hop_distance', avg_hop_distance)
        print('pd_spy: ', pd_spy, 'avg_spy_hop_distance', avg_spy_hop_distance)
        print('pd_lei: ', pd_lei, 'avg_lei_hop_distance', avg_lei_hop_distance)
        results = (pd_ml, avg_hop_distance, pd_spy, avg_spy_hop_distance, pd_lei, avg_lei_hop_distance)

    # return avg_num_infected, pd_ml, pd_rand_leaf
    return avg_num_infected, results
