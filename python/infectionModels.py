# node infection
import estimation
import estimation_spies
import random
from scipy import stats
import numpy as np
import utilities
import networkx as nx
import time
# import cProfile

#Semi-distributed adaptive diffusion over regular trees (uses two timesteps)
def infect_nodes_adaptive_tree(source, adjacency, max_degree, max_time, alpha):
    '''infect_nodes_adaptive_irregular_tree runs the spreading model over a regular
    tree.
    
    Arguments:
        source: The ID of the source node (usually 0)
        adjacency: The adjacency relations of the underlying tree
        max_degree: The maximum number of nodes an infected node will infect in any timestep
        max_time: The maximum number of timesteps to run the algorithm
        alpha: Probability with which the virtual source moves
      
    Returns:
        adjacency: (updated)
        num_infected: Number of infected nodes
    '''
    num_nodes = len(adjacency)
    timesteps = 1;
    
    # initially the virtual source and the true source are the same
    virtual_source = source
    virtual_source_candidate = virtual_source
    
    blocked = False
    
    while timesteps < max_time:
        print('time', timesteps)
    
        # in odd timesteps, choose a direction to expand in
        if timesteps == 1:
            adjacency = utilities.infect_nodes_infinite_tree(source, 1, adjacency)
            virtual_source_candidate = 1
            timesteps += 1
            continue
            
        if timesteps % 2 == 1:
            current_neighbors = [k for k in adjacency[virtual_source]]
            # If blocked, break out of the loop
            if len(current_neighbors) < 1:
                blocked = True
                break
                
            virtual_source_candidate, current_neighbors = pick_random_elements(current_neighbors,1)
            virtual_source_candidate = virtual_source_candidate[0]
            adjacency = pass_branch_message_infinite_tree(virtual_source, virtual_source_candidate, adjacency, max_degree)
            
            
        # in even timesteps, choose whether to move the virtual source
        else:
            # With probability alpha, move the virtual source
            if random.random() < utilities.compute_alpha(m,timesteps,max_infection):
                adjacency = pass_branch_message_infinite_tree(virtual_source, virtual_source_candidate, adjacency, max_degree)
                virtual_source = virtual_source_candidate
            # Otherwise, keep it where it is, and spread symmetrically
            else:
                adjacency = pass_branch_message_infinite_tree(virtual_source_candidate, virtual_source, adjacency, max_degree)
    
        num_infected = len(adjacency)
        print('num infected', num_infected)
        timesteps += 1
        
    return adjacency, num_infected

#Distributed diffusion over random trees
def infect_nodes_diffusion_irregular_tree(source, max_time, degrees_rv, p = 0.5,
                                          spy_probability = 0.0, est_times = None):
    '''infect_nodes_diffusion_tree runs the spreading model over a regular
    tree.
    
    Arguments:
        source: The ID of the source node (0)
        degrees_rv: The random variable describing the degree distribution
        p: probability of passing message
        max_time: The maximum number of timesteps to run the algorithm
      
    Returns:
        adjacency: (updated)
        timestamps: the times at which each node got infected
        num_infected: Number of infected nodes
    '''
    timesteps = 1;
    
    ml_correct = []
    spy_correct = []
    hop_distances = []
    spy_hop_distances = []
    tot_num_infected = []
    num_infected = 0
    boundary_nodes = [source]
    who_infected = [[]]
    degrees = degrees_rv.rvs(size=1).tolist()
    timestamps = [0]
    active_nodes = [0] # marks which nodes are candidate sources >=0 => valid, <0 => not valid
    spies = []
    active_spies = []
    while boundary_nodes:
    # while timesteps <= max_time:
    # while (1 in [active_nodes[item] for item in boundary_nodes]) and timesteps < 10:
        # print('time', timesteps)
    
        num_candidates = len(boundary_nodes)
        node = boundary_nodes.pop(0)
        # print('degrees[node]', node, len(degrees), len(who_infected))
        num_uninfected_neighbors = degrees[node] - len(who_infected[node])
        # num_to_infect = np.random.binomial(num_uninfected_neighbors, p)
        num_to_infect = num_uninfected_neighbors
        # print('num to infect', num_to_infect, 'out of ', num_uninfected_neighbors,'\n')
        # print('who_infected before: ', who_infected)
        if num_to_infect > 0:
            # timestamps += [timesteps for j in range(num_to_infect)]
            # timestamps += [j + timestamps[node] for j in np.random.exponential(p, num_to_infect)]
            # neighbor_times = [max(j,0) + timestamps[node] for j in np.random.normal(1.0/p, 0.5, num_to_infect)]
            neighbor_times = [max(j,0) + timestamps[node] for j in np.random.normal(1.0/p, 1.7321, num_to_infect)] # Truncated gaussian
            # neighbor_times = [max(j,0) + timestamps[node] for j in np.random.geometric(p, num_to_infect)]
            # neighbor_times = [max(j,0) + timestamps[node] for j in np.random.normal(20, 5, num_to_infect)]
            # neighbors = [neighbors[k] for k in range(num_to_infect) if neighbor_times[k] <= max_time]
            neighbor_times = [neighbor_times[k] for k in range(num_to_infect) if neighbor_times[k] <= max_time]
            neighbors = ([k+len(degrees) for k in range(len(neighbor_times))])
            
            timestamps += neighbor_times
            
            degrees, who_infected = utilities.infect_nodes_randtree(node, neighbors, degrees, degrees_rv, who_infected)[:2]
            new_spies = utilities.update_spies_diffusion(neighbors, spy_probability=spy_probability)
            # new_spies = [spy for spy in new_spies if timestamps[spy] <= max_time]
            spies += new_spies
            # mark whether the new additions are possible candidates
            if active_nodes[node] < 0:
                active_nodes += [active_nodes[node] for item in neighbor_times]
            else:
                active_nodes += neighbor_times
                for (neighbor,neighbor_time) in zip(neighbors,neighbor_times):
                    if neighbor in new_spies:
                        active_nodes[neighbor] = -neighbor_time
                        active_spies += [neighbor]
            
            # boundary_nodes += [neighbor for neighbor in neighbors if timestamps[neighbor] <= max_time]
            boundary_nodes += neighbors
            if num_to_infect < num_uninfected_neighbors:
                boundary_nodes.append(node)
        else:
            boundary_nodes.append(node)
        # print('who_infected after: ', who_infected)
                
        
        num_infected = len(who_infected)
        # print('infected',num_infected)
        # print('who infected',who_infected)
        # print('bondarey nodes',boundary_nodes)
        # CHANGE THIS!
        # print('timestamps: ', timestamps, len(timestamps))
        # print('spies are', spies)
        # timesteps += 1
        
    print('num_infected: ', num_infected)
    
    
        
    # Remove the nodes with time greater than current timestamp
        
    adjacency = [set(item) for item in who_infected]
    nodes = [i for i in range(num_infected)]

    # for node in nodes:
        # print('Node ',node,': ',timestamps[node])
    for est_time in est_times:
        
        reached_spies = [spy for spy in active_spies if timestamps[spy] <= est_time]
        print('reached ',len(reached_spies), ' spies')
        
            
        spies_timestamps = [timestamps[j] for j in active_spies if timestamps[j] <= est_time]
        
        # current_active_nodes = [item if ((abs(item) <= est_time) and (n not in active_spies)) else abs(item) for (item,n) in zip(active_nodes,nodes)]
        if len(reached_spies) == 0:
            current_active_nodes = [1 if n not in spies else -1 for n in nodes]
        else:
            current_active_nodes = active_nodes
        print('est time',est_time)
        # print('original active_nodes',active_nodes)
        # print('new_active_nodes', current_active_nodes)
        # spies_timestamps = [timestamps[j] for j in spies]
        # print('spies:',spies)
        # print('active spies:', active_spies)
        # print('active spies timestamps', spies_timestamps)
        # print('adjacnecy', [adjacency[i] for i in active_spies])
        # print('active nodes: ', [item for item in range(len(active_nodes)) if active_nodes[item]==1])
        
        # spy_active_nodes = [1 if i not in spies else -1 for i in range(num_infected)]
        # print('Considering ',len(reached_spies), ' spies: ',reached_spies)
  
  
        # estimator = estimation_spies.OptimalEstimator(adjacency, spies, spies_timestamps, active_nodes=spy_active_nodes)
        # spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, spies, spies_timestamps, active_nodes=spy_active_nodes)
        # estimator = estimation_spies.OptimalEstimator(adjacency, spies, spies_timestamps, active_nodes=active_nodes)
        # spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, spies, spies_timestamps, active_nodes=active_nodes)
        # estimator = estimation_spies.OptimalEstimator(adjacency, active_spies, spies_timestamps, active_nodes=current_active_nodes)
        # spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, active_spies, spies_timestamps, active_nodes=current_active_nodes)
        estimator = estimation_spies.OptimalEstimator(adjacency, reached_spies, spies_timestamps, active_nodes=current_active_nodes)
        spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, reached_spies, spies_timestamps, active_nodes=current_active_nodes)
        old_spy_estimate = -1
        
        # for spy in reached_spies:
            # print('spy',spy,' has time ',timestamps[spy],' and is ',nx.shortest_path_length(estimator.graph,0,spy),' hops away')
        
        if active_spies:
            start = time.time()
            ml_estimate = estimator.estimate_source()
            # print('Surrounded? ',not any([current_active_nodes[n]>0 for n in estimator.graph if len(estimator.adjacency[n])==1]))
            print('ml est', ml_estimate)
            # print('spies:        ',reached_spies)
            # print('all spies:    ',spies)
            # print('active spies: ',active_spies)
            # print('timestamps',spies_timestamps)
            # estimator.draw_graph()
            
            spy_estimate = spy_estimator.estimate_source()
            if spy_estimate != old_spy_estimate and old_spy_estimate != -1:
                print('CHANGED THE ESTIMATE!')
                exit(0)
            old_spy_estimate = spy_estimate
            end = time.time()
            print('elapsed:',end-start)
        else:
            # choose a random node
            ml_estimate = random.randint(0,num_infected - 1)
            spy_estimate = ml_estimate
        hop_distance = nx.shortest_path_length(estimator.graph, source, ml_estimate)
        spy_hop_distance = nx.shortest_path_length(estimator.graph, source, spy_estimate)

        print('True source: ', source, ' estimate: ', ml_estimate)
        ml_correct.append(ml_estimate == source)
        hop_distances.append(hop_distance)
        spy_correct.append(spy_estimate == source)
        spy_hop_distances.append(spy_hop_distance)
        tot_num_infected.append(num_infected)
    # max_time = max(spies_timestamps)
    # # print('spies_timestamps:', spies_timestamps)
    # print('spy timestamps', spies_timestamps)
    # for t_o in range(int(max_time)):
        # t_o_timestamps = [spies_timestamps[i] for i in range(len(active_spies)) if spies_timestamps[i] < t_o]
        # t_o_spies = [active_spies[i] for i in range(len(active_spies)) if spies_timestamps[i] < t_o]
        # if t_o_spies:
            # # print('spy timestamps', spies_timestamps)
            # estimator = estimation_spies.OptimalEstimator(adjacency, t_o_spies, t_o_timestamps, active_nodes)
            # ml_estimate = estimator.estimate_source()
        # else:
            # ml_estimate = -1
        # print('True source is ', source)
        # ml_correct.append(ml_estimate == source)
        # tot_num_infected.append(num_infected)
        
    print('Num spies are: ', len(spies), ' out of ', num_infected)
    print('\n\n')
    infection_details = (tot_num_infected, who_infected, hop_distances, spy_hop_distances)
    results = (ml_correct, spy_correct)
    return infection_details, results
    
# Semi-distributed adaptive diffusion over irregular trees (uses 1 timestep, this
# is the version we used in our simulations)
def infect_nodes_adaptive_irregular_tree(source, max_time, max_infection,
                                              degrees_rv, alt = False,
                                              additional_time = 0):
    '''infect_nodes_adaptive_irregular_tree runs the spreading model over an irregular
    tree.
    
    Arguments:
        source: The ID of the source node (usually 0)
        max_time: The maximum amount of timesteps to run the algorithm
        max_infection: Defines the (nominal degree-1) that is assumed in order to choose alpha
        degrees_rv: A random variable that describes the graph degree
        alt: Tells whether you would like to use the alternative spreading (weights virtual
            sources by their degree)
      
    Returns:
        infection_details: A list of the characteristics of the infection:
            - tot_num_infected: Total number of infected nodes at each timestep
            - infection_pattern: Which nodes got infected in which order
            - who_infected: Adjacency matrix of infected subgraph
        ml_correct: The ML estimate of the true source'''
    
    timesteps = 0
    
    # initially the virtual source and the true source are the same
    virtual_source = source

    # ML estimate
    ml_correct = [0 for i in range(max_time)]
    tot_num_infected = [0 for i in range(max_time)]
    num_infected = 0
    
    who_infected = [[]]
    degrees = degrees_rv.rvs(size=1).tolist()
    while timesteps < max_time:
    
            
        if timesteps == 0:
            virtual_source = 1
            previous_vs = 0
            
            # infect twice in one direction, always
            degrees, who_infected = utilities.infect_nodes_randtree(source, [virtual_source], degrees, degrees_rv, who_infected)[:2]
            infection_pattern, who_infected = pass_branch_message_randtree(source, virtual_source, degrees, degrees_rv, who_infected)[:2]
            m = 1       # the virtual source is always going to be 1 hop away from the true source
            
        else:
            current_neighbors = [k for k in who_infected[virtual_source]]

            # with probability alpha, spread symmetrically (keep the virtual source where it is)
            if random.random() < utilities.compute_alpha(m,timesteps,max_infection):
                # branch once in every direction
                for neighbor in current_neighbors:
                    degrees, who_infected = pass_branch_message_randtree(virtual_source, neighbor, degrees, degrees_rv, who_infected)[:2]
            
            # Otherwise, spread symmetrically
            else:
                # find a direction to move
                if alt: # Weight new virtual sources by their degree
                    current_neighbors.remove(previous_vs)
                    weights = [degrees[i] for i in current_neighbors]
                    # print('degree', degrees[virtual_source])
                    # print(weights)
                    weights = [i/sum(weights) for i in weights]
                    # print(weights)
                    virtual_source_rv = stats.rv_discrete(name='rv_discrete', values=(current_neighbors, weights))
                    virtual_source_candidate = virtual_source_rv.rvs(size=1)
                    # print('alt: ', timesteps, degrees[virtual_source_candidate])
                else: # Choose uniformly between all virtual sources
                    virtual_source_candidate = [previous_vs]
                    while virtual_source_candidate[0] == previous_vs:
                        virtual_source_candidate, current_neighbors, new_vs_likelihood = pick_random_elements(current_neighbors,1)
                    virtual_source_candidate = virtual_source_candidate[0]
                    # print('regular: ', timesteps, degrees[virtual_source_candidate])
                previous_vs = virtual_source
                
                # the virtual source moves one more hop away from the true source
                m += 1;
                
                # branch twice in one direction
                degrees, who_infected = pass_branch_message_randtree(virtual_source, virtual_source_candidate, degrees, degrees_rv, who_infected)[:2]
                degrees, who_infected = pass_branch_message_randtree(virtual_source, virtual_source_candidate, degrees, degrees_rv, who_infected)[:2]
                
                virtual_source = virtual_source_candidate
            
        num_infected = len(who_infected)

        # Estimating the error
        # ML estimate
        if alt: # Call the weighted estimator
            ml_estimate = estimation.ml_estimate_irregular_trees(max_infection, max_time, virtual_source, degrees, who_infected, degrees_rv, 1)
        else: # Call the regular estimator
            ml_estimate = estimation.ml_estimate_irregular_trees(max_infection, max_time, virtual_source, degrees, who_infected)
            # print('who_infected: ', who_infected)
            # print('source: ', source, "ml_estimate: ", ml_estimate)
        ml_correct[timesteps] = (ml_estimate == source)
        tot_num_infected[timesteps] = num_infected
        
        timesteps += 1
        # print('num infected', num_infected)
    
    original_m = m
    additional_hops = []
    additional_pd = []
    for i in range(additional_time):
        hopped = (random.random() < utilities.compute_alpha(m, timesteps, max_infection + 1))
        additional_hops.append(hopped)
        additional_pd.append(estimation.pd_additional_time(max_time, num_infected, additional_hops, max_infection + 1, original_m))
        if hopped:
            m += 1
        timesteps += 1
    
    infection_details = (tot_num_infected, infection_pattern, who_infected, additional_hops)
    return infection_details, ml_correct, additional_pd

# Semi-distributed adaptive diffusion over predefined random tree    
def infect_nodes_adaptive_planned_irregular_tree(source, max_time, max_infection, degrees_rv, known_degrees):
    '''infect_nodes_adaptive_planned_irregular_tree runs the spreading model over an irregular
    tree when you already know the structure of the random tree.
    
    Arguments:
        source: The ID of the source node (usually 0)
        max_time: The maximum amount of timesteps to run the algorithm
        max_infection: The maximum number of nodes an infected node will infect in any timestep
        degrees_rv: A random variable that describes the graph degree
        known_degrees: Tells the degrees of nodes in the irregular tree
      
    Returns:
        infection_details: A list of the characteristics of the infection:
            - tot_num_infected: Total number of infected nodes at each timestep
            - infection_pattern: Which nodes got infected in which order
            - who_infected: Adjacency matrix of infected subgraph
        ml_correct: List of whether the ML estimate was right in each timestep
        rand_leaf_correct: List of whether the random leaf estimate was right in each timestep
        known_degrees: The degrees of nodes in the graph that got added during infection'''
    
    timesteps = 0
    
    # initially the virtual source and the true source are the same
    virtual_source = source

    # ML estimate
    ml_correct = [0 for i in range(max_time)]
    rand_leaf_correct = [0 for i in range(max_time)]
    tot_num_infected = [0 for i in range(max_time)]
    num_infected = 0
    
    who_infected = [[]]
    degrees = [known_degrees.pop(0)]
    
    while timesteps < max_time:
    
            
        if timesteps == 0:
            virtual_source = 1
            previous_vs = 0
            
            # infect twice in one direction, always
            degrees, who_infected, known_degrees = utilities.infect_nodes_randtree(source, [virtual_source], degrees, degrees_rv, who_infected, known_degrees)
            infection_pattern, who_infected, known_degrees = pass_branch_message_randtree(source, virtual_source, degrees, degrees_rv, who_infected, known_degrees)
            m = 1       # the virtual source is always going to be 1 hop away from the true source
            
        else:
            current_neighbors = [k for k in who_infected[virtual_source]]

            if random.random() < utilities.compute_alpha(m,timesteps,max_infection):     # with probability alpha, spread symmetrically (keep the virtual source where it is)
                # branch once in every direction
                for neighbor in current_neighbors:
                    degrees, who_infected, known_degrees = pass_branch_message_randtree(virtual_source, neighbor, degrees, degrees_rv, who_infected, known_degrees)
            
            else:           # spread asymmetrically
                # find a direction to move
                virtual_source_candidate = previous_vs
                while virtual_source_candidate == previous_vs:
                    # virtual_source_candidate, current_neighbors, new_vs_likelihood = pick_random_elements(current_neighbors,1)
                    virtual_source_candidate = current_neighbors.pop(0)
                # virtual_source_candidate = virtual_source_candidate[0]
                previous_vs = virtual_source
                
                # the virtual source moves one more hop away from the true source
                m += 1;
                
                # branch twice in one direction
                degrees, who_infected, known_degrees = pass_branch_message_randtree(virtual_source, virtual_source_candidate, degrees, degrees_rv, who_infected, known_degrees)
                degrees, who_infected, known_degrees = pass_branch_message_randtree(virtual_source, virtual_source_candidate, degrees, degrees_rv, who_infected, known_degrees)
                
                virtual_source = virtual_source_candidate
            
        num_infected = len(who_infected)

        # Estimating the error
        # ML estimate
        ml_estimate = estimation.ml_estimate_irregular_trees(max_infection, max_time, virtual_source, degrees, who_infected)
        rand_leaf_estimate = estimation.rand_leaf_estimate(who_infected, degrees, timesteps)
        ml_correct[timesteps] = (ml_estimate == source)
        rand_leaf_correct[timesteps] = (rand_leaf_estimate == source)
        tot_num_infected[timesteps] = num_infected
        
        # print('who_infected',who_infected)
        # print('degrees',degrees)
        
        timesteps += 1
        
    infection_details = (tot_num_infected, infection_pattern, who_infected)
    return infection_details, ml_correct, rand_leaf_correct, known_degrees

# Fully distributed adaptive diffusion spreading algorithm for a line, applied to random trees
def infect_nodes_line_adaptive(source, max_time, max_infection, degrees_rv):

    timesteps = 0
    
    # initially the virtual source and the true source are the same
    virtual_source = source

    tot_num_infected = [0 for i in range(max_time)]
    
    who_infected = [[]]
    degrees = degrees_rv.rvs(size=1).tolist()
            
    
    leaves = [source]
    leaves_dist_from_source = [0];
    num_infected = 1
    
    jordan_correct = [0 for i in range(max_time)]
    rumor_correct = [0 for i in range(max_time)]
    # We don't know how to compute the ML in this case
    
    blocked = False
    # print('trial')
    while timesteps < max_time:
        new_leaves = []
        new_leaves_dist_from_source = []
        # For each edge emanating from each leaf, flip a coin
        for leaf, leaf_dist_from_source in zip(leaves, leaves_dist_from_source):
            neighbors = [i for i in range(num_infected, num_infected + max_infection + 1 - len(who_infected[leaf]))]
            
            coin_flips = [random.random() for i in range(len(neighbors))]
            ratio = (float(leaf_dist_from_source + 1) / (timesteps+2))
            # Figure out how many leaves actually got infected
            infections = sum([i < ratio for i in coin_flips])
            infected_neighbors = neighbors[:infections]
            degrees, who_infected = utilities.infect_nodes_randtree(leaf, infected_neighbors, degrees, degrees_rv, who_infected)[:2]
            
            # Update the new leaves, and their distances from the source
            new_leaves += infected_neighbors
            new_leaves_dist_from_source += [(leaf_dist_from_source + 1) for i in range(infections)]
            num_infected += infections
            # Make sure the original leaf is still a "leaf" if all its neighbors did not get infected
            if len(who_infected[leaf]) < (max_infection + 1):
                new_leaves += [leaf]
                new_leaves_dist_from_source += [leaf_dist_from_source]

        # Sort the new leaves
        new_leaves, new_leaves_dist_from_source = zip(*sorted(zip(new_leaves, new_leaves_dist_from_source)))
        leaves = new_leaves
        leaves_dist_from_source = new_leaves_dist_from_source
        
        # Jordan-centrality estimate
        jordan_estimate = estimation.jordan_centrality(who_infected)
        jordan_correct[timesteps] = (jordan_estimate == source)
        
        # Rumor centrality estimate
        rumor_estimate = estimation.rumor_centrality(who_infected)
        rumor_correct[timesteps] = (rumor_estimate == source)
                
        results = (jordan_correct, rumor_correct)
        tot_num_infected[timesteps] = num_infected
        timesteps += 1
        
    infection_details = (tot_num_infected, who_infected)
    return infection_details, results

# def do_cprofile(func):
    # def profiled_func(*args, **kwargs):
        # profile = cProfile.Profile()
        # try:
            # profile.enable()
            # result = func(*args, **kwargs)
            # profile.disable()
            # return result
        # finally:
            # profile.print_stats()
    # return profiled_func

# @do_cprofile
# Adaptive, semi-distributed diffusion with known adjacency matrix (used for datasets)
def infect_nodes_adaptive_diff(source, adjacency, max_time, max_infection):
    num_nodes = len(adjacency)
    timesteps = 0
    
    # initially the virtual source and the true source are the same
    virtual_source = source
    virtual_source_candidate = virtual_source
    
    infection_pattern = [0]*num_nodes
    infection_pattern[source] = 1
    dist_from_source = [-1]*num_nodes
    dist_from_source[source] = 0
    who_infected = [[] for i in range(num_nodes)]
    jordan_correct = [0 for i in range(max_time)]
    rumor_correct = [0 for i in range(max_time)]
    ml_correct = [0 for i in range(max_time)]
    ml_distances = [[] for i in range(max_time)]
        
    num_infected = []
    
    blocked = False
    
    while timesteps < max_time:
        # print('time', timesteps)
    
        # in odd timesteps, choose a direction to expand in
        if timesteps == 0:
            current_neighbors = [k for k in adjacency[source]]
            virtual_source_candidate, current_neighbors, source_likelihood = pick_random_elements(current_neighbors,1)
            previous_vs = virtual_source
            
            # infect twice in one direction, always
            infection_pattern, who_infected, dist_from_source = utilities.infect_nodes(source, virtual_source_candidate, infection_pattern, who_infected, dist_from_source)
            virtual_source_candidate = virtual_source_candidate[0]
            infection_pattern, who_infected, dist_from_source = pass_branch_message(source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected, dist_from_source)
            virtual_source = virtual_source_candidate
            m = 1       # the virtual source is always going to be 1 hop away from the true source at T=1
            
        else:
            current_neighbors = [k for k in who_infected[virtual_source]]
            # if (len(current_neighbors) < 2) or (random.random() < utilities.compute_alpha(m,timesteps,max_infection)):     # with probability alpha, spread symmetrically (keep the virtual source where it is)
            if (len(current_neighbors) < 2):
                # if there is nowhere for the virtual source to move, keep it where it is
                if len(current_neighbors) < 1:
                    blocked = True
                    print('Blocked. Exiting.')
                    break
                    
                # branch once in every direction
                for neighbor in current_neighbors:
                    infection_pattern, who_infected, dist_from_source = pass_branch_message(virtual_source, neighbor, infection_pattern, adjacency, max_infection, who_infected, dist_from_source)
                # take care of getting stuck
                if len(current_neighbors) == 1:
                    infection_pattern, who_infected, dist_from_source = pass_branch_message(virtual_source, current_neighbors[0], infection_pattern, adjacency, max_infection, who_infected, dist_from_source)
                    previous_vs = virtual_source
                    virtual_source = current_neighbors[0]
            
            else:           # spread asymmetrically
                # find a direction to move
                virtual_source_candidate = [previous_vs]
                while virtual_source_candidate[0] == previous_vs:
                    virtual_source_candidate, current_neighbors, new_vs_likelihood = pick_random_elements(current_neighbors,1)
                virtual_source_candidate = virtual_source_candidate[0]
                previous_vs = virtual_source
                # the virtual source moves one more hop away from the true source
                m += 1;
            
                # branch twice in one direction
                infection_pattern, who_infected, dist_from_source = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected, dist_from_source)
                infection_pattern, who_infected, dist_from_source = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected, dist_from_source)
                
                virtual_source = virtual_source_candidate
        # print('Adjacency at time ', timesteps)
        # utilities.print_adjacency(who_infected, adjacency)
        # print('\n')
        num_infected = num_infected + [sum(infection_pattern)]

        # Jordan-centrality estimate
        # jordan_estimate = estimation.jordan_centrality(who_infected)
        # jordan_correct[timesteps] = (jordan_estimate == source)
        jordan_correct[timesteps] = 0

        # Rumor centrality estimate
        # rumor_estimate = estimation.rumor_centrality(who_infected)
        # rumor_correct[timesteps] = (rumor_estimate == source)
        rumor_correct[timesteps] = 0
        
        # ML estimate
        ml_leaf, likelihoods, ml_distance = estimation.max_likelihood(who_infected, virtual_source, adjacency, max_infection, dist_from_source, source)
        ml_correct[timesteps] = (ml_leaf == source)
        ml_distances[timesteps] = ml_distance
        
        results = (jordan_correct, rumor_correct, ml_correct, ml_distances)
        timesteps += 1
    return num_infected, infection_pattern, who_infected, results
    
def compute_permutation_likelihood(T,m):
    return float(prob_G_given_m_and_T(T,m)) / pow(d,m)
    
def prob_G_given_m_and_T(T,m):
    if m < 1:
        print('MAJOR ERROR, m cannot be less than 1')
        return -1
    prob = float(1-utilities.compute_alpha(1,T,d))/pow(d,m)
    for i in range(2,m+1):
        prob = prob * utilities.compute_alpha(i,T,d)
    return prob

    
# Pramod's deterministic tree-shaped spreading algorithm.
def infect_nodes_deterministic(source, adjacency, greedy = False):
    num_nodes = len(adjacency)
    num_infected = [];
    infection_pattern = [0]*num_nodes;
    infection_pattern[source] = 2;
    infection_height = [0]*num_nodes;
    infection_height[source] = 0;
    
    
    old_infecting_nodes = [source]
    new_infecting_nodes = []
    
    blocked = False
    infected_nodes = 1
    while infected_nodes < num_nodes:
        print('old infecting', old_infecting_nodes)
        for node in old_infecting_nodes:
            current_neighbors = [k for k in adjacency[node] if infection_pattern[k]==0]
            num_neighbors = len(current_neighbors)
            # if len(current_neighbors) < 1:
                # blocked = True
                # break
            # if you're an up node, create an up and a down node
            if infection_pattern[node] == 1:
                if num_neighbors < 1:
                    print("Blocked. Neighbors are: ", adjacency[node])
                    blocked = True
                    break
                if greedy:
                    up_elements, current_neighbors = pick_random_elements_greedy(current_neighbors,1,adjacency,infection_pattern,node)[:2]
                else:
                    up_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                for item in up_elements:
                    infection_pattern[item] = 1
                    infection_height[item] = infection_height[node] + 1
                down_elements = []
                if num_neighbors > 1:
                    down_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                    for item in down_elements:
                        infection_pattern[item] = -1
                        infection_height[item] = infection_height[node] - 1
                    infected_nodes += 2
                else:
                    infected_nodes += 1
                new_infecting_nodes = new_infecting_nodes + up_elements + down_elements
            # if you're a down node, create two down nodes
            elif (infection_pattern[node] == -1) and (infection_height[node] > 0):
                if num_neighbors < 3:
                    down_elements, current_neighbors = pick_random_elements(current_neighbors,num_neighbors)[:2]
                    infected_nodes += num_neighbors
                else: 
                    down_elements, current_neighbors = pick_random_elements(current_neighbors,2)[:2]
                    infected_nodes += 2
                for item in down_elements:
                    infection_pattern[item] = -1
                    infection_height[item] = infection_height[node] - 1
                new_infecting_nodes = new_infecting_nodes + down_elements
            # if you're the true source, just infect one up node
            elif infection_pattern[node]==2:
                if num_neighbors < 1:
                    print("Blocked. Neighbors are: ", adjacency[node])
                    blocked = True
                    break
                up_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                for item in up_elements:
                    infection_pattern[item] = 1
                    infection_height[item] = 1
                infected_nodes += 1
                new_infecting_nodes = new_infecting_nodes + up_elements
        num_infected = num_infected + [infected_nodes]
        if blocked:
            break
            
        # swap the arrays
        print('new infecting', new_infecting_nodes)
        print('infection status: ', [infection_pattern[rr] for rr in new_infecting_nodes])
        print('infection height: ', [infection_height[rr] for rr in new_infecting_nodes])
        old_infecting_nodes, new_infecting_nodes = new_infecting_nodes, old_infecting_nodes
        new_infecting_nodes = []
    print("Infected ", infected_nodes, " / ", num_nodes, " total nodes.")
    return num_infected, infection_pattern

# Pramod's deterministic tree-shaped spreading algorithm.
def infect_nodes_deterministic_reinfect(source, adjacency):
    num_nodes = len(adjacency)
    num_infected = [];
    infection_pattern = [0]*num_nodes;
    infection_pattern[source] = 2;
    infection_height = [0]*num_nodes;
    infection_height[source] = 0;
    
    
    old_infecting_nodes = [source]
    new_infecting_nodes = []
    
    blocked = False
    infected_nodes = 1
    counter = 0
    while (infected_nodes < num_nodes) and (counter < 30):
        counter += 1
        print(counter)
        # print('old infecting', old_infecting_nodes)
        for node in old_infecting_nodes:
            current_neighbors = [k for k in adjacency[node] if ((infection_height[k] < 1) or (infection_height[k] != max(infection_height)))
                                                            and (k not in old_infecting_nodes)]
            num_neighbors = len(current_neighbors)
            if num_neighbors < 1:
                print("Blocked. Neighbors are: ", adjacency[node], [(infection_pattern[k],infection_height[k]) for k in adjacency[node]])
                blocked = True
                break
            # if you're an up node, create an up and a down node
            if infection_pattern[node] == 1:
                up_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                
                # print('up_elements', up_elements)
                for item in up_elements:
                    if (infection_pattern[item] == 0):
                        infected_nodes += 1
                    infection_pattern[item] = 1
                    # print("ones!")
                    infection_height[item] = infection_height[node] + 1
                down_elements = []
                if num_neighbors > 1:
                    # down_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                    down_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                    for item in down_elements:
                        if (infection_pattern[item] == 0):
                            infected_nodes += 1
                        infection_pattern[item] = -1
                        infection_height[item] = infection_height[node] - 1
                new_infecting_nodes = new_infecting_nodes + up_elements + down_elements
            # if you're a down node, create two down nodes
            elif (infection_pattern[node] == -1) and (infection_height[node] > 0):
                if num_neighbors < 2:
                    down_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                else: 
                    down_elements, current_neighbors = pick_random_elements(current_neighbors,2)[:2]
                for item in down_elements:
                    if (infection_pattern[item] == 0):
                        infected_nodes += 1
                    infection_pattern[item] = -1
                    infection_height[item] = infection_height[node] - 1
                new_infecting_nodes = new_infecting_nodes + down_elements
            # if you're the true source, just infect one up node
            elif infection_pattern[node]==2:
                if num_neighbors < 1:
                    print("Blocked. Neighbors are: ", adjacency[node])
                    blocked = True
                    break
                up_elements, current_neighbors = pick_random_elements(current_neighbors,1)[:2]
                for item in up_elements:
                    if (infection_pattern[item] == 0):
                        infected_nodes += 1
                    infection_pattern[item] = 1
                    infection_height[item] = 1
                new_infecting_nodes = new_infecting_nodes + up_elements
        num_infected = num_infected + [infected_nodes]
        if blocked:
            break
            
        # swap the arrays
        # print('new infecting', new_infecting_nodes)
        # print('infection status: ', [infection_pattern[rr] for rr in new_infecting_nodes])
        # print('infection height: ', [infection_height[rr] for rr in new_infecting_nodes])
        old_infecting_nodes, new_infecting_nodes = new_infecting_nodes, old_infecting_nodes
        new_infecting_nodes = []
    print("Infected ", infected_nodes, " / ", num_nodes, " total nodes.")
    return num_infected, infection_pattern
    
def pick_random_elements(neighbors,num_neighbors):
    # remove 'num_neighbors' random elements from the set neighbors, and return them along with the shortened list
    # Inputs
    #       neighbors:          the list of neighbors
    #       num_neighbors:      the number of elements to pick from the list
    #
    # Outputs
    #       random_elements     the neighbors chosen
    #       neighbors           the updated neighbors, without random_elements
    
    random_elements = []
    for i in range(num_neighbors):
        random_element = random.choice(neighbors)
        neighbors.remove(random_element)
        random_elements.append(random_element)
    if neighbors:
        likelihood = float(num_neighbors) / len(neighbors)
    else:
        likelihood = 1
    random_elements.sort()
    return random_elements, neighbors, likelihood

def pick_random_elements_greedy(neighbors,num_neighbors,adjacency,infection_pattern,node):
    # remove 'num_neighbors' random elements from the set neighbors, and return them along with the shortened list.
    # This time, pick a neighbor weighted by the number of uninfected neighbors.
    # Inputs
    #       neighbors:          the list of neighbors
    #       num_neighbors:      the number of elements to pick from the list
    #
    # Outputs
    #       random_elements     the neighbors chosen
    #       neighbors           the updated neighbors, without random_elements
    
    random_elements = []
    # Find which neighbor has the fewest infected friends, and give it the up chain
    surroundedness = []
    for neighbor in neighbors:
        infected_neighbors = sum([1 for k in adjacency[neighbor] if ((k!=node) and (infection_pattern[k]==0))])# / len(adjacency[neighbor])
        # surroundedness = surroundedness + [infected_neighbors]
        mutual_friends = sum([1 for k in adjacency[neighbor] if ((k!=node) and (k in neighbors))]) #/ len(adjacency[neighbor])
        surroundedness = surroundedness + [mutual_friends * max(infected_neighbors,1)]
    # sorted_neighbors = [x for (y,x) in sorted(zip(surroundedness,neighbors))]
    
    print("Sorted neighbors", sorted(zip(surroundedness,neighbors)))
    sorted_neighbors = [x for (y,x) in sorted(zip(surroundedness,neighbors))]
    if sum(surroundedness) == 0:
        total = 1
    else:
        total = sum(surroundedness)
    surroundedness = [1 - i/float(total) for i in surroundedness]
    virtual_source_rv = stats.rv_discrete(name='rv_discrete', values=(neighbors, surroundedness))
        
    for i in range(num_neighbors):
        element = virtual_source_rv.rvs(size=1)
        neighbors.remove(element)
        # element = sorted_neighbors.pop()
        random_elements.append(element)
    if neighbors:
        likelihood = float(num_neighbors) / len(neighbors)
    else:
        likelihood = 1
    random_elements.sort()
    # return random_elements, sorted_neighbors, likelihood
    return random_elements, neighbors, likelihood
    
def pass_branch_message_infinite_tree(source, recipient, adjacency, max_degree):
    # pass an instruction to branch from the source to the leaves overn an infinite tree
    # Inputs
    #       source:             source of the infection
    #       recipient:          array of child ids
    #       adjacency:          adjacency relations for the infected subgraph
    #       max_degree:         degree of each node in the regular tree     
    #
    # Outputs
    #       adjacency          (updated)

    leaf = True
    neighbors = adjacency[recipient]
    for neighbor in neighbors:
        if not neighbor == source:
            leaf = False
            adjacency =  pass_branch_message_infinite_tree(recipient, neighbor, adjacency, max_degree)
            
    if leaf:
        adjacency = utilities.infect_nodes_infinite_tree(recipient, max_degree-1, adjacency)
    return adjacency

def pass_branch_message_randtree(source, recipient, degrees, degrees_rv, who_infected, known_degrees = []):
    # pass an instruction to branch from the source to the leaves
    # Inputs
    #       source:             source of the infection
    #       recipient:          array of child ids
    #       degrees:            the degree of each node in the underlying irregular tree
    #       degrees_rv:          a random variable describing the degree distribution
    #       who_infected:       adjacency relations for the infected subgraph
    #
    # Outputs
    #       degrees             (updated)
    #       who_infected        (updated)
    
    leaf = True
    
    # pass to the neighbors who are your neighbors in the true infection graph
    neighbors = [k for k in who_infected[recipient] if (not k == source)]
    neighbors.sort() # this should be commented for a randomized distribution!
    
    for neighbor in neighbors:
        leaf = False
        degrees, who_infected, known_degrees =  pass_branch_message_randtree(recipient, neighbor, degrees, degrees_rv, who_infected, known_degrees)
            
    if leaf:
        to_infect = degrees[recipient] - 1
        neighbors = ([k+len(degrees) for k in range(to_infect)])
        degrees,who_infected, known_degrees = utilities.infect_nodes_randtree(recipient, neighbors, degrees, degrees_rv, who_infected, known_degrees)
    return degrees, who_infected, known_degrees
    
def pass_branch_message(source, recipient, infection_pattern, adjacency, max_infection, who_infected, dist_from_source):
    # pass an instruction to branch from the source to the leaves
    # Inputs
    #       source:             source of the infection
    #       recipient:          array of child ids
    #       infection_pattern:  binary array describing whether each node is already infected or not
    #       adjacency:          adjacency relations for the underlying network
    #       max_infection:      maximum number of people who can be infected by a single node     
    #       who_infected:       adjacency relations for the infected subgraph
    #
    # Outputs
    #       infection_pattern   (updated)
    #       who_infected        (updated)
    
    leaf = True
    
    # pass to the neighbors who are your neighbors in the true infection graph
    neighbors = [k for k in who_infected[recipient] if (not k == source)]
    neighbors.sort() # this should be commented for a randomized distribution!
    
    for neighbor in neighbors:
        leaf = False
        infection_pattern, who_infected, dist_from_source =  pass_branch_message(recipient, neighbor, infection_pattern, adjacency, max_infection, who_infected, dist_from_source)
            
    if leaf:
        neighbors = [k for k in adjacency[recipient] if infection_pattern[k]==0]
        if len(neighbors) > max_infection:
            neighbors, remains, leaf_likelihood = pick_random_elements(neighbors,max_infection)
        infection_pattern,who_infected, dist_from_source = utilities.infect_nodes(recipient, neighbors, infection_pattern, who_infected, dist_from_source)
    return infection_pattern, who_infected, dist_from_source
    
