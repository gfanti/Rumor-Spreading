# node infection
import estimation
import estimation_spies
import random
from scipy import stats
import numpy as np
import utilities
import infectionUtils
import networkx as nx
import time
from infectionUtils import *
from spies import *
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
            adjacency = infectionUtils.infect_nodes_infinite_tree(source, 1, adjacency)
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

#Distributed diffusion over random trees, with iid spies in the network
def infect_nodes_diffusion_irregular_tree(source, max_time, degrees_rv, q = 0.5,
                                          spy_probability = 0.0, est_times = None,
                                          diffusion = False):
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
    
    num_infected = 0
    boundary_nodes = [source]
    who_infected = [[]]
    # degrees = degrees_rv.rvs(size=1).tolist()
    degrees = degrees_rv.draw_values(1)
    timestamps, active_nodes = [0],[0] # marks which nodes are candidate sources >=0 => valid, <0 => not valid
    spies, active_spies, infectors = [],[],[]

    while boundary_nodes:
    # while timesteps <= max_time:
    # while (1 in [active_nodes[item] for item in boundary_nodes]) and timesteps < 10:
        # print('time', timesteps)
    
        num_candidates = len(boundary_nodes)
        node = boundary_nodes.pop(0)
        # print('degrees[node]', node, len(degrees), len(who_infected))
        num_uninfected_neighbors = degrees[node] - len(who_infected[node])
        # num_to_infect = np.random.binomial(num_uninfected_neighbors, q)
        num_to_infect = num_uninfected_neighbors
        # print('num to infect', num_to_infect, 'out of ', num_uninfected_neighbors,'\n')
        # print('who_infected before: ', who_infected)
        if num_to_infect > 0:
            # neighbor_times = [max(j,0) + timestamps[node] for j in np.random.normal(1.0/q, 0.5, num_to_infect)]
            # neighbor_times = [max(j,0) + timestamps[node] for j in np.random.normal(1.0/q, 1.7321, num_to_infect)] # Truncated gaussian
            neighbor_times = [max(j,0) + timestamps[node] for j in np.random.geometric(q, num_to_infect)]
            # neighbor_times = [max(j,0) + timestamps[node] for j in np.random.normal(20, 5, num_to_infect)]
            neighbor_times = [neighbor_times[k] for k in range(num_to_infect) if neighbor_times[k] <= max_time]
            neighbors = ([k+len(degrees) for k in range(len(neighbor_times))])
            
            timestamps += neighbor_times
            
            degrees, who_infected = infectionUtils.infect_nodes_randtree(node, neighbors, degrees, degrees_rv, who_infected)[:2]
            # Choose spies with probability spy_probability
            new_spies = utilities.update_spies_diffusion(neighbors, spy_probability=spy_probability)
            spies += new_spies
            # update the infectors
            infectors += [node for item in neighbors if item in spies]
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
        
    print('num_infected: ', num_infected)
    
    infection_details, results = estimate_source_spies(max_time, est_times, source, who_infected,
                                                       num_infected, timestamps, spies,
                                                       active_nodes, active_spies, infectors, diffusion, q)
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
    # degrees = degrees_rv.rvs(size=1).tolist()
    degrees = degrees_rv.draw_values(1)
    
    while timesteps < max_time:
    
            
        if timesteps == 0:
            virtual_source = 1
            previous_vs = 0
            
            # infect twice in one direction, always
            degrees, who_infected = infectionUtils.infect_nodes_randtree(source, [virtual_source], degrees, degrees_rv, who_infected)[:2]
            infection_pattern, who_infected = pass_branch_message_randtree(source, virtual_source, degrees, degrees_rv, who_infected)[:2]
            m = 1       # the virtual source is always going to be 1 hop away from the true source
            
            # Account for any spies in the network
            # if spy_probability > 0:
                
            
        else:
            current_neighbors = [k for k in who_infected[virtual_source]]
            # with probability alpha, spread symmetrically (keep the virtual source where it is)
            # (when we have spies, we only spread asymmetrically)
            if random.random() < utilities.compute_alpha(m,timesteps,max_infection):
                # branch once in every direction
                for neighbor in current_neighbors:
                    degrees, who_infected = pass_branch_message_randtree(virtual_source, neighbor, degrees, degrees_rv, who_infected)[:2]
            # Otherwise, spread asymmetrically
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
        # The estimation for spies comes at the end, whereas we can do snapshot estimation on the fly
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
    
    # Deal with any additional snapshots
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
            degrees, who_infected, known_degrees = infectionUtils.infect_nodes_randtree(source, [virtual_source], degrees, degrees_rv, who_infected, known_degrees)
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
    # degrees = degrees_rv.rvs(size=1).tolist()
    degrees = degrees_rv.draw_values(1)
            
    
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
            degrees, who_infected = infectionUtils.infect_nodes_randtree(leaf, infected_neighbors, degrees, degrees_rv, who_infected)[:2]
            
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
            infection_pattern, who_infected, dist_from_source = infectionUtils.infect_nodes(source, virtual_source_candidate, infection_pattern, who_infected, dist_from_source)
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



