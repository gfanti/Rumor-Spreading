# node infection
import estimation
import random
from scipy import stats
import numpy as np
import utilities

# Semi-distributed adaptive diffusion over regular trees
def infect_nodes_adaptive_tree(source, adjacency, max_degree, max_time, alpha):
    '''infect_nodes_adaptive_irregular_tree runs the spreading model over an irregular
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

# Semi-distributed adaptive diffusion over irregular trees
def infect_nodes_adaptive_irregular_tree(source, max_time, max_infection,
                                              degrees_rv, alt = False):
    '''infect_nodes_adaptive_irregular_tree runs the spreading model over an irregular
    tree.
    
    Arguments:
        source: The ID of the source node (usually 0)
        max_time: The maximum amount of timesteps to run the algorithm
        max_infection: The maximum number of nodes an infected node will infect in any timestep
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
                    weights = [i/sum(weights) for i in weights]
                    virtual_source_rv = stats.rv_discrete(name='rv_discrete', values=(current_neighbors, weights))
                    virtual_source_candidate = virtual_source_rv.rvs(size=1)
                else: # Choose uniformly between all virtual sources
                    virtual_source_candidate = [previous_vs]
                    while virtual_source_candidate[0] == previous_vs:
                        virtual_source_candidate, current_neighbors, new_vs_likelihood = pick_random_elements(current_neighbors,1)
                    virtual_source_candidate = virtual_source_candidate[0]
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
        ml_correct[timesteps] = (ml_estimate == source)
        tot_num_infected[timesteps] = num_infected
        
        timesteps += 1
        
    infection_details = (tot_num_infected, infection_pattern, who_infected)
    return infection_details, ml_correct

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

# Adaptive, semi-distributed diffusion with known adjacency matrix (used for datasets)
def infect_nodes_adaptive_diff(source, adjacency, max_time, max_infection, stepsize):
    num_nodes = len(adjacency)
    timesteps = 0
    
    # initially the virtual source and the true source are the same
    virtual_source = source
    virtual_source_candidate = virtual_source
    
    infection_pattern = [0]*num_nodes;
    infection_pattern[source] = 1;
    who_infected = [[] for i in range(num_nodes)]
    jordan_correct = [0 for i in range(max_time)]
    rumor_correct = [0 for i in range(max_time)]
    ml_correct = [0 for i in range(max_time)]
    
    # compute the bins
    max_nodes = int(N_nodes(max_time, max_infection+1))
    bins = [i for i in range(1,max_nodes+stepsize+1,stepsize)]
    jordan_detected = [0 for i in range(len(bins))]
    rumor_detected = [0 for i in range(len(bins))]
    instances = [0 for i in range(len(bins))]
    
    num_infected = 0
    
    blocked = False
    
    while timesteps < max_time:
        # print('time', timesteps)
    
        # in odd timesteps, choose a direction to expand in
        if timesteps == 0:
            current_neighbors = [k for k in adjacency[source]]
            virtual_source_candidate, current_neighbors, source_likelihood = pick_random_elements(current_neighbors,1)
            previous_vs = virtual_source
            
            # infect twice in one direction, always
            infection_pattern, who_infected = utilities.infect_nodes(source, virtual_source_candidate, infection_pattern, who_infected)
            virtual_source_candidate = virtual_source_candidate[0]
            infection_pattern, who_infected = pass_branch_message(source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
            virtual_source = virtual_source_candidate
            m = 1       # the virtual source is always going to be 1 hop away from the true source
            
        else:
            current_neighbors = [k for k in who_infected[virtual_source]]

            if random.random() < utilities.compute_alpha(m,timesteps,max_infection):     # with probability alpha, spread symmetrically (keep the virtual source where it is)
                # print('stay')
                # if there is nowhere for the virtual source to move, keep it where it is
                if len(current_neighbors) < 1:
                    blocked = True
                    break
                    
                # branch once in every direction
                for neighbor in current_neighbors:
                    infection_pattern, who_infected = pass_branch_message(virtual_source, neighbor, infection_pattern, adjacency, max_infection, who_infected)
                # take care of getting stuck
                if len(current_neighbors) == 1:
                    infection_pattern, who_infected = pass_branch_message(virtual_source, current_neighbors[0], infection_pattern, adjacency, max_infection, who_infected)
                    previous_vs = virtual_source
                    virtual_source = current_neighbors[0]
                    m -= 1
            
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
                infection_pattern, who_infected = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
                infection_pattern, who_infected = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
                
                virtual_source = virtual_source_candidate
            
        num_infected = sum(infection_pattern)

        # Estimating the error
        bin_idx = next(x[0] for x in enumerate(bins) if x[1] > num_infected)
        instances[bin_idx] += 1

        # Jordan-centrality estimate
        jordan_estimate = estimation.jordan_centrality(who_infected)
        jordan_correct[timesteps] = (jordan_estimate == source)
        jordan_detected[bin_idx] += (jordan_estimate == source)
        jordan_results = (bins,instances,jordan_detected, jordan_correct)

        # Rumor centrality estimate
        rumor_estimate = estimation.rumor_centrality(who_infected)
        rumor_correct[timesteps] = (rumor_estimate == source)
        rumor_detected[bin_idx] += (rumor_estimate == source)
        rumor_results = (bins,instances,rumor_detected, rumor_correct)
        
        # ML estimate
        ml_estimate, likelihoods = estimation.max_likelihood(who_infected, virtual_source, adjacency, max_infection, source)
        ml_correct[timesteps] = (ml_estimate == source)
        
        results = (jordan_results, rumor_results, ml_correct)
        
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
def infect_nodes_deterministic(source, adjacency):
    num_nodes = len(adjacency)
    num_infected = 0;
    infection_pattern = [0]*num_nodes;
    infection_pattern[source] = 1;
    
    
    old_infecting_nodes = [source]
    new_infecting_nodes = []
    
    blocked = False
    while num_infected < num_nodes:
        
        for node in old_infecting_nodes:
            current_neighbors = [k for k in adjacency[node] if infection_pattern[k]==0]
            if len(current_neighbors) < 2:
                blocked = True
                break
            # if you're an up node, create an up and a down node
            if infection_pattern[node] == 1:
                up_elements, current_neighbors = pick_random_elements(current_neighbors,1)
                for item in up_elements:
                    infection_pattern[item] = 1
                down_elements, current_neighbors = pick_random_elements(current_neighbors,1)
                for item in down_elements:
                    infection_pattern[item] = -1
                new_infecting_nodes = new_infecting_nodes + up_elements + down_elements
            # if you're a down node, create two down nodes
            elif infection_pattern[node] == -1:
                down_elements, current_neighbors = pick_random_elements(current_neighbors,2)
                for item in down_elements:
                    infection_pattern[item] = -1
                new_infecting_nodes = new_infecting_nodes + down_elements
            # otherwise, something's wrong
            else:
                print('ERROR')
                
            num_infected += 2
        if blocked:
            break
            
        # swap the arrays
        old_infecting_nodes, new_infecting_nodes = new_infecting_nodes, old_infecting_nodes
    
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
    
def pass_branch_message(source, recipient, infection_pattern, adjacency, max_infection, who_infected):
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
        infection_pattern, who_infected =  pass_branch_message(recipient, neighbor, infection_pattern, adjacency, max_infection, who_infected)
            
    if leaf:
        # print('leaf')
        neighbors = [k for k in adjacency[recipient] if infection_pattern[k]==0]
        if len(neighbors) > max_infection:
            neighbors, remains, leaf_likelihood = pick_random_elements(neighbors,max_infection)
        infection_pattern,who_infected = utilities.infect_nodes(recipient, neighbors, infection_pattern, who_infected)
    return infection_pattern, who_infected
    
