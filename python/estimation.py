import numpy as np
import random 
import utilities
import infectionUtils
import math
import scipy.misc as sc


def rumor_centrality_up(up_messages, who_infected, calling_node, called_node):
    if called_node == calling_node:
        for i in who_infected[called_node]:
                up_messages = rumor_centrality_up(up_messages, who_infected, called_node, i)
    elif len(who_infected[called_node]) == 1:   # leaf node
        up_messages[calling_node] += 1 # check
    else:
        for i in who_infected[called_node]:
            if i != calling_node:
                up_messages = rumor_centrality_up(up_messages, who_infected, called_node, i)
        up_messages[calling_node] += up_messages[called_node]
    return up_messages        

def rumor_centrality_down(down_messages, up_messages, who_infected, calling_node, called_node):
    if called_node == calling_node:
        for i in who_infected[called_node]:
            down_messages = rumor_centrality_down(down_messages, up_messages, who_infected, called_node, i)   
    else:
        down_messages[called_node] = down_messages[calling_node]*(float(up_messages[called_node])/(len(who_infected)-up_messages[called_node]))
        for i in who_infected[called_node]:
            if i != calling_node:
                down_messages = rumor_centrality_down(down_messages, up_messages, who_infected, called_node, i)
    return down_messages 


def rumor_centrality(who_infected):
    # computes the estimate of the source based on rumor centrality
    initial_node = 0       # can use arbitrary initial index
    up_messages = [1]*len(who_infected) 
    down_messages = [1]*len(who_infected)
    up_messages = rumor_centrality_up(up_messages,who_infected,initial_node,initial_node)
    down_messages = rumor_centrality_down(down_messages,up_messages,who_infected,initial_node,initial_node)
    max_down = max(down_messages)
    max_down_ind = [i for i, j in enumerate(down_messages) if j == max_down]
    return max_down_ind[random.randrange(0,len(max_down_ind),1)]

def jordan_centrality(adjacency):
    # computes the estimate of the source based on jordan centrality
    num_nodes = len(adjacency)
    
    if num_nodes == 1:
        return 0
    elif num_nodes == 2:
        return int(random.random() < 0.5)

    # choose a root that is not a leaf
    while True:
        root = random.randint(0,num_nodes-1)
        if len(adjacency[root]) > 1:
            break
    # initialize lists
    up_count = [0 for i in range(num_nodes)]
    up_count_2nd = [0 for i in range(num_nodes)]
    up_id = [-1 for i in range(num_nodes)]
    up_id_2nd = [-1 for i in range(num_nodes)]
    
    # start passing messages
    up_count, up_count_2nd, up_id, up_id_2nd = get_messages(root, root, adjacency, up_count, up_count_2nd, up_id, up_id_2nd)
    
    # now find the jordan center
    idx_jordan = find_center(root, up_count, up_count_2nd, up_id, up_id_2nd)
    return idx_jordan
    
def get_messages(parent,node,adjacency,up_count,up_count_2nd,up_id,up_id_2nd):
    # obtains the message from node to parent
    
    children = adjacency[node]
    
    # deal with leaves
    if len(children) > 1:    
        # forward messages appropriately
        for child in children:     
            
            if child == parent:
                continue
            up_count, up_count_2nd, up_id, up_id_2nd = get_messages(node, child, adjacency, up_count, up_count_2nd, up_id, up_id_2nd)
            up_count, up_count_2nd, up_id, up_id_2nd = update_counts(node, child, up_count, up_count_2nd, up_id, up_id_2nd)
    
    return up_count, up_count_2nd, up_id, up_id_2nd
    
def update_counts(parent, node, up_count, up_count_2nd, up_id, up_id_2nd):
    # updates the message arrays for jordan centrality
    
    candidates = np.array([up_count[parent], up_count_2nd[parent], up_count[node]+1])
        
    sorted_vals = np.sort(candidates)
    indices = np.argsort(candidates)
    
    up_count[parent] = sorted_vals[-1]
    up_count_2nd[parent] = sorted_vals[-2]
    if (indices[-1] == 2): # the incoming max is the largest
        up_id_2nd[parent] = up_id[parent]
        up_id[parent] = node
        
    elif (indices[-2] == 2): # the incoming max is 2nd largest
        # up_id stays the same
        up_id_2nd[parent] = node
            
    return up_count, up_count_2nd, up_id, up_id_2nd
    
def find_center(root, up_count, up_count_2nd, up_id, up_id_2nd):
    # finds the center of the graph using jordan centrality metric, after messages have been passed
    idx = root
    l1 = up_count[idx]
    l2 = up_count_2nd[idx]
    breakflag = False
    while True:
        if (l1-l2) == 0:
            idx_jordan = idx
            break
        elif (abs(l1-l2)) == 1:
            if breakflag:
                # break ties randomly
                if random.random() < 0.5:
                    idx_jordan = idx
                break
            breakflag = True
            idx_jordan = idx
        elif breakflag:
            # there's no tie, so use the previously found index
            break
            
        idx = up_id[idx]
        l2 = l2 + 1
        l1 = l1 - 1

        if (l1 <= 1):
            idx_jordan = idx
            break
    return idx_jordan
    
# ML Estimate
def max_likelihood(who_infected, virtual_source, adjacency, max_infection, dist_from_source, source):
    # compute the maximum likelihood estimate of the source using only leaves
    # Inputs
    #       who_infected:       adjacency relations for the infected subgraph
    #       virtual_source:     virtual source of the infection (centroid of the infection subgraph)
    #       adjacency:          adjacency relations for the underlying network
    #       max_infection:      maximum number of people who can be infected by a single node     
    #       source:             true source of the rumor
    #
    # Outputs
    #       ml_estimate         maximum likelihood node
    #       likelihoods         likelihoods for all the nodes in teh graph
    
    paths = [[] for i in range(len(who_infected))]
    paths[virtual_source].append(virtual_source)
    paths = get_paths(paths, who_infected, adjacency, virtual_source, virtual_source)
    # maximum length from leaf to virtual source
    max_length = max([len(k) for k in paths])
    likelihoods = [float('-inf') for i in range(len(who_infected))]
    for i in range(len(who_infected)):
        # this algorithm only works for leaves in the infected subgraph, and they should be at the boundary of the graph
        if (len(who_infected[i]) == 1) and len(paths[i]) == max_length:
            vs_path = [k for k in paths[i]]
            likelihoods[i] = compute_graph_likelihood(i, who_infected, adjacency, vs_path, max_infection)
            if likelihoods[i] == float('-inf'):
                print('negative inf', i,vs_path, adjacency[i])
                print('likelihood: ',math.log(1.0/len(adjacency[i])))
            # else:
                # print('not inf', i,vs_path)
                
    max_likelihood = max(likelihoods)
    indices = [i for i, x in enumerate(likelihoods) if x == max_likelihood]
    ml_estimate = random.choice(indices)
    # print('Finding the tree distance')
    tree_dist = dist_from_source[ml_estimate]
    # print('Tree distance: ', tree_dist)
    # print('Tree distance: ', tree_dist, '. Now finding real dist. ')
    # real_dist = get_estimate_dist(source, ml_estimate, adjacency)
    # print('Real distance: ', real_dist)
    # distances = [tree_dist, real_dist]
    distances = tree_dist
    return ml_estimate, likelihoods, distances
    
def compute_graph_likelihood(source, who_infected, adjacency, vs_path, max_infection):
    # compute the likelihood of the infected subgraph starting from node source
    # Inputs
    #       source:             assumed source of the rumor
    #       who_infected:       adjacency relations for the infected subgraph
    #       adjacency:          adjacency relations for the underlying network
    #       vs_path:            path that the virtual source takes in the graph
    #       max_infection:      maximum number of people who can be infected by a single node     
    #
    # Outputs
    #       likelihood          likelihood of the graph given source 
    
    if len(vs_path) == 1:
        print('The vs path is 1 hop!', source, vs_path)
        utilities.print_adjacency(who_infected, adjacency)
        return float('-inf')
    
    vs_path.pop()
    nodes = range(len(who_infected))
    new_infection_pattern = [0 for i in nodes]
    new_infection_pattern[source] = 1
    new_who_infected = [[] for i in nodes]

    # first element is just the source itself
    current_vs = vs_path.pop(0)    
    # log likelihood of the 1st node
    likelihood = math.log(1.0/len(adjacency[source]))
    
    # get the new vs
    current_vs = path.pop(0)
    new_infection_pattern, new_who_infected, tmp = infectionUtils.infect_nodes(source, [current_vs], new_infection_pattern, new_who_infected)
    
    # infect the neighbors of the new vs
    infected = [i for i in who_infected[current_vs]]
    infected.remove(source)
    
    new_infection_pattern, new_who_infected, tmp = infectionUtils.infect_nodes(current_vs, infected, new_infection_pattern, new_who_infected)
    likelihood += infect_set_likelihood(infected, adjacency[current_vs], new_infection_pattern, max_infection)
    
    while vs_path:
        new_infection_pattern, new_who_infected, likelihood = pass_branch_message_likelihood(current_vs, vs_path[0], new_infection_pattern, adjacency, max_infection, new_who_infected, who_infected, likelihood)
        current_vs = vs_path.pop(0)
        
    return likelihood

def pass_branch_message_likelihood(source, recipient, new_infection_pattern, adjacency, max_infection, new_who_infected, who_infected, likelihood):
    # pass an instruction to branch from the source to the leaves using the who_infected pattern, and compute the likelihood
    # Inputs
    #       source:             source of the infection
    #       recipient:          array of child ids
    #       new_infection_pattern:  binary array describing whether each node is already infected or not in the graph being built up
    #       adjacency:          adjacency relations for the underlying network
    #       max_infection:      maximum number of people who can be infected by a single node     
    #       new_who_infected:   adjacency relations for the infected subgraph being built up
    #       who_infected:       adjacency relations for the original infected subgraph
    #       likelihood:         the likelihood for this graph
    #
    # Outputs
    #       new_infection_pattern   (updated)
    #       new_who_infected        (updated)
    #       likelihood              (updated)
    
    leaf = True
    
    # pass to the neighbors who are your neighbors in the true infection graph
    neighbors = [k for k in new_who_infected[recipient] if (not k == source)]
    neighbors.sort()
    
    for neighbor in neighbors:
        leaf = False
        new_infection_pattern, new_who_infected, likelihood =  pass_branch_message_likelihood(recipient, neighbor, new_infection_pattern, adjacency, max_infection, new_who_infected, who_infected, likelihood)
            
    if leaf:
        neighbors = [k for k in who_infected[recipient] if not k==source]
        likelihood += infect_set_likelihood(neighbors, adjacency[recipient], new_infection_pattern, max_infection)
        new_infection_pattern, new_who_infected, tmp = infectionUtils.infect_nodes(recipient, neighbors, new_infection_pattern, new_who_infected)
        
    return new_infection_pattern, new_who_infected, likelihood
    
def infect_set_likelihood(infected, adjacency_choices, new_infection_pattern, max_infection):
    # likelihood of infecting a particular set given a state of the current infected subgraph
    # Inputs
    #       infected:           set of infected neighbor nodes
    #       adjacency_choices:  nodes available in the adjacency graph (infected and not)
    #       new_infection_pattern:  binary array describing whether each node is already infected or not in the graph being built up
    #       max_infection:      maximum number of people who can be infected by a single node     
    #
    # Outputs
    #       likelihood          numeric likelihood of this event happening
    choices = [i for i in adjacency_choices if new_infection_pattern[i] == 0]
    
    di = len(choices)
    num_infected = len(infected)
    
    if (di >= max_infection) and (num_infected < max_infection):
        likelihood = float('-inf')
    elif (di < max_infection) and (num_infected < di):
        likelihood = float('-inf')
    elif di < num_infected:
        likelihood = float('-inf')
    elif num_infected == 0:
        likelihood = 1
    else:
        likelihood = math.log(1.0 / utilities.nCk(di, num_infected))
        
    return likelihood

def get_paths(paths, who_infected, adjacency, called_node, calling_node):
    # gets the path from every node in the infected subtree to the virtual source of that tree
    # Inputs
    #       paths:              list of paths from each node to the V.S.
    #       who_infected:       adjacency relations for the infected subgraph
    #       adjacency:          adjacency relations for the underlying network
    #       called_node:        node to which the messages get passed
    #       calling_node:       the node that called this function
    #
    # Outputs
    #       paths              (updated)

    if not called_node == calling_node:
        paths[called_node] = [called_node] + paths[calling_node] 
    neighbors = who_infected[called_node]
    if len(neighbors) > 1:
        for neighbor in neighbors:
            if not neighbor == calling_node:
                paths = get_paths(paths, who_infected, adjacency, neighbor, called_node) 
    return paths

# ML over irregular infinite trees
def ml_message_passing_irregular_trees(d, depth, messages, degrees, who_infected, called, calling): 
    # d:            d+1 is the assumed regular degree (aka d_o)
    # depth:        distance from the virtual source
    # messages:     the likelihood messages stored at each node
    # degrees:      degree of each node in the irregular graph
    # who_infected: adjacency matrix of the infected subgraph
    # called:       which node received the message
    # calling:      which node passed the message

    if called == calling:
        for i in who_infected[called]:
            # (Prob. of choosing virtual source) x (Prob. infecting infected nodes)
            messages[i] = messages[calling] * (1.0 / (degrees[i])) * (d + 1)
            # print('msgs', messages)
            messages = ml_message_passing_irregular_trees(d, depth+1, messages, degrees, who_infected, i, called) 

    elif len(who_infected[called]) != 1: #if not a leaf
        for i in who_infected[called]:
            if i != calling:
                messages[i] = messages[called] * d * \
                             (float(degrees[called])/(degrees[called]-1)) * (1.0/(degrees[i]))
                # print('new messages in message passing: ', messages)
                messages = ml_message_passing_irregular_trees(d, depth+1, messages, degrees, who_infected, i, called) 
    return messages 

# ML over irregular infinite trees using weighted spreading
def ml_message_passing_irregular_trees_alt(d, depth, messages, degrees, who_infected, called, calling, degrees_rv, prev_prob=0): 
    prob = [float(degrees[calling]),sum([degrees[k] for k in who_infected[called]])]
    if called == calling:
        for i in who_infected[called]:
            # probability of choosing the vs instead of another neighbor
            messages[i] = messages[calling] * d * prob[0] / prob[1]
            messages = ml_message_passing_irregular_trees_alt(d, depth+1, messages, degrees, who_infected, i, called,degrees_rv, prob) 

    elif len(who_infected[called]) != 1: #if not a leaf
        correction_factor = float(prev_prob[1]) / (prev_prob[1]-degrees[called])
        for i in who_infected[called]:
             if i != calling:
                if len(who_infected[i])==1:
                    prob[1] = sum([degrees_rv.draw_values(1)[0] for k in range(degrees[i])])
                messages[i] = messages[called] * (d-1) * correction_factor * prob[0] / prob[1]
                messages = ml_message_passing_irregular_trees_alt(d, depth+1, messages, degrees, who_infected, i, called, degrees_rv,prob) 
    return messages 
    
    
def ml_estimate_irregular_trees(d, T, virtual_source, infected_nodes_degree, who_infected, degrees_rv=[], mode=0):        
    # Returns the ML estimate of the source over an irregular tree 
    # Inputs:
    #
    #       d:                      assumed regular degree - 1 (for computing the alphas)
    #       T:                      number of timesteps to run
    #       virtual_source:         the virtual source of the infected subtree
    #       infected_nodes_degree:  the degrees of all infected nodes
    #       who_infected:           adjacency matrix of the infected subgraph
    #       mode:                   which type of estimator to use (mode=0 -> use the regular estimator, mode=1-> use the alternative estimator
    #
    # Outputs:
    #
    #       ml_estimate:            the node with the highest likelihood of being the true source 

    
    # compute the probability of not passing the virtual source token
    p = 1.0
    # initializing the messages vector to likelihood 1
    messages = [p]*len(who_infected)
    # print('who infected', who_infected)

    # computing the likelihood of all the nodes via a message passing algorithm
    if mode==0:
        messages = ml_message_passing_irregular_trees(d, 0, messages, infected_nodes_degree, who_infected, virtual_source, virtual_source)
    else:
        messages = ml_message_passing_irregular_trees_alt(d, 0, messages, infected_nodes_degree, who_infected, virtual_source, virtual_source, degrees_rv)
    # print('MESSAGES ARE : ', messages)
    # print('DEGREES ARE : ', infected_nodes_degree)
    
    # the likelihood of the virtual source is equal to zero
    # print('messages', messages)
    messages[virtual_source] = 0
    # finding the likelihood of the ML estimate
    max_message = max(messages)
    # finding the indices of most likely nodes
    max_message_ind = [i for i, j in enumerate(messages) if j == max_message]
    # ml_estimate = max_message_ind[random.randrange(0,len(max_message_ind),1)]
    # print('choices', max_message_ind)
    ml_estimate = random.choice(max_message_ind)
    return ml_estimate
    
# estimates the ML depth of the true source using subsequent observations of the graph
# this currently only works on a line
def pd_additional_time(T, N_T, additional_hops, degree, true_depth):
    '''Compute the probability of detection using the ML estimator
        Args:
            T: The time at which the first snapshot was taken
            additional_hops: The hop/not hop pattern following
            degree: The degree of infection at each timestep
            true_depth: The depth of the true source in the tree at time T
        Outputs the probability of detection given this information
    '''

    # The depth of the final tree is T
    # From each source candidate, find the likelihood of additional hops
    depth_likelihood = []
    for m in range(1,T+1):
        likelihood = 1
        t = T
        new_m = m
        for hop in additional_hops:
            if hop:
                likelihood *= (1.0 - utilities.compute_alpha(new_m, t, degree)) / (degree - 1)
                new_m += 1
            else:
                likelihood *= utilities.compute_alpha(new_m, t, degree)
            t += 1
        depth_likelihood.append(likelihood)
    # Get the most likely depth(s) of the true source in the original tree
    max_likelihood = max(depth_likelihood)
    depth = [i+1 for i, x in enumerate(depth_likelihood) if x == max_likelihood]
    # Now estimate pd
    pd = 1.0 / (N_T - 1)
    if true_depth in depth:
        if True in additional_hops:
            # Find how many nodes were at that depth
            options = 0
            for d in depth:
                options += pow(degree - 1, d)
        else:
            options = 0
            for d in depth:
                options += (degree) * pow(degree - 1, d - 1)
        pd = 1.0 / options
    return pd

def rand_leaf_estimate(who_infected, degrees, T):
    n = len(who_infected)
    candidates = [i for i in range(n) if len(who_infected[i])==1]
    # print('num candidates',len(candidates))
    # if (len(candidates) != 7) and ( T > 0):
        # print('who infected',who_infected)
        # print('degrees', degrees)
    # if len(who_infected[0])>1:
        # print('0 is not a leaf!',who_infected)  

    rand_leaf_estimate = random.choice(candidates)
    return rand_leaf_estimate

def get_estimate_dist(source, destination, adjacency):
    ''' Finds the minimum distances between source and ml_estimate, both on the
    infection tree (who_infected) and on the underlying social network
    (adjacency).
    '''
    visited, queue = [], [[source,0]]
    while queue:
        vertex = queue.pop(0)
        if vertex[0] == destination:
            return vertex[1]
        if vertex[0] not in visited:
            visited.append(vertex[0])
            newnodes = [[i, vertex[1] + 1] for i in adjacency[vertex[0]] if i not in visited]
            queue.extend(newnodes)
    return -1