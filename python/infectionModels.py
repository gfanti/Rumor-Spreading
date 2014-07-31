# node infection
import estimation
import random

def infect_nodes_adaptive_diff_tree(source, adjacency, max_degree, max_time, alpha, beta):
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
            adjacency = infect_nodes_infinite_tree(source, 1, adjacency)
            virtual_source_candidate = 1
            timesteps += 1
            continue
            
        if timesteps % 2 == 1:
            current_neighbors = [k for k in adjacency[virtual_source]]
            if len(current_neighbors) < 1:
                # print('Blocked. Total neighbors: ',adjacency[node])
                # print('Available neighbors: ',current_neighbors)
                blocked = True
                break
                
            virtual_source_candidate, current_neighbors = pick_random_elements(current_neighbors,1)
            virtual_source_candidate = virtual_source_candidate[0]
            # include beta here
            # print('passing from ',virtual_source, ' to ', virtual_source_candidate)
            adjacency = pass_branch_message_infinite_tree(virtual_source, virtual_source_candidate, adjacency, max_degree)
            
            
        # in even timesteps, choose whether to move the virtual source
        else:
            if random.random() < alpha:
                # print('passing from ',virtual_source, ' to ', virtual_source_candidate)
                adjacency = pass_branch_message_infinite_tree(virtual_source, virtual_source_candidate, adjacency, max_degree)
                virtual_source = virtual_source_candidate
            else:
                # print('passing from ',virtual_source_candidate, ' to ', virtual_source)
                adjacency = pass_branch_message_infinite_tree(virtual_source_candidate, virtual_source, adjacency, max_degree)
    
        num_infected = len(adjacency)
        print('num infected', num_infected)
        timesteps += 1
        
    return adjacency, num_infected
    
def infect_nodes_adaptive_diff(source, adjacency, max_time, max_infection):
    num_nodes = len(adjacency)
    timesteps = 0
    
    # initially the virtual source and the true source are the same
    virtual_source = source
    virtual_source_candidate = virtual_source
    
    infection_pattern = [0]*num_nodes;
    infection_pattern[source] = 1;
    who_infected = [[] for i in range(num_nodes)]
    jordan_correct = [0 for i in range(max_time)]
    num_infected = 0
    
    blocked = False
    
    while timesteps < max_time:
        # print('time', timesteps)
    
        # in odd timesteps, choose a direction to expand in
        if timesteps == 0:
            current_neighbors = [k for k in adjacency[source]]
            virtual_source_candidate, current_neighbors = pick_random_elements(current_neighbors,1)
            previous_vs = virtual_source
            
            # infect twice in one direction, always
            infection_pattern, who_infected = infect_nodes(source, virtual_source_candidate, infection_pattern, adjacency, who_infected)
            virtual_source_candidate = virtual_source_candidate[0]
            infection_pattern, who_infected = pass_branch_message(source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
            virtual_source = virtual_source_candidate
            m = 1       # the virtual source is always going to be 1 hop away from the true source
            
        else:
            current_neighbors = [k for k in who_infected[virtual_source]]

            if random.random() < compute_alpha(m,timesteps,max_infection):     # with probability alpha, spread symmetrically (keep the virtual source where it is)
                
                # if there is nowhere for the virtual source to move, keep it where it is
                if len(current_neighbors) < 1:
                    blocked = True
                    break
                    
                # branch once in every direction
                for neighbor in current_neighbors:
                    infection_pattern, who_infected = pass_branch_message(virtual_source, neighbor, infection_pattern, adjacency, max_infection, who_infected)
                if len(current_neighbors) == 1:
                    infection_pattern, who_infected = pass_branch_message(current_neighbors[0], virtual_source, infection_pattern, adjacency, max_infection, who_infected)
            
            else:           # spread asymmetrically
                
                # find a direction to move
                virtual_source_candidate = [previous_vs]
                while virtual_source_candidate[0] == previous_vs:
                    virtual_source_candidate, current_neighbors = pick_random_elements(current_neighbors,1)
                virtual_source_candidate = virtual_source_candidate[0]
                visited_nodes = [virtual_source, virtual_source_candidate]
                previous_vs = virtual_source
                
                # the virtual source moves one more hop away from the true source
                m += 1;
                
                # branch twice in one direction
                infection_pattern, who_infected = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
                infection_pattern, who_infected = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
                
                
            
            
        num_infected = sum(infection_pattern)
        
        # print('infected subgraph', who_infected,'timesteps',timesteps)
        jordan_estimate = estimation.jordan_centrality(who_infected)
        jordan_correct[timesteps] = (jordan_estimate == source)
        if jordan_estimate == source:
            print('num infected', num_infected)    
            # print('source neighbors',who_infected[source])
            
        
        timesteps += 1
        
    return num_infected, infection_pattern, who_infected, jordan_correct
    
def compute_alpha(m,T,d):
    # Compute the probability of keeping the virtual source
    # Inputs
    #       m:          distance from the virtual source to the true source, in hops
    #       T:          time
    #       d:          degree of the d-regular tree
    #
    # Outputs
    #       alpha:      the probability of keeping the virtual source where it is
    
    alpha1 = N(T,d) / (N(T+1,d))
    if m == 1:
        return alpha1
    elif m == 2:
        alpha = d/(d-1) * alpha1 - 1/(d-1)
    else:
        # alpha = alpha1 + compute_alpha(m-1,T,d)/(d-1) - 1/(d-1) 
        alpha = ((1-alpha1)/(d-2))/pow(d-1,m-1) + (alpha1*(d-1)-1)/(d-2)
    return alpha

def N(T,d):
    # Compute the number of nodes at time T in a d-regular graph
    # Inputs
    #       T:          time
    #       d:          degree of the d-regular tree
    #
    # Outputs
    #       n           the number of nodes at time T
    
    n = d / (d-2) * (pow(d-1,T)-1)
    return n
    
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
    return random_elements, neighbors
        
        
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
        adjacency = infect_nodes_infinite_tree(recipient, max_degree-1, adjacency)
    return adjacency
    
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
    
    for neighbor in neighbors:
        leaf = False
        infection_pattern, who_infected =  pass_branch_message(recipient, neighbor, infection_pattern, adjacency, max_infection, who_infected)
            
    if leaf:
        neighbors = [k for k in adjacency[recipient] if infection_pattern[k]==0]
        if len(neighbors) > max_infection:
            neighbors, remains = pick_random_elements(neighbors,max_infection)
        infection_pattern,who_infected = infect_nodes(recipient, neighbors, infection_pattern, adjacency, who_infected)
    return infection_pattern, who_infected
        
def infect_nodes_infinite_tree(node, num_children, adjacency):
    adjacency[node] = adjacency[node] + [i for i in range(len(adjacency),len(adjacency)+num_children)]
    adjacency = adjacency + [[node] for i in range(num_children)]
    return adjacency
    
def infect_nodes(node, children, infection_pattern, adjacency, who_infected):
    # infect_nodes infects the nodes listed in children from 'node'
    # Inputs
    #       node:               source of the infection
    #       children:           array of child ids
    #       infection_pattern:  binary array describing whether each node is already infected or not
    #       adjacency:          adjacency relations for the underlying network
    #       who_infected:       adjacency relations for the infected subgraph
    #
    # Outputs
    #       infection_pattern   (updated)
    #       who_infected        (updated)

    who_infected[node] += children
    for child in children:
        infection_pattern[child] = 1
        who_infected[child] += [node]
    return infection_pattern, who_infected