# node infection

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
    timesteps = 1
    
    # initially the virtual source and the true source are the same
    virtual_source = source
    virtual_source_candidate = virtual_source
    
    infection_pattern = [0]*num_nodes;
    infection_pattern[source] = 1;
    who_infected = [[] for i in range(num_nodes)]
    num_infected = 0
    # infected_adjacency = [[] for i in range(num_nodes)]
    
    blocked = False
    
    while timesteps <= max_time:
        # print('time', timesteps)
    
        # in odd timesteps, choose a direction to expand in
        if timesteps == 1:
            current_neighbors = [k for k in adjacency[source] if infection_pattern[k]==0]
            virtual_source_candidate, current_neighbors = pick_random_elements(current_neighbors,1)
            
            # infect twice in one direction, always
            infection_pattern, who_infected = infect_nodes(source, virtual_source_candidate, infection_pattern, adjacency, who_infected)
            virtual_source_candidate = virtual_source_candidate[0]
            infection_pattern, who_infected = pass_branch_message(source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
            virtual_source = virtual_source_candidate
            m = 1       # the virtual source is always going to be 1 hop away from the true source
            
        else:
            current_neighbors = [k for k in who_infected[virtual_source]]

            if random.random() < compute_alpha(m,timesteps,max_infection):     # with probability alpha, we should expand in all directions (keep the virtual source where it is)
                if len(current_neighbors) < 1:
                    # print('Blocked. Total neighbors: ',adjacency[node])
                    # print('Available neighbors: ',current_neighbors)
                    blocked = True
                    break
                    
                # branch once in every direction
                # print('spreading symmetrically')
                for neighbor in current_neighbors:
                    infection_pattern, who_infected = pass_branch_message(virtual_source, neighbor, infection_pattern, adjacency, max_infection, who_infected)
                if len(current_neighbors) == 1:
                    infection_pattern, who_infected = pass_branch_message(current_neighbors[0], virtual_source, infection_pattern, adjacency, max_infection, who_infected)
            else:
                
                # find a direction to move
                previous_vs_candidate = virtual_source_candidate
                while previous_vs_candidate == virtual_source_candidate:
                    virtual_source_candidate, current_neighbors = pick_random_elements(current_neighbors,1)
                virtual_source_candidate = virtual_source_candidate[0]
                visited_nodes = [virtual_source, virtual_source_candidate]
                
                # the virtual source moves one more hop away from the true source
                m += 1;
                
                # print('spreading asymmetrically')
                # branch twice in one direction
                infection_pattern, who_infected = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
                infection_pattern, who_infected = pass_branch_message(virtual_source, virtual_source_candidate, infection_pattern, adjacency, max_infection, who_infected)
                
                
            
            
        num_infected = sum(infection_pattern)
        # print('num infected', num_infected)
        timesteps += 1
        
    return num_infected, infection_pattern
    
def compute_alpha(m,T,d):
    alpha1 = N(T,d) / (N(T+1,d))
    if m == 1:
        return alpha1
    elif m == 2:
        alpha = d/(d-1) * alpha1 - 1/(d-1)
    else:
        alpha = alpha1 + compute_alpha(m-1,T,d)/(d-1) - 1/(d-1) 
    return alpha

def N(T,d):
    n = 0
    for i in range(T):
        n += d * pow(d-1,i)
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
            # print(infection_pattern,node)
            # print('neighbors',adjacency[node])
            current_neighbors = [k for k in adjacency[node] if infection_pattern[k]==0]
            # print('valid neighbors: ',current_neighbors)
            if len(current_neighbors) < 2:
                # print('Blocked. Total neighbors: ',adjacency[node])
                # print('Available neighbors: ',current_neighbors)
                blocked = True
                break
            if infection_pattern[node] == 1:
                up_elements, current_neighbors = pick_random_elements(current_neighbors,1)
                for item in up_elements:
                    infection_pattern[item] = 1
                down_elements, current_neighbors = pick_random_elements(current_neighbors,1)
                for item in down_elements:
                    infection_pattern[item] = -1
                new_infecting_nodes = new_infecting_nodes + up_elements + down_elements
            elif infection_pattern[node] == -1:
                down_elements, current_neighbors = pick_random_elements(current_neighbors,2)
                for item in down_elements:
                    infection_pattern[item] = -1
                new_infecting_nodes = new_infecting_nodes + down_elements
            else:
                print('ERROR')
                
            num_infected += 2
        if blocked:
            break
            
        # swap the arrays
        old_infecting_nodes, new_infecting_nodes = new_infecting_nodes, old_infecting_nodes
    
    return num_infected, infection_pattern
        
def pick_random_elements(neighbors,num_neighbors):
    random_elements = []
    for i in range(num_neighbors):
        random_element = random.choice(neighbors)
        neighbors.remove(random_element)
        random_elements.append(random_element)
    return random_elements, neighbors
        
        
def pass_branch_message_infinite_tree(source, recipient, adjacency, max_degree):
    leaf = True
    neighbors = adjacency[recipient]
    # print('neighbors are ',neighbors)
    for neighbor in neighbors:
        if not neighbor == source:
            leaf = False
            adjacency =  pass_branch_message_infinite_tree(recipient, neighbor, adjacency, max_degree)
            
    if leaf:
        adjacency = infect_nodes_infinite_tree(recipient, max_degree-1, adjacency)
    return adjacency
    
def pass_branch_message(source, recipient, infection_pattern, adjacency, max_infection, who_infected):
    leaf = True
    
    # print('adjacency [recipient]',adjacency[recipient])
    
    # pass only to already infected nodes that are not the source
    # neighbors = [k for k in adjacency[recipient] if (not k in visited_nodes) and (infection_pattern[k] == 1)]
    # print('neighbors',neighbors)
    
    # pass to the neighbors who are your neighbors in the true infection graph
    neighbors = [k for k in who_infected[recipient] if (not k == source)]
    
    # print(neighbors,len(adjacency[recipient]))
    
    for neighbor in neighbors:
        leaf = False
        # visited_nodes.append(neighbor)
        infection_pattern, who_infected =  pass_branch_message(recipient, neighbor, infection_pattern, adjacency, max_infection, who_infected)
        # print(neighbors, sum(infection_pattern),source)
        # input("Press Enter to continue...")
            
    if leaf:
        neighbors = [k for k in adjacency[recipient] if infection_pattern[k]==0]
        if len(neighbors) > max_infection:
            neighbors, remains = pick_random_elements(neighbors,max_infection)
        # print('infecting ',len(neighbors),' neighbors')
        infection_pattern,who_infected = infect_nodes(recipient, neighbors, infection_pattern, adjacency, who_infected)
        # print('num nodes infected = ',sum(infection_pattern))
    return infection_pattern, who_infected
        
def infect_nodes_infinite_tree(node, num_children, adjacency):
    adjacency[node] = adjacency[node] + [i for i in range(len(adjacency),len(adjacency)+num_children)]
    # print('new item',[[node] for i in range(num_children)])
    adjacency = adjacency + [[node] for i in range(num_children)]
    
    # print('new adjacency', adjacency)
    return adjacency
    
def infect_nodes(node, children, infection_pattern, adjacency, who_infected):
    who_infected[node] += children
    # children = array of child ids
    for child in children:
        infection_pattern[child] = 1
        who_infected[child] += [node]
    return infection_pattern, who_infected