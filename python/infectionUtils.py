#infectionUtils.py useful functions for infection
import utilities
import random
from scipy import stats
import numpy as np
import time
import networkx as nx
import estimation_spies     # TODO: get rid of this and function estimate_source_spies!


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
        adjacency = infect_nodes_infinite_tree(recipient, max_degree-1, adjacency)
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
        degrees,who_infected, known_degrees = infect_nodes_randtree(recipient, neighbors, degrees, degrees_rv, who_infected, known_degrees)
    return degrees, who_infected, known_degrees
    
        
def infect_nodes_up_down(node, num_to_infect, boundary, levels, degrees, degrees_rv, who_infected):
    '''
    infect_nodes infects the nodes listed in children from 'node' over a random tree
    Arguments:
        node:               source of the infection
        children:           array of child ids
        degrees:            array of infected nodes' degrees
        degrees_rv:         random variable describing the degree distribution
        who_infected:       adjacency relations for the infected subgraph
    
    Outputs:
        infection_pattern   (updated)
        who_infected        (updated)
    '''
    num_infected = len(who_infected)
    children = [num_infected + i for i in range(num_to_infect)]
    who_infected[node] += children
    for child in children:
        who_infected.append([node])
    degrees += degrees_rv.rvs(size=len(children)).tolist()
    if levels[node] is None:
        levels += [1 for child in children]
    elif levels[node] > 0:
        up_node = random.choice(children)
        levels += [-(levels[node]-1) for child in children]
        levels[up_node] = -levels[up_node] + 2
    elif levels[node] < 0:
        levels += [levels[node]+1 for child in children]
    boundary += [child for child in children if levels[child] != 0]
    return degrees, who_infected, boundary, levels, children
    
def infect_nodes_infinite_tree(node, num_children, adjacency):
    '''
    Infect nodes in an infinite tree given the adjacency matrix
    
    Arguments:
        node: The infecting node
        num_children: The number of children to infect
        adjacency: The adjacency matrix of the infected subgraph so far.
    
    Outputs:
        The updated adjacency matrix.
    '''

    adjacency[node] = adjacency[node] + [i for i in range(len(adjacency),len(adjacency)+num_children)]
    adjacency = adjacency + [[node] for i in range(num_children)]
    return adjacency

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
        infection_pattern,who_infected, dist_from_source = infect_nodes(recipient, neighbors, infection_pattern, who_infected, dist_from_source)
    return infection_pattern, who_infected, dist_from_source
 

def infect_nodes(node, children, infection_pattern, who_infected, dist_from_source = []):
    '''
    Infects the nodes listed in children
    Inputs
          node:               source of the infection
          children:           array of child ids
          infection_pattern:  binary array describing whether each node is already infected or not
          adjacency:          adjacency relations for the underlying network
          who_infected:       adjacency relations for the infected subgraph
    
    Outputs
          infection_pattern   (updated)
          who_infected        (updated)
    '''
    
    who_infected[node] += children
    for child in children:
        infection_pattern[child] = 1
        who_infected[child] += [node]
        if dist_from_source:
            dist_from_source[child] = dist_from_source[node] + 1
    return infection_pattern, who_infected, dist_from_source
 
def infect_nodes_randtree(node, children, degrees, degrees_rv, who_infected, known_degrees=[]):
    '''
    infect_nodes infects the nodes listed in children from 'node' over a random tree
    Arguments:
        node:               source of the infection
        children:           array of child ids
        degrees:            array of infected nodes' degrees
        degrees_rv:         random variable describing the degree distribution
        who_infected:       adjacency relations for the infected subgraph
    
    Outputs:
        infection_pattern   (updated)
        who_infected        (updated)
    '''
    
    who_infected[node] += children
    for child in children:
        who_infected.append([node])
    if not known_degrees:
        # degrees += degrees_rv.rvs(size=len(children)).tolist()
        degrees += degrees_rv.draw_values(len(children))    
    elif len(known_degrees) >= len(children):
        degrees += known_degrees[0:len(children)]
        known_degrees[0:len(children)] = []
    else:
        slack = len(children) - len(known_degrees)
        degrees += known_degrees
        known_degrees = []
        # degrees += degrees_rv.rvs(size=slack).tolist()
        degrees += degrees_rv.draw_values(slack)
        
    return degrees, who_infected, known_degrees
    
    
    
def estimate_source_spies(max_time, est_times, source, who_infected, num_infected, timestamps, spies,
                          active_nodes, active_spies, infectors=None, diffusion=False, q=None):
    
        
    ml_correct, spy_correct, lei_correct = [],[],[]
    hop_distances, spy_hop_distances,lei_hop_distances = [],[],[]
    tot_num_infected = []
    
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
        
        if diffusion:   # diffusion
            # estimator = estimation_spies.EpflEstimator(adjacency, spies, spies_timestamps, active_nodes=spy_active_nodes)
            # spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, spies, spies_timestamps, active_nodes=spy_active_nodes)
            # estimator = estimation_spies.EpflEstimator(adjacency, spies, spies_timestamps, active_nodes=active_nodes)
            # spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, spies, spies_timestamps, active_nodes=active_nodes)
            # estimator = estimation_spies.EpflEstimator(adjacency, active_spies, spies_timestamps, active_nodes=current_active_nodes)
            # spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, active_spies, spies_timestamps, active_nodes=current_active_nodes)
            estimator = estimation_spies.EpflEstimator(adjacency, reached_spies, spies_timestamps, active_nodes=current_active_nodes)
            spy_estimator = estimation_spies.FirstSpyEstimator(adjacency, reached_spies, spies_timestamps, active_nodes=current_active_nodes)
            lei_estimator = estimation_spies.LeiYingEstimator(adjacency, reached_spies, spies_timestamps, infectors=infectors, active_nodes=current_active_nodes)
        else:
            estimator = estimation_spies.AdaptiveEstimator(adjacency, reached_spies, spies_timestamps, active_nodes=current_active_nodes)
        
        # for spy in reached_spies:
            # print('spy',spy,' has time ',timestamps[spy],' and is ',nx.shortest_path_length(estimator.graph,0,spy),' hops away')
        
        if active_spies:
            start = time.time()
            if diffusion:
                mean = 1.0 / q
                ml_estimate = estimator.estimate_source(mu_delay=mean)
            else:
                ml_estimate = estimator.estimate_source()
            # print('Surrounded? ',not any([current_active_nodes[n]>0 for n in estimator.graph if len(estimator.adjacency[n])==1]))
            print('ml est', ml_estimate)
            end = time.time()
            print('ml elapsed:',end-start)
            # estimator.draw_graph()
            
            if diffusion:
                spy_estimate = spy_estimator.estimate_source()
                start = time.time()
                print('spy elapsed:',start-end)
                # lei_estimate = lei_estimator.estimate_source(mu_delay=mean)
                lei_estimate = 0
                end = time.time()
                print('lei elapsed:',end - start)
            
            
        else:
            # choose a random node
            ml_estimate = random.randint(0,num_infected - 1)
            if diffusion:
                spy_estimate = ml_estimate
                lei_estimate = ml_estimate
        hop_distance = nx.shortest_path_length(estimator.graph, source, ml_estimate)
        print('True source: ', source, ' EPFL estimate: ', ml_estimate)
        if diffusion:
            print('True source: ', source, ' spy estimate: ', spy_estimate)
            print('True source: ', source, ' Lei Ying estimate: ', lei_estimate)
            spy_hop_distance = nx.shortest_path_length(estimator.graph, source, spy_estimate)
            lei_hop_distance = nx.shortest_path_length(estimator.graph, source, lei_estimate)
            spy_correct.append(spy_estimate == source)
            spy_hop_distances.append(spy_hop_distance)
            lei_correct.append(lei_estimate == source)
            lei_hop_distances.append(lei_hop_distance)

        
        ml_correct.append(ml_estimate == source)
        hop_distances.append(hop_distance)
        tot_num_infected.append(num_infected)
        
    print('Num spies are: ', len(spies), ' out of ', num_infected)
    print('\n\n')
    if diffusion:
        infection_details = (tot_num_infected, who_infected, hop_distances, spy_hop_distances, lei_hop_distances)
        results = (ml_correct, spy_correct, lei_correct)
    else:
        infection_details = (tot_num_infected, who_infected, hop_distances)
        results = (ml_correct)
    return infection_details, results
    