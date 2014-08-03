import numpy as np
import random 

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
def max_likelihood(who_infected, virtual_source, adjacency):
    paths = get_paths(who_infected, virtual_source)
    for i in range(len(who_infected)):
        # this algorithm only works for leaves in the infected subgraph
        if len(who_infected[i]) == 1:
            likelihoods = compute_graph_likelihood(i, who_infected, adjacency, paths[i], likelihoods)
    indices = [i for i, x in enumerate(likelihoods) if x == max(likelihoods)]
    ml_estimate = indices[random.choice(indices)]
    return ml_estimate
    
def compute_graph_likelihood(source, who_infected, adjacency, vs_path, likelihoods):
    likelihood = 1/len(adjacency[source])
    vs_path.pop(0)
    next_vs = vs_path[0]