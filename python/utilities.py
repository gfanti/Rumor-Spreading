# utilities
import random
from scipy import stats


def nCk(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0
        
def compute_alpha(m,T,d):
    # Compute the probability of keeping the virtual source
    # Inputs
    #       m:          distance from the virtual source to the true source, in hops
    #       T:          time
    #       d:          degree of the d-regular tree
    #
    # Outputs
    #       alpha:      the probability of keeping the virtual source where it is
    
    alpha1 = float(N(T,d)) / (N(T+1,d))
    if m == 1:
        return alpha1
    else:
        # alpha = alpha1 + compute_alpha(m-1,T,d)/(d-1) - 1/(d-1) 
        if d > 2:
            alpha = (float(1-alpha1)/(d-2))/pow(d-1,m-1) + float(alpha1*(d-1)-1)/(d-2)
        else:
            alpha = float(T-m) / T
    return alpha

def N(T,d):
    # Compute the number of graphs that can appear at time T in a d-regular graph
    # Inputs
    #       T:          time
    #       d:          degree of the d-regular tree
    #
    # Outputs
    #       n           the number of nodes at time T
    
    if d > 2:
        n = float(d) / (d-2) * (pow(d-1,T)-1)
    else:
        n = 1 + 2*T
    return n
    
def N_nodes(T,d):
    # Compute the number of nodes that appear in a graph at time T in a d-regular graph
    return N(T,d) + 1
    
def infect_nodes_infinite_tree(node, num_children, adjacency):
    adjacency[node] = adjacency[node] + [i for i in range(len(adjacency),len(adjacency)+num_children)]
    adjacency = adjacency + [[node] for i in range(num_children)]
    return adjacency
    
def infect_nodes(node, children, infection_pattern, who_infected):
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
    
def infect_nodes_randtree(node, children, degrees, degrees_rv, who_infected, known_degrees=[]):
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
        who_infected.append([node])
    if not known_degrees:
        degrees += degrees_rv.rvs(size=len(children)).tolist()
    elif len(known_degrees) >= len(children):
        degrees += known_degrees[0:len(children)]
        known_degrees[0:len(children)] = []
    else:
        slack = len(children) - len(known_degrees)
        degrees += known_degrees
        known_degrees = []
        degrees += degrees_rv.rvs(size=slack).tolist()
        
    return degrees, who_infected, known_degrees