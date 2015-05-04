# utilities
import random
from scipy import stats    
import argparse

def parse_args(args):
    '''
    Parses the arguments used to call trial simulations.
    '''
    parser = argparse.ArgumentParser(description = "Run adaptive diffusion over trees.")
    parser.add_argument("-t", '--trials', help = "Number of trials to run (default: 1)", type = int, default = 1)
    parser.add_argument("-w", '--write_results', help = "Write results to file? (default: false)", action = 'store_true')
    parser.add_argument("-d", '--diffusion', help = "Spread using regular diffusion? (default: false)", action = 'store_true')
    parser.add_argument("-s", '--spy_probability', help = "Probability of a node being a spy (default: 0)", type = float, default = 0.0)
    parser.add_argument("-a", '--alt', help = "Use alternative spreading model? (default: false)", action = 'store_true')
    parser.add_argument("--db", nargs = '?', help = "Which database to use(fb=facebook, pg=power grid)", type = str, default = 'none')
    parser.add_argument("-r", '--run', help = "Which run number to save as", type = int, default = 0)
    args = parser.parse_args()
    
    if args.trials:
        trials = int(args.trials)
    else:
        trials = 1
    if args.write_results:
        write_results = bool(args.write_results)
    else:
        write_results = False # Do not write results to file
    if args.diffusion:
        diffusion = bool(args.diffusion)
        if args.spy_probability:
            spy_probability = float(args.spy_probability)
        return {'trials':trials, 'write_results':write_results, 'diffusion':diffusion,'spy_probability':spy_probability}
    if args.alt:
        alt = bool(args.alt)
    else:
        alt = False # Use the alternative method of spreading virtual sources?
    if not (args.db == 'none'):
        database = args.db
        print("The parameters are:\nDataset = ", database,"\nTrials = ",trials,"\nwrite_results = ",write_results,"\nalt = ",alt,"\n")
        if args.run:
            run = args.run
        else:
            run = 0
        return {'trials':trials, 'write_results':write_results, 'database':database, 'run':run} 
        
    print("The parameters are:\nTrials = ",trials,"\nwrite_results = ",write_results,"\nalt = ",alt,"\n")
    return {'trials':trials, 'write_results':write_results, 'alt':alt}

def nCk(n, k):
    '''
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    '''
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
    '''
    Compute the probability of keeping the virtual source
    Arguments:
          m:          distance from the virtual source to the true source, in hops
          T:          time
          d:          degree of the d-regular tree
    
    Outputs:
          alpha:      the probability of keeping the virtual source where it is
    '''
    
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
    '''
    Compute the number of graphs that can appear at time T in a d-regular graph
    Arguments:
          T:          time
          d:          degree of the d-regular tree
    
    Outputs:
          n           the number of nodes at time T
    '''
    
    if d > 2:
        n = float(d) / (d-2) * (pow(d-1,T)-1)
    else:
        n = 1 + 2*T
    return n
    
def N_nodes(T,d):
    '''
    Compute the number of nodes that appear in a graph at time T in a d-regular graph
    '''
    return N(T,d) + 1

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
    
def print_adjacency(adj, true):
    for i in range(len(adj)):
        if len(adj[i]) > 0:
            print(i, ': ', adj[i], ' (', true[i],')')
            
def update_spies_diffusion(candidates, spy_probability = 0.3):
    spies = []
    for candidate in candidates:
        if random.random() < spy_probability:
            spies.append(candidate)
    return spies
        