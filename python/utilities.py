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
    parser.add_argument("-q", '--delay_parameter', help = "Probability of passing the message in a timestep (default: 0)", type = float, default = 0.5)
    parser.add_argument("-a", '--alt', help = "Use alternative spreading model? (default: false)", action = 'store_true')
    parser.add_argument("--db", nargs = '?', help = "Which database to use(fb=facebook, pg=power grid)", type = str, default = 'none')
    parser.add_argument("-r", '--run', help = "Which run number to save as", type = int, default = 0)
    args = parser.parse_args()

    # Num trials
    if args.trials:
        trials = int(args.trials)
    else:
        trials = 1
    # Write results
    if args.write_results:
        write_results = bool(args.write_results)
    else:
        write_results = False # Do not write results to file
    # Run number
    if args.run:
            run = args.run
    else:
        run = 0
    # Spy probability
    if args.spy_probability:
        spy_probability = float(args.spy_probability)
        spies = True
    else: 
        spy_probability = 0.0
        spies = False
    if args.diffusion:
        diffusion = bool(args.diffusion)
        # return {'trials':trials, 'write_results':write_results, 'diffusion':diffusion,
                # 'spy_probability':spy_probability,'run':run, 'delay_parameter':delay_parameter}
    else:
        diffusion = False
    if args.delay_parameter:
        delay_parameter = float(args.delay_parameter)
    else:
        delay_parameter = 0.5
    if args.alt:
        alt = bool(args.alt)
    else:
        alt = False # Use the alternative method of spreading virtual sources?
    if not (args.db == 'none'):
        database = args.db
        print("The parameters are:\nDataset = ", database,"\nTrials = ",trials,"\nwrite_results = ",write_results,"\nalt = ",alt,"\n")
        
        # return {'trials':trials, 'write_results':write_results, 'database':database, 'run':run, 'spy_probability':spy_probability} 
    else:
        database = None
        
    print("The parameters are:\nTrials = ",trials,"\nwrite_results = ",write_results,"\nalt = ",alt,"\n")
    return {'trials':trials, 'write_results':write_results, 'alt':alt,'spy_probability':spy_probability,'run':run, 
            'diffusion':diffusion, 'spies':spies, 'delay_parameter':delay_parameter, 'database':database}

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
    

    if d > 2:
        return float(pow((d-1),(T-m+1))-1)/(pow(d-1,T+1)-1)
    elif d == 2:
        return float(T-m+1)/(T+1)

    # alpha1 = float(N(T,d)) / (N(T+1,d))
    # if m == 1:
    #     return alpha1
    # else:
    #     # alpha = alpha1 + compute_alpha(m-1,T,d)/(d-1) - 1/(d-1) 
    #     if d > 2:
    #         alpha = (float(1-alpha1)/(d-2))/pow(d-1,m-1) + float(alpha1*(d-1)-1)/(d-2)
    #     else:
    #         alpha = float(T-m) / T
    # return alpha

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
    
# def N_nodes(T,d):
#     '''
#     Compute the number of nodes that appear in a graph at time T in a d-regular graph
#     '''
#     return N(T,d) + 1

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
        