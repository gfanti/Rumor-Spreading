# estimation.py
import random
import math
import numpy as np
from numpy.linalg import inv, pinv
import heapq
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from spies import *
import infectionUtils
import estimation

class Estimator(object):
    def __init__(self, adjacency, malicious_nodes, timestamps=None, infectors = None, active_nodes = None, source = None):
        self.adjacency = adjacency
        self.malicious_nodes = malicious_nodes
        self.timestamps = timestamps
        self.infectors = infectors
        self.source = source
        self.graph = nx.Graph()
        
        # Populate the active nodes
        if active_nodes is None:
            self.active_nodes = [1 for i in range(len(adjacency))]
        else:
            self.active_nodes = active_nodes
        # Populate the graph
        for idx in range(len(self.adjacency)):
            edges = self.adjacency[idx]
            for e in edges:
                self.graph.add_edge(idx, e)
        
    def estimate_source(self):
        pass
    
    def prune_graph(self):
        # prunes the graph to include only possible sources
        num_spies = len(self.malicious_nodes)
        G_old = self.graph.copy()
        #Prune edges
        for spy_idx in range(num_spies):
            for e in nx.edges(G_old,[self.malicious_nodes[spy_idx]]):
                if self.infectors[spy_idx] not in e:
                    G_old.remove_edge(e[0],e[1])

        # Keep only the connected component with the first spy
        H = list(nx.connected_component_subgraphs(G_old))
        points = zip(self.timestamps,self.malicious_nodes,self.infectors)
        points = sorted(points)
        times = [item[0] for item in points]
        spies = [item[1] for item in points]
        infectors = [item[2] for item in points]
        members = [spies[0] in comp.nodes() for comp in H]
        idx = members.index(1)
        H = H[idx]

        return H
    
    def get_diameter(self):
        ''' Returns the diameter of the graph'''
        # computes the diameter of the adjacency matrix
        return nx.diameter(self.graph)
        
    def draw_graph(self):
        G = self.graph
        pos = nx.spring_layout(G)
        # print('self active_nodes',self.active_nodes)
        nl = [x for x in G.nodes() if self.active_nodes[x] >= 0]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="#A0CBE2")
        nl = [x for x in G.nodes() if self.active_nodes[x] < 0]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="yellow")
        nl = [x for x in G.nodes() if x in self.malicious_nodes]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="red")
        nl = [0]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="green")
        nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
        labels = {}
        for x in G.nodes(data=True):
            labels[x[0]] = x[0]
        nx.draw_networkx_labels(G, pos, labels)
        plt.show()
        
    def get_spanning_tree(self, node, use_infectors = False):
        ''' Returns a networkx spanning tree of the adjacency matrix
        rooted at node'''
        G = nx.bfs_tree(self.graph, node).to_undirected()
        if not nx.is_connected(G):
            return None
        return G
        
                
class EpflEstimator(Estimator):
    def __init__(self, adjacency, malicious_nodes, timestamps, infectors = None, active_nodes = None, source = None):
        super(EpflEstimator, self).__init__(adjacency, malicious_nodes, timestamps, infectors, active_nodes, source)
        
    def estimate_source(self, use_directions=True, mu_delay=2.0, pd=False):
        '''Estimates the source based on the observed information.
        Arguments:
            use_directions          tells whether to use information about who sent the the message
            mu                      the mean spreading delay
            pd                      return the pd or a real estimate?
            
        Outputs:
            estimate                the index of the node that is estimated'''
            
        # Sums the distance to the unvisited nodes and visited nodes at time_t
        max_likelihood = None
        max_indices = []
        
        if use_directions and (self.infectors is not None):
            G = self.prune_graph()
        else:
            G = self.graph

        points = zip(self.timestamps,self.malicious_nodes)
        pruned_points = [p for p in points if p[1] in G.nodes()]
        times = [item[0] for item in pruned_points]
        spies = [item[1] for item in pruned_points]
        
        num_spies = len(spies)
        
        if num_spies <= 1:
            try:
                return random.choice([node for node in G.nodes() if self.active_nodes[node] >= 0])
            except:
                return random.choice([node for node in range(len(self.adjacency)) if self.active_nodes[node] >= 0])
        
        d = np.array([times[k+1] - times[0] for k in range(num_spies - 1)])
        
        # Cycle through all the nodes in the graph that are active
        for node in G.nodes():
            if (node in spies) or (self.active_nodes[node] < 0):
                continue
            # spanning_tree = nx.bfs_tree(G, node).to_undirected()
            spanning_tree = G # only bc we're on trees, the spanning tree is the graph itself!
            
                    
            mu = np.array([mu_delay*(nx.shortest_path_length(spanning_tree, node, spies[k+1]) - 
                                nx.shortest_path_length(spanning_tree, node, spies[0])) for k in range(num_spies-1)])
                                
            mu.shape = (1,len(mu))
            Lambda_inv, Lambda = self.compute_lambda_inv(node, spies, spanning_tree)
            # subtract distance from nodes that have seen the message already
            d_norm = []
            for idx in range(len(d)):
                d_norm.append(d[idx] - 0.5 * mu[0,idx])
            d_norm = np.transpose(np.array(d_norm))
            likelihood = float(np.dot(np.dot(mu, Lambda_inv), d_norm))
            # likelihood = np.exp(-0.5 * np.dot(np.dot(Lambda_inv, d_norm), d_norm)) / pow(np.linalg.det(Lambda), 0.5)
            
            if (max_likelihood is None) or (max_likelihood < likelihood):
                max_likelihood = likelihood
                max_indices = [node]
            elif (max_likelihood == likelihood):
                max_indices.append(node)
        print('the candidates are ', max_indices)
        # print('the spies are ', spies)
        if max_indices:
            if pd:
                if self.source in max_indices:
                    return 1.0 / len(max_indices)
                else:
                    return 0.0
            else:
                return random.choice(max_indices)
        else:
            print('No valid nodes!')
            if pd:
                return 1.0 / len(self.adjacency)
            else:
                return random.choice([node for node in range(len(self.adjacency)) if self.active_nodes[node] >= 0])
        
    def compute_lambda_inv(self, node, spies, spanning_tree = None):
        num_spies = len(spies)
        Lambda = np.matrix(np.zeros((num_spies-1, num_spies-1)))
        if spanning_tree is None:
            spanning_tree = self.get_spanning_tree(node)
                
        paths = []
        for i in range(num_spies-1):
            source = spies[0]
            destination = spies[i+1]
            path = nx.shortest_path(spanning_tree, source, destination)
            paths.append(path)
        for i in range(num_spies-1):
            for j in range(num_spies-1):
                if i == j:
                    Lambda[i,j] = nx.shortest_path_length(spanning_tree,spies[0],spies[i+1])
                else:
                    p_i = zip(paths[i], paths[i][1:])
                    p_j = zip(paths[j], paths[j][1:])

                    count = 0
                    for e1 in p_i:
                        for e2 in p_j:
                            if (e1[0] == e2[0] and e1[1] == e2[1]) or (e1[1] == e2[0] and e1[0] == e2[1]):
                                count += 1

                    Lambda[i,j] = count
                    Lambda[j,i] = Lambda[i,j]
        Lambda = 0.5**2 * Lambda
        try:
            Lambda_inv = inv(Lambda)
        except:
            # print('matrix was not invertible.')
            Lambda_inv = pinv(Lambda)
        return Lambda_inv, Lambda


class LeiYingEstimator(Estimator):
    def __init__(self, adjacency, malicious_nodes, timestamps, infectors = None, active_nodes = None):
        super(LeiYingEstimator, self).__init__(adjacency, malicious_nodes, timestamps, infectors, active_nodes)
        
    def estimate_source(self, use_directions=True, mu_delay=2):
        '''Estimates the source based on the observed information.
        Arguments:
            use_directions          tells whether to use information about who sent the the message
            
        Outputs:
            estimate                the index of the node that is estimated'''
            
        min_cost = float("inf")
        min_indices = []
        
        # if use_directions and (self.infectors is not None):
            # G = self.prune_graph()
        # else:
        G = self.graph
            
        points = zip(self.timestamps,self.malicious_nodes, self.infectors)
        pruned_points = [p for p in points if p[1] in G.nodes()]
        pruned_points = sorted(pruned_points)
        times = [item[0] for item in pruned_points]
        spies = [item[1] for item in pruned_points]
        infectors = [item[2] for item in pruned_points]
        
        num_spies = len(spies)
        
        if num_spies <= 1:
            try:
                return random.choice([node for node in G.nodes() if self.active_nodes[node] >= 0])
            except:
                return random.choice([node for node in range(len(self.adjacency)) if self.active_nodes[node] >= 0])
        
        # Cycle through all the nodes in the graph that are active
        for node in G.nodes():
            if (node in spies) or (self.active_nodes[node] < 0):
                continue
            
            # Initialize the spreading tree
            spreading_tree = nx.Graph();
            spreading_tree.add_node(node)
            cost = 0
            
            remaining_spies = [spy for spy in spies]
            remaining_times = [timestamp for timestamp in times]
            remaining_infectors = [infector for infector in infectors]
            first_spy = remaining_spies.pop(0)
            first_timestamp = remaining_times.pop(0)
            first_infector = remaining_infectors.pop(0)
            
            # Compute the path from node to the first spy
            modified_graph = G.copy()
            modified_graph.remove_nodes_from(remaining_spies)
            path = nx.shortest_path(modified_graph, node, first_spy)
            spreading_tree.add_path(path)
            # Update the timing info
            spreading_times = {n:(first_timestamp - mu_delay * (len(path)-i)) for (n,i) in zip(path,range(len(path)))}
            # This adds no cost bc we assume all delays are mu
            invalid_candidate = False
            while remaining_spies:
                spy = remaining_spies.pop(0)
                timestamp = remaining_times.pop(0)
                infector = remaining_infectors.pop(0)
                
                # self.draw_progress(spreading_tree,spy,spies)

                min_path_cost, min_path = self.find_modified_path(spy, remaining_spies, spreading_tree, spies, timestamp, spreading_times, mu_delay)
                
                
                # If there are no valid paths, skip to the next candidate
                if not min_path:
                    invalid_candidate = True
                    break
                
                cost += min_path_cost
                spreading_tree.add_path(min_path)
                spreading_times = self.update_timestamps(spreading_times, min_path, spy, timestamp)
            # self.draw_progress(spreading_tree,spy,spies,spreading_tree)
            if invalid_candidate:
                continue
                
            # Update the minimum cost node
            if (cost < min_cost):
                min_cost = cost
                min_indices = [node]
            elif (cost == min_cost):
                min_indices.append(node)
        # print('the Lei Ying candidates are ', min_indices)
        if min_indices:
            return random.choice(min_indices)
        else:
            print('No valid nodes!')
            return random.choice([node for node in range(len(self.adjacency)) if self.active_nodes[node] >= 0])
     
    def find_modified_path(self, spy, remaining_spies, spreading_tree, spies, timestamp, spreading_times, mu_delay):
        # Create a modified graph without the spreading tree and other spies
        spyless_graph = self.graph.copy()
        spyless_graph.remove_nodes_from(remaining_spies)
        
        min_path_cost = float("inf")
        min_path = []
        for n in spreading_tree.nodes():
            if n in spies:
                continue
            modified_graph = spyless_graph.copy()
            invalid_candidates = [i for i in spreading_tree.nodes() if i != n]
            modified_graph.remove_nodes_from(invalid_candidates)
            
            # Compute the modified path to that spy
            try:
                path = nx.shortest_path(modified_graph, n, spy)
            except:
                # There is no valid path from n to the next spy, so we ignore it
                continue
            pathlength = len(path) - 1
            path_cost = pathlength * ((timestamp - spreading_times[n]) / pathlength - mu_delay)**2;
            if path_cost < min_path_cost:
                min_path_cost = path_cost
                min_path = path
        # print('min_cost', min_path_cost, 'min_path', min_path)
        return min_path_cost, min_path
     
    def update_timestamps(self, spreading_times, path, spy, timestamp):
        pathlength = len(path) - 1
        # Update timing information
        mstar = path.pop(0)
        # remove the last item in the list, which is the spy in question, and assign its timestamp
        path.pop() 
        spreading_times[spy] = timestamp
        
        # Assign timestamps for the rest of the items in the list
        hop_length = 1
        for hop in path:
            spreading_times[hop] = spreading_times[mstar] + (hop_length - 1) * ((timestamp - spreading_times[mstar]) / pathlength) 
            hop_length += 1
        return spreading_times
        
    def draw_progress(self,spreading_tree,spy,spies, modified):
        # Initial plot
        f1 = plt.figure()
        G = self.graph
        pos = nx.spring_layout(G)
        # print('self active_nodes',self.active_nodes)
        nl = [x for x in G.nodes()]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="#A0CBE2")
        nl = [x for x in spreading_tree.nodes()]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="yellow")
        nl = [x for x in G.nodes() if x in self.malicious_nodes]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="red")
        nl = [spy]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="orange")
        # nl = [0]
        # nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="green")
        nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
        labels = {}
        for x in G.nodes(data=True):
            labels[x[0]] = x[0]
        nx.draw_networkx_labels(G, pos, labels)
        
        # Modified plot
        f2 = plt.figure()
        G = modified
        pos = nx.spring_layout(G)
        nl = [x for x in G.nodes()]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="#A0CBE2")
        nl = [x for x in spreading_tree.nodes() if x in G.nodes()]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="yellow")
        # nl = [0]
        # nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="green")
        nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
        labels = {}
        for x in G.nodes(data=True):
            labels[x[0]] = x[0]
        nx.draw_networkx_labels(G, pos, labels)
        
        #Show
        plt.show()

    
class FirstSpyEstimator(Estimator):
    def __init__(self, adjacency, malicious_nodes, timestamps, infectors = None, active_nodes = None, source = None):
        super(FirstSpyEstimator, self).__init__(adjacency, malicious_nodes, timestamps, infectors, active_nodes, source)
        
    def estimate_source(self, use_directions=True):
        '''Estimates the source based on the observed information.
        Arguments:
            use_directions          tells whether to use information about who sent the the message
            
        Outputs:
            estimate                the index of the node that is estimated'''
            
        # Picks a random neighbor of the first spy to receive the message
        estimate = random.randint(0, len(self.adjacency)-1)
        
        if use_directions:
            if self.infectors is not None:
                points = zip(self.timestamps,self.malicious_nodes,self.infectors)
                points = sorted(points)
                times = [item[0] for item in points]
                spies = [item[1] for item in points]
                infectors = [item[2] for item in points]
                for (spy,infector) in zip(spies,infectors):
                    if infector not in spies:
                        estimate = infector
                        return estimate
            else:
                points = zip(self.timestamps,self.malicious_nodes)
                points = sorted(points)
                times = [item[0] for item in points]
                spies = [item[1] for item in points]
                if (self.active_nodes is not None):
                    for spy in spies:
                        options = [option for option in self.adjacency[spy] if self.active_nodes[option] >= 0]
                        if options:
                            estimate = random.choice(options)
                            return estimate
        for spy in spies:
            options = [option for option in self.adjacency[spy] if self.active_nodes[option] >= 0]
            if options:
                estimate = random.choice(options)
                return estimate
        return estimate
        
        
class AdaptiveEstimator(Estimator):
    def __init__(self, adjacency, malicious_nodes, active_nodes = None, degrees_rv = None, source=None):
        super(AdaptiveEstimator, self).__init__(adjacency, malicious_nodes, active_nodes, source=source)
        self.degrees_rv = degrees_rv
        
    def estimate_source(self):
        '''Estimates the source based on the observed information.
        Arguments:
            use_directions          tells whether to use information about who sent the the message
            mu                      the mean spreading delay
            
        Outputs:
            estimate                the index of the node that is estimated'''
            
        spies = self.malicious_nodes
        up_spies = [spy for spy in spies if spy.up_node] # level is greater than zero
        if not up_spies:
            print('NO UPS')
            cnt = 0
            paths = []
            while not paths:
                # There is no pivot, so just take the farthest-apart spies in time
                result = self.compute_min_pivot_down(spies[cnt],spies[-1])
                hops, pivot, bad_neighbors = result
                min_pivot = pivot.node
                paths = self.compute_candidate_paths(hops, min_pivot, bad_neighbors)
                cnt += 1
        else:
            up_level = min([spy.level for spy in up_spies])
            up_spy = None
            for spy in up_spies:
                if spy.level == up_level:
                    up_spy = spy
                    break
            #  find the pivot
            # print('Spies are', [spy.node for spy in self.malicious_nodes])
            # print('Up Spy is', up_spy.node, ' with level ', up_spy.level)
            # print('Num nodes: ', len(self.graph.nodes()))
            # print('Num spies: ', len(self.malicious_nodes))
            result = self.compute_min_pivot(up_spy)
            hops, min_pivot, bad_neighbors = result
            paths = self.compute_candidate_paths(hops, min_pivot, bad_neighbors)
        # print ('min pivot',min_pivot)
        # print ('were looking ',hops,' hops away')
        
        # estimate = random.choice(paths)
        
        candidates = self.get_max_likelihoods(paths)
        # FOr regular trees, the true source is always among the ML estimates
        pd = 1.0 / len(candidates)
        # For irregular trees, sometimes not
        if self.source not in candidates:
            pd = 0.0
        

        # print('Positive levels',[level for level in self.levels if level > 0])
        
        return pd
        
    def compute_min_pivot(self, up_spy):
        ''' Given the up-spy and the '''
        spies = self.malicious_nodes
        
        
        pivots = [-1 for i in range(up_spy.level)]
        bad_neighbors = [set() for i in range(up_spy.level)]
        for spy in [ item for item in spies if item.node != up_spy.node]:
            # print('graph nodes:', self.graph.nodes())
            # print('path from x to y:', spy.node, up_spy.node)
            path = nx.shortest_path(self.graph, spy.node, up_spy.node)
            # print('path:', path)
            if up_spy.infector in path:
                # print('nodes',spy.node,'end',up_spy.node,'path', path)
                pivot = self.compute_pivot(up_spy, spy, path)
                # print('found pivot',pivot.node)
                # print('level',pivot.level, pivots)
                if pivot.level < up_spy.level:
                    pivots[pivot.level] = pivot.node
                    # pivot_path = path[:pivot.level - spy.level + 1]
                    # print('path:',path)
                    # print('pivot:',pivot.node)
                    # print('levels:',pivot.level - spy.level - 1,pivot.level - spy.level+1)
                    # print('pivot: ', pivot, 'pivot path', pivot_path)
                    bad_neighbors[pivot.level].add(path[pivot.level - spy.level - 1])
                    bad_neighbors[pivot.level].add(path[pivot.level - spy.level + 1])
                
        # print('all pivots',pivots)
        hops = 0
        for pivot,neighbor_list in zip(pivots, bad_neighbors):
            if pivot >= 0:
                return (hops, pivot, neighbor_list) 
            hops += 1
        print('There were no viable pivots.',len(self.malicious_nodes))
        return None
        
    def compute_min_pivot_down(self, down_spy1, down_spy2):
        ''' Given the farthest down_spies, guess the pivot '''
        spies = self.malicious_nodes
        
        
        bad_neighbors = set()
        path = nx.shortest_path(self.graph, down_spy1.node, down_spy2.node)
        pivot = self.compute_pivot(down_spy2, down_spy1, path, True) #happens to be the same for (up-down) and (down-down)
        for spy in [ item for item in spies if item.node != pivot.node]:
            pivot_path = nx.shortest_path(self.graph, spy.node, pivot.node)
            # print('path:',pivot_path)
            # print('pivot:',pivot.node)
            # print('to omit:',pivot_path[-2])
            bad_neighbors.add(pivot_path[-2])
                
        hops = pivot.level
        return (hops, pivot, bad_neighbors) 
            
        
    def compute_candidate_paths(self, hops, min_pivot, bad_neighbors):
    
        # Get the feasible candidate paths
        # print('hops',hops,'min_pivot',min_pivot)
        all_paths = nx.single_source_shortest_path_length(self.graph, min_pivot, cutoff=hops)
        # print('all paths before',all_paths)
        all_paths = [(k,v) for k,v in all_paths.items() if v == hops]
        
        # print('all paths',all_paths)
        # print('bad neighbors',bad_neighbors)
        
        # Prune out those passing through bad neighbors
        feasible_paths = []
        for candidate in all_paths:
            path = nx.shortest_path(self.graph, candidate[0], min_pivot)
            if path[-2] not in bad_neighbors:
                feasible_paths += [path]
        
        return feasible_paths
        
    def compute_pivot(self, up_spy, spy, path, both_down = False):
        Pl = len(path) - 1  # path length between spies
        dt = up_spy.timestamp - spy.timestamp
        h1 = int(0.5 * (Pl - dt))
        if both_down:
            dm = up_spy.level - spy.level
            h2 = int(0.5 * (dt + dm))
            h3 = int(0.5 * (Pl - dm))
            print(h1,h2,h3)
            if h2 < 1:
                pivot = Spy(path[h1 + abs(h2)],level=up_spy.level + h3,up_node=True)
                return pivot
        
        if spy.level < 0:
            print('MEGA PROBLEM IN COMPUTE_PIVOT')
        pivot = Spy(path[h1],level=spy.level + h1,up_node=True)
        
        return pivot
        
    def get_max_likelihoods(self, paths):
    
        if len(self.degrees_rv.xk) == 1:
            return [path[0] for path in paths]
    
        likelihoods = []
        for path in paths:
            degrees = list(self.graph.degree(path[1:]).values())
            degrees = self.degrees_rv.draw_values(1) + degrees
            degrees = [1.0/(degrees[i]-1) if i != 0 else 1.0/(degrees[i]) for i in range(len(degrees)) ]
            likelihoods.append(np.prod(degrees))
        # Get the ML candidates
        candidates = []
        for path,likelihood in zip(paths, likelihoods):
            if likelihood == max(likelihoods):
                candidates.append(path[0])
            
        return candidates
            
        
class DatasetAdaptiveEstimator(AdaptiveEstimator):
    def __init__(self, adjacency, who_infected, malicious_nodes, active_nodes = None, source=None, max_infection=3):
        super(DatasetAdaptiveEstimator, self).__init__(who_infected, malicious_nodes, active_nodes,source=source)
        self.contact_adjacency = adjacency
        self.max_infection = max_infection
        
        self.contact_graph = nx.Graph()
        # Populate the graph
        for idx in range(len(adjacency)):
            self.contact_graph.add_node(idx)
            edges = self.contact_adjacency[idx]
            for e in edges:
                self.contact_graph.add_edge(idx, e)
        
        
    def get_max_likelihoods(self, paths):
        ''' Need to overwrite the method in AdaptiveEstimator that uses degree distributions'''
        # return [path[0] for path in paths]
        # paths = [[] for i in range(len(who_infected))]
        # paths[virtual_source].append(virtual_source)
        # paths = get_paths(paths, who_infected, adjacency, virtual_source, virtual_source)
        
        # maximum length from leaf to virtual source
        likelihoods = []
        for path in paths:
            likelihood = self.compute_graph_likelihood(path)
            likelihoods.append(likelihood)
            
        # Get the ML candidates
        candidates = []
        max_likelihood = max(likelihoods)
        for path,likelihood in zip(paths, likelihoods):
            if likelihood == max_likelihood:
                candidates.append(path[0])
            
        return candidates
        
    def compute_graph_likelihood(self, path):
        if len(path) == 1:
            print('ERROR: The source is a spy!', self.source, path)
            return [self.source]
    
        # print('THe path is ',path)
        # print('THe source is ',self.source)
        path.pop()
        nodes = self.contact_graph.nodes()
        num_nodes = self.contact_graph.number_of_nodes()
        # print('There are ',num_nodes,'nodes here')
        # print('And the other graph has ',self.graph.number_of_nodes(),'nodes')
        new_infection_pattern = [0 for i in nodes]
        new_infection_pattern[self.source] = 1
        new_who_infected = [[] for i in nodes]

        # first element is just the source itself
        current_vs = path.pop(0)
        path_source = current_vs
        # log likelihood of the 1st passage to the VS
        likelihood = math.log(1.0/self.contact_graph.degree(self.source))
        if not path:
            path.append(path_source)
            # print('path',path)
            return likelihood
        # get the new vs
        current_vs = path.pop(0)
        new_infection_pattern, new_who_infected, tmp = infectionUtils.infect_nodes(self.source, [current_vs], new_infection_pattern, new_who_infected)
        
        # infect the neighbors of the new vs
        infected = [i for i in self.adjacency[current_vs]]
        infected.remove(path_source)
        
        new_infection_pattern, new_who_infected, tmp = infectionUtils.infect_nodes(current_vs, infected, new_infection_pattern, new_who_infected)
        likelihood += estimation.infect_set_likelihood(infected, self.contact_adjacency[current_vs], new_infection_pattern, self.max_infection)
        
        while path:
            new_infection_pattern, new_who_infected, likelihood = estimation.pass_branch_message_likelihood(current_vs, path[0], new_infection_pattern, 
                                                                    self.contact_adjacency, self.max_infection, new_who_infected, self.adjacency, likelihood)
            current_vs = path.pop(0)
            
        path.append(path_source)
        
        # print('path',path)
        return likelihood
        
        