# estimation.py
import random
import math
import numpy as np
from numpy.linalg import inv, pinv
import heapq
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

class Estimator(object):
    def __init__(self, adjacency, malicious_nodes, timestamps, infectors = None, active_nodes = None):
        self.adjacency = adjacency
        self.malicious_nodes = malicious_nodes
        self.timestamps = timestamps
        self.infectors = infectors
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
        
                
class OptimalEstimator(Estimator):
    def __init__(self, adjacency, malicious_nodes, timestamps, infectors = None, active_nodes = None):
        super(OptimalEstimator, self).__init__(adjacency, malicious_nodes, timestamps, infectors, active_nodes)

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
        
    def estimate_source(self, use_directions=True):
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
        
        for node in G.nodes():
            # print('Node ',node,' : active',self.active_nodes[node])
            if (node in spies) or (self.active_nodes[node] < 0):
                continue
            spanning_tree = nx.bfs_tree(G, node).to_undirected()
            
                    
            mu = np.array([2.0*(nx.shortest_path_length(spanning_tree, node, spies[k+1]) - 
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
            return random.choice(max_indices)
        else:
            print('No valid nodes!')
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

       
class FirstSpyEstimator(Estimator):
    def __init__(self, adjacency, malicious_nodes, timestamps, infectors = None, active_nodes = None):
        super(FirstSpyEstimator, self).__init__(adjacency, malicious_nodes, timestamps, infectors, active_nodes)
        
    def estimate_source(self, use_directions=True):
        # Picks a random neighbor of the first spy to receive the message
        # print(self.timestamps)
        # print('timestamps 0',self.timestamps[0],self.malicious_nodes[0],'adj:',self.adjacency[self.malicious_nodes[0]])
        estimate = random.randint(0, len(self.adjacency)-1)
        
        if use_directions:
            if self.infectors is not None:
                # print('infectors', self.infectors)
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
        print('no spies!')
        return estimate