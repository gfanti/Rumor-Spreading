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
    def __init__(self, adjacency, malicious_nodes, timestamps, active_nodes = None):
        self.adjacency = adjacency
        self.malicious_nodes = malicious_nodes
        self.timestamps = timestamps
        self.graph = nx.Graph()
        
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
        
    def get_distances(self, source):
        ''' Returns a vector with the distance from source to every node in the
        graph.'''
        visited, queue = set(), [source]
        distances = [0 for i in range(len(self.adjacency))]
        counter = 0
        # while (0 in [distances[j] for j in self.malicious_nodes]):
        while (0 in distances):
            vertex = queue.pop(0)
            if vertex not in visited:
                visited.add(vertex)
                queue.extend(self.adjacency[vertex] - visited)
                vertex_dist = distances[vertex]
                for i in (self.adjacency[vertex] - visited):
                    distances[i] = vertex_dist + 1
            if not queue:
                break
        return distances
        
    def get_diameter(self):
        ''' Returns the diameter of the graph'''
        # computes the diameter of the adjacency matrix
        return nx.diameter(self.graph)
        
    def get_spanning_tree(self, node):
        ''' Returns a networkx spanning tree of the adjacency matrix
        rooted at node'''
        num_nodes = len(self.adjacency)
        # sp_adjacency = [set() for i in range(num_nodes)]
        G = nx.Graph()
        for vertex in range(len(self.adjacency)):
            G.add_node(vertex)
        nodes = set([i for i in range(num_nodes)])
        visited, queue = set(), [node]
        while queue and nodes:
            vertex = queue.pop(0)
            if vertex not in visited:
                nodes.remove(vertex)
                visited.add(vertex)
                queue.extend(self.adjacency[vertex] - visited)
                # Fix the adjacency matrix here!
                for i in (self.adjacency[vertex] - visited):
                    # sp_adjacency[vertex].add(i)
                    # sp_adjacency[i].add(vertex)
                    G.add_edge(vertex, i)
        # return sp_adjacency
        return G
        
    def draw_graph(self):
        G = self.graph
        pos = nx.spring_layout(G)
        nl = [x for x in G.nodes() if x not in self.malicious_nodes]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="#A0CBE2")
        nl = [x for x in G.nodes() if x == 0]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="green")
        nl = [x for x in G.nodes() if x in self.malicious_nodes]
        nx.draw_networkx_nodes(G,pos,nodelist=nl,node_color="red")
        nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
        labels = {}
        for x in G.nodes(data=True):
            labels[x[0]] = x[0]
        nx.draw_networkx_labels(G, pos, labels)
        plt.show()
                
class OptimalEstimator(Estimator):
    def estimate_source(self):
        # Sums the distance to the unvisited nodes and visited nodes at time_t
        max_likelihood = None
        max_indices = []
        num_spies = len(self.malicious_nodes)
        # print('spies', self.malicious_nodes)
        # print('timestamps', self.timestamps)
        # d = np.diff(self.timestamps)
        d = np.array([self.timestamps[k+1] - self.timestamps[0] for k in range(num_spies - 1)])
        # First compute the paths between spy 1 and the rest
                            
        # print('adjacency',self.adjacency)
        for node in range(len(self.adjacency)):
            if (node in self.malicious_nodes) or (self.active_nodes[node] == -1):
                continue
            # distances = self.get_distances(node)
            # 2 is the mean delay if a message gets forwarded
            mu = np.array([2.0*(nx.shortest_path_length(self.graph,node, self.malicious_nodes[k+1]) - 
                                nx.shortest_path_length(self.graph,node, self.malicious_nodes[0])) for k in range(num_spies-1)])
            mu.shape = (1,len(mu))
            # print('timestamps are ', self.timestamps)
            # print('mu is ', mu, 'd is ',d)
            Lambda_inv = self.compute_lambda_inv(node)
            # subtract distance from nodes that have seen the message already
            # print('d is',d, 'mu is ',mu)
            # print('diff is', [item_d - 0.5*item_mu for (item_d, item_mu) in zip(d, mu)])
            # d_norm = np.array([item_d - 0.5*item_mu for (item_d, item_mu) in zip(d, mu)])
            d_norm = []
            for idx in range(len(d)):
                d_norm.append(d[idx] - 0.5 * mu[0,idx])
            d_norm = np.transpose(np.array(d_norm))
            # print('d_norm', d_norm, 'mu', mu)
            # print('lambda inv', Lambda_inv)
            likelihood = float(np.dot(np.dot(mu, Lambda_inv), d_norm))
            # print('Node ', node,': likelihood is ', likelihood)
            if (max_likelihood is None) or (max_likelihood < likelihood):
                max_likelihood = likelihood
                max_indices = [node]
            elif (max_likelihood == likelihood):
                max_indices.append(node)
            
        # print('the candidates are ', max_indices)
        # print('with likelihood ', max_likelihood)
        # if 0 not in max_indices:
            # self.draw_graph()
        # print('the spies are ', self.malicious_nodes)
        estimate = random.choice(max_indices)
        # print(self.adjacency[estimate])
        return estimate
        
    def compute_lambda_inv(self, node):
        num_spies = len(self.malicious_nodes)
        Lambda = np.matrix(np.zeros((num_spies-1, num_spies-1)))
        spanning_tree = self.get_spanning_tree(node)
        # distances = self.get_distances(self.malicious_nodes[0])
        # spy_distances = [distances[i] for i in self.malicious_nodes]
                
        paths = []
        for i in range(num_spies-1):
            source = self.malicious_nodes[0]
            destination = self.malicious_nodes[i+1]
            # path = self.dijkstra(source, destination, spanning_tree)
            path = nx.shortest_path(spanning_tree, source, destination)
            path.pop(0)
            # print('path is ', path)
            # print('original adjacency is ', self.adjacency)
            # print('dijstra result from', source, ' to ', destination, ' gives ', path)
            paths.append(set(path))
        for i in range(num_spies-1):
            for j in range(num_spies-1):
                if i == j:
                    # Lambda[i,j] = spy_distances[i+1]
                    Lambda[i,j] = nx.shortest_path_length(spanning_tree,self.malicious_nodes[0],self.malicious_nodes[i+1])
                else:
                    Lambda[i,j] = len(paths[i].intersection(paths[j]))
                    Lambda[j,i] = Lambda[i,j]
        # print('Adjacency: ', self.adjacency)
        # print('Spies: ', self.malicious_nodes)
        # print('Spy_times: ', self.timestamps)
        # print('spy distances', spy_distances)
        # print('Lambda: ', Lambda)
        try:
            Lambda_inv = inv(Lambda)
        except:
            print('matrix was not invertible.')
            # self.draw_graph()
            # return max_index
            Lambda_inv = pinv(Lambda)
        return Lambda_inv
       
class FirstSpyEstimator(Estimator):
    def estimate_source(self):
        # Picks a random neighbor of the first spy to receive the message
        estimate = random.randint(0, len(self.adjacency)-1)
        for spy in [x for (y,x) in sorted(zip(self.timestamps,self.malicious_nodes))]:
            options = [option for option in self.adjacency[spy] if self.active_nodes[option] == 1]
            # print('options', options)
            if options:
                estimate = random.choice(options)
                break
        return estimate