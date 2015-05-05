# estimation.py
import random
import math
import numpy as np
from numpy.linalg import inv, pinv
import heapq
import networkx

class Estimator(object):
    def __init__(self, adjacency, malicious_nodes, timestamps, active_nodes = None):
        self.adjacency = adjacency
        self.malicious_nodes = malicious_nodes
        self.timestamps = timestamps
        if active_nodes is None:
            self.active_nodes = [1 for i in range(len(adjacency))]
        else:
            self.active_nodes = active_nodes
        
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
        max_dist = 0
        for node in range(len(self.adjacency)):
            distances = self.get_distances(node)
            if max(distances) > max_dist:
                max_dist = max(distances)
        return max_dist
        
    def get_spanning_tree(self, node):
        ''' Returns a networkx spanning tree of the adjacency matrix
        rooted at node'''
        num_nodes = len(self.adjacency)
        # sp_adjacency = [set() for i in range(num_nodes)]
        G = networkx.Graph()
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
                
class OptimalEstimator(Estimator):
    def estimate_source(self):
        # Sums the distance to the unvisited nodes and visited nodes at time_t
        max_likelihood = None
        max_indices = []
        num_spies = len(self.malicious_nodes)
        print('spies', self.malicious_nodes)
        print('timestamps', self.timestamps)
        # d = np.diff(self.timestamps)
        d = np.array([self.timestamps[k+1] - self.timestamps[0] for k in range(num_spies - 1)])
        # First compute the paths between spy 1 and the rest
                            
        for node in range(len(self.adjacency)):
            if (node in self.malicious_nodes) or (self.active_nodes[node] == -1):
                continue
            distances = self.get_distances(node)
            # 2 is the mean delay if a message gets forwarded
            mu = np.array([2*(distances[self.malicious_nodes[k+1]] - distances[self.malicious_nodes[0]]) for k in range(num_spies-1)])
            mu.shape = (1,len(mu))
            # print('timestamps are ', self.timestamps)
            # print('mu is ', mu, 'd is ',d)
            Lambda_inv = self.compute_lambda_inv(node)
            # subtract distance from nodes that have seen the message already
            d_norm = np.array([item_d - 0.5*item_mu for (item_d, item_mu) in zip(d, mu)])
            d_norm = np.transpose(d_norm)
            likelihood = float(np.dot(np.dot(mu, Lambda_inv), d_norm))
            print('Node ', node,': likelihood is ', likelihood)
            if (max_likelihood is None) or (max_likelihood < likelihood):
                max_likelihood = likelihood
                max_indices = [node]
            elif (max_likelihood == likelihood):
                max_indices.append(node)
            
        print('the candidates are ', max_indices)
        print('with likelihood ', max_likelihood)
        # print('the spies are ', self.malicious_nodes)
        estimate = random.choice(max_indices)
        print(self.adjacency[estimate])
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
            path = networkx.shortest_path(spanning_tree, source, destination)
            path.pop(0)
            # print('path is ', path)
            # print('original adjacency is ', self.adjacency)
            # print('dijstra result from', source, ' to ', destination, ' gives ', path)
            paths.append(set(path))
        for i in range(num_spies-1):
            for j in range(num_spies-1):
                if i == j:
                    # Lambda[i,j] = spy_distances[i+1]
                    Lambda[i,j] = networkx.shortest_path_length(spanning_tree,self.malicious_nodes[0],self.malicious_nodes[i+1])
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
            # print('matrix was not invertible.')
            # return max_index
            Lambda_inv = pinv(Lambda)
        return Lambda_inv
       
class FirstSpyEstimator(Estimator):
    def estimate_source(self):
        # Picks a random neighbor of the first spy to receive the message
        estimate = random.randint(0, len(self.adjacency)-1)
        for spy in self.malicious_nodes:
            options = [option for option in self.adjacency[spy] if self.active_nodes[option] == 1]
            # print('options', options)
            if options:
                estimate = random.choice(options)
                break
        return estimate
       
class SumDistanceEstimator(Estimator):
    def estimate_source(self, time_t):
        # Sums the distance to the unvisited nodes and visited nodes at time_t
        max_sum_distance = -100000
        max_index = -1
        for node in range(len(self.adjacency)):
            sum_distance = 0
            distances = self.get_distances(node)
            # subtract distance from nodes that have seen the message already
            sum_distance -= sum([distances[i] for i in range(len(self.malicious_nodes)) if self.timestamps[i] <= time_t])
            # add the distance to nodes that did not see the message yet
            sum_distance += sum([distances[i] for i in range(len(self.malicious_nodes)) if self.timestamps[i] > time_t])
            if max_sum_distance < sum_distance:
                max_index = node
                max_sum_distance = sum_distance
        return max_index
       
class EntropyEstimator(Estimator):
    def estimate_source(self):
        # Sums the distance to the unvisited nodes and visited nodes at time_t
        max_entropy = -1
        min_index = -1
        for node in range(len(self.adjacency)):
            distances = self.get_distances(node)
            # subtract distance from nodes that have seen the message already
            distribution = [distances[i]/self.timestamps[i] for i in range(len(self.malicious_nodes))]
            distribution = [item / sum(distribution) for item in distribution if item > 0]
            entropy = self.compute_entropy(distribution)
            if max_entropy < entropy:
                max_index = node
                max_entropy = entropy
        return max_index
        
    def compute_entropy(self, dist):
        # computes the entropy of distribution dist
        return sum([-p*math.log(p) for p in dist])
        
class JordanEstimator(Estimator):
    
    def estimate_source(self):
        # Computes the jordan estimate of the source
        jordan_dist = -1
        jordan_center = -1
        for node in range(len(self.adjacency)):
            dist = self.compute_jordan_centrality(node)
            if (jordan_dist == -1) or (dist < jordan_dist):
                jordan_dist = dist
                jordan_center = node
        print("The best estimate is ", jordan_center, " with a centrality of ", jordan_dist)
        return jordan_center
    
    def compute_jordan_centrality(self, source):
        # Computes the jordan centrality of a single node 'source', as viewed by the malicious nodes
        
        # compute the distances to all malicious nodes
        distances = self.get_distances(source)
        
        # Extract the longest distance to any one of the malicious nodes
        jordan_dist = max([distances[j] for j in self.malicious_nodes])
        return jordan_dist
        
    