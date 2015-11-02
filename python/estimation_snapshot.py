import numpy as np
import random 
import utilities
import networkx as nx
import matplotlib.pyplot as plt
import pygraphviz
from networkx.algorithms.traversal.depth_first_search import dfs_tree


class Estimator(object):
    def __init__(self, who_infected, source):
        self.who_infected = who_infected
        self.source = source
        self.graph = nx.Graph()
        
        # Populate the graph
        for idx in range(len(self.who_infected)):
            edges = self.who_infected[idx]
            for e in edges:
                self.graph.add_edge(idx, e)

    def get_diameter(self):
        ''' Returns the diameter of the graph'''
        # computes the diameter of the adjacency matrix
        return nx.diameter(self.graph)


class EstimatorRandTree(Estimator):
    def __init__(self, who_infected, source, degrees, degrees_rv):
        super(EstimatorRandTree, self).__init__(who_infected, source)
        self.degrees = degrees
        self.degrees_rv = degrees_rv

    # def estimate_source(self, virtual_source, src_neighbors):
    #     '''Estimates the source based on the observed information.
    #     Arguments:
    #         virtual_source          the virtual_source of the infected graph
    #         src_neighbors           the uninfected neighbors of the source
            
    #     Outputs:
    #         estimate                the index of the node that is estimated'''

    #     self.src_neighbors = src_neighbors

    #     # initialize the messages vector to likelihood 1
    #     p = 1.0
    #     self.messages = [p]*len(self.who_infected)

    #     self.pass_messages(virtual_source, virtual_source)
    #     messages = [self.messages[i] if len(self.who_infected[i])==1 else 0 for i in range(len(self.messages))]

    #     # finding the likelihood of the ML estimate
    #     max_message = max(messages)
    #     # finding the indices of most likely nodes
    #     max_message_ind = [i for i, j in enumerate(messages) if j == max_message]
    #     # ml_estimate = max_message_ind[random.randrange(0,len(max_message_ind),1)]
    #     # print('choices', max_message_ind)
    #     ml_estimate = random.choice(max_message_ind)
    #     return ml_estimate



class ADEstimatorRandTree(EstimatorRandTree):
    def __init__(self, who_infected, source, degrees, degrees_rv, d_o = 1000):
        super(ADEstimatorRandTree, self).__init__(who_infected, source, degrees, degrees_rv)
        self.d_o = d_o

    def estimate_source(self, virtual_source):
        '''Estimates the source based on the observed information.
        Arguments:
            virtual_source          the virtual_source of the infected graph
            src_neighbors           the uninfected neighbors of the source
            
        Outputs:
            estimate                the index of the node that is estimated'''

        # initialize the messages vector to likelihood 1
        p = 1.0
        self.messages = [p]*len(self.who_infected)

        self.pass_messages(virtual_source, virtual_source)
        messages = [self.messages[i] if (i != virtual_source) else 0 for i in range(len(self.messages))]
        # print("messages", messages)
        # # Plot with pygraphviz
        # A = nx.to_agraph(self.graph)
        # A.layout('dot', args='-Nfontsize=10 -Nwidth=".2" -Nheight=".2" -Nmargin=0 -Gfontsize=8')
        # A.draw('current.png')
        # print("degrees", self.degrees)

        # finding the likelihood of the ML estimate
        max_message = max(messages)
        # finding the indices of most likely nodes
        max_message_ind = [i for i, j in enumerate(messages) if j == max_message]
        # ml_estimate = max_message_ind[random.randrange(0,len(max_message_ind),1)]
        
        if self.source in max_message_ind:
            ml_estimate = 1.0 / len(max_message_ind)
        else:
            ml_estimate = 0.0
            
        # ml_estimate = random.choice(max_message_ind)

        return ml_estimate

    def pass_messages(self, calling, called):
    # ML over irregular infinite trees
    # def ml_message_passing_irregular_trees(d, depth, messages, degrees, who_infected, called, calling): 
        # d:            d+1 is the assumed regular degree (aka d_o)
        # depth:        distance from the virtual source
        # messages:     the likelihood messages stored at each node
        # degrees:      degree of each node in the irregular graph
        # who_infected: adjacency matrix of the infected subgraph
        # called:       which node received the message
        # calling:      which node passed the message

        if len(self.who_infected[called]) != 1: #if not a leaf
            for i in self.who_infected[called]:
                if i == calling:
                    continue
                if called == calling:
                    # (Prob. of choosing virtual source) x (Prob. infecting infected nodes)
                    self.messages[i] = self.messages[calling] * (1.0 / (self.degrees[i])) * self.d_o
                else:
                    self.messages[i] = self.messages[called] * (self.d_o - 1) * \
                        (float(self.degrees[called])/(self.degrees[called]-1)) * (1.0/(self.degrees[i]))
                self.pass_messages(called, i) 

    

class PAADEstimatorRandTree(EstimatorRandTree):
    def __init__(self, who_infected, source, degrees, degrees_rv, num_hops_pa = 1):
        super(PAADEstimatorRandTree, self).__init__(who_infected, source, degrees, degrees_rv)
        self.num_hops_pa = num_hops_pa

    def estimate_source(self, virtual_source, local_neighborhood):
        '''Estimates the source based on the observed information.
        Arguments:
            virtual_source          the virtual_source of the infected graph
            local_neighborhood           the uninfected neighbors of the source
            
        Outputs:
            estimate                the index of the node that is estimated'''

        self.local_neighborhood = local_neighborhood

        # initialize the messages vector to likelihood 1
        p = 1.0
        self.messages = [p]*len(self.who_infected)

        self.pass_messages(virtual_source, virtual_source)
        messages = [self.messages[i] if len(self.who_infected[i])==1 else 0 for i in range(len(self.messages))]
        # print("messages", messages)
        # # Plot with pygraphviz
        # A = nx.to_agraph(self.local_neighborhood)
        # A.layout('dot', args='-Nfontsize=10 -Nwidth=".2" -Nheight=".2" -Nmargin=0 -Gfontsize=8')
        # A.draw('current.png')

        # finding the likelihood of the ML estimate
        max_message = max(messages)
        # finding the indices of most likely nodes
        max_message_ind = [i for i, j in enumerate(messages) if j == max_message]
        # ml_estimate = random.choice(max_message_ind)
        # ml_estimate = 1.0 / len(max_message_ind)
        if self.source in max_message_ind:
            ml_estimate = 1.0 / len(max_message_ind)
        else:
            ml_estimate = 0.0
        return ml_estimate

    def pass_messages(self, calling, called, prev_prob = []):

        local_neighborhood = self.local_neighborhood

        correction_factor = 1

        if len(self.who_infected[called]) != 1: #if not a leaf
            # correction_factor = float(prev_prob[1]) / (prev_prob[1]-degrees[called])
            for i in self.who_infected[called]:
                if i != calling:
                    # print("getting called deg")
                    called_deg = float(local_neighborhood.get_hop_degrees(i, called, self.num_hops_pa-1)) 
                    # prob = [called_deg, sum([self.degrees[k]-1 for k in self.who_infected[i]])]
                    # print("computing prob")
                    prob = [called_deg, sum([local_neighborhood.get_hop_degrees(i, k, self.num_hops_pa-1) \
                            for k in local_neighborhood.neighbors(i) ]) ]
                    if called != calling:
                        correction_factor = float(prev_prob[1]) / (prev_prob[1] - local_neighborhood.get_hop_degrees(called, i, self.num_hops_pa-1) )
                    
                    self.messages[i] = self.messages[called] * correction_factor * prob[0] / prob[1]
                    # print("message at ",i," is ", messages[i])
                    self.pass_messages(called, i, prob) 
