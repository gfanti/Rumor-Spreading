import utilities
import random

class SpiesInformation(object):
    def __init__(self, timestamps=[], spies = [], active_nodes = [], infectors = [], levels=[]):
        self.timestamps = timestamps
        self.infectors = infectors
        self.active_nodes = active_nodes
        self.spies = spies
        self.levels = levels
        self.active_spies = [] # list of spies spy:(timestamp,level, infector)
        
    def add_nodes(self, node, neighbors, levels, spy_probability = 0.0):
        self.timestamps += [self.timestamps[node]+1 for neighbor in neighbors]
        new_spies = utilities.update_spies_diffusion(neighbors, spy_probability=spy_probability)
        self.spies += new_spies
        # update the infectors
        neighbor_times = [self.active_nodes[node] + 1 for neighbor in neighbors]
        if self.active_nodes[node] < 0:
                self.active_nodes += [self.active_nodes[node] for item in neighbor_times]
        else:
            self.active_nodes += neighbor_times
            for (neighbor,neighbor_time) in zip(neighbors,neighbor_times):
                if neighbor in new_spies:
                    self.active_nodes[neighbor] = -neighbor_time
                    new_spy = Spy(neighbor,self.timestamps[neighbor], abs(levels[neighbor]), levels[neighbor]>0, node)
                    # print(new_spy.node,': ',new_spy.level,new_spy.up_node)
                    self.active_spies.append(new_spy)
                    # self.active_spies[neighbor] = (self.timestamps[neighbor],levels[neighbor],node)
                    # self.levels += [levels[neighbor]]
                    # self.infectors += [node]
                    
class DatasetSpiesInformation(SpiesInformation):
    def __init__(self, num_nodes, source):
        timestamps = [-1 for i in range(num_nodes)]
        timestamps[source] = 0
        active_nodes = [-1 for i in range(num_nodes)]
        active_nodes[source] = 0
        super(DatasetSpiesInformation,self).__init__(timestamps,[],active_nodes,infectors=[],levels=[])
    
    def add_nodes(self, node, neighbors, levels, directions, spy_probability = 0.0):
        for neighbor in neighbors:
            self.timestamps[neighbor] = self.timestamps[node] + 1
        new_spies = utilities.update_spies_diffusion(neighbors, spy_probability=spy_probability)
        self.spies += new_spies
        # update the infectors
        if self.active_nodes[node] < 0:
            for neighbor in neighbors:
                self.active_nodes[neighbor] = self.active_nodes[node]
        else:
            for neighbor in neighbors:
                if neighbor in new_spies:
                    self.active_nodes[neighbor] = -(self.active_nodes[node] + 1)
                    new_spy = Spy(neighbor,self.timestamps[neighbor], levels[neighbor], directions[neighbor], node)
                    self.active_spies.append(new_spy)    
                    # print('Spy ', neighbor, ' was infected by node ',node, ' and is active')
                else:
                    self.active_nodes[neighbor] = self.active_nodes[node] + 1
    
class DatasetDiffusionSpiesInformation(DatasetSpiesInformation):
    def __init__(self, num_nodes, source):
        super(DatasetDiffusionSpiesInformation, self).__init__(num_nodes, source)
        
    def add_nodes(self, node, neighbors, spy_probability = 0.0):
        for neighbor in neighbors:
            self.timestamps[neighbor] = self.timestamps[node] + 1
        new_spies = utilities.update_spies_diffusion(neighbors, spy_probability=spy_probability)
        self.spies += new_spies
        # update the infectors
        if self.active_nodes[node] < 0:
            for neighbor in neighbors:
                self.active_nodes[neighbor] = self.active_nodes[node]
        else:
            for neighbor in neighbors:
                if neighbor in new_spies:
                    self.active_nodes[neighbor] = -(self.active_nodes[node] + 1)
                    new_spy = Spy(neighbor,self.timestamps[neighbor], infector=node)
                    self.active_spies.append(new_spy)    
                    # print('Spy ', neighbor, ' was infected by node ',node, ' and is active')
                else:
                    self.active_nodes[neighbor] = self.active_nodes[node] + 1
    
class Spy(object):
    def __init__(self, node, timestamp=None, level=None, up_node = None, infector=None):
        self.node = node
        self.timestamp = timestamp
        self.level = level
        self.up_node = up_node
        self.infector = infector