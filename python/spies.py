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
        if 3 not in new_spies and max(neighbors) > 3:
            new_spies.append(3)
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
        
        
class Spy(object):
    def __init__(self, node, timestamp=None, level=None, up_node = None, infector=None):
        self.node = node
        self.timestamp = timestamp
        self.level = level
        self.up_node = up_node
        self.infector = infector