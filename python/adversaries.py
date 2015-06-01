import estimation_spies
import time
import networkx as nx
import random
import numpy as np

class Adversary(object):
    def __init__(self):
        pass
        
        
class UpDownAdversary(Adversary):
    
    def __init__(self, source, spies_info, who_infected, degrees_rv):
        super(UpDownAdversary, self).__init__()
        
        self.source = source
        self.spies_info = spies_info
        self.who_infected = who_infected
        self.num_infected = len(who_infected)
        self.degrees_rv = degrees_rv
        
    def get_estimates(self, max_time, est_times=None):
        '''Estimate the source  '''
        if est_times is None:
            est_times = [max_time]
            
        ml_correct, hop_distances = [],[]
        
        adjacency = [set(item) for item in self.who_infected]
        nodes = range(self.num_infected)
        spies_info = self.spies_info

        # for node in nodes:
            # print('Node ',node,': ',timestamps[node])
        for est_time in est_times:
            
            # reached_spy_points = [(spy,level,infector) for spy,level,infector in 
                # zip(spies_info.active_spies, spies_info.levels, spies_info.infectors) if spies_info.timestamps[spy] <= est_time]
            # reached_spies = [item[0] for item in reached_spy_points]
            # reached_levels = [item[1] for item in reached_spy_points]
            # reached_infectors = [item[2] for item in reached_spy_points]
            spies_timestamps = [spy.timestamp for spy in spies_info.active_spies if spy.timestamp <= est_time]
            
            # reached_spies = [spy for spy in spies_info.active_spies]
            reached_spies = [spy for spy in spies_info.active_spies if spy.timestamp <= est_time]
            # print('rr2',reached_spies)
            
            # print('reached spies:',reached_spies)
            # print('reached levels:',reached_levels)
            # print('infected by:',reached_infectors)
            
                
            
            # current_active_nodes = [item if ((abs(item) <= est_time) and (n not in active_spies)) else abs(item) for (item,n) in zip(active_nodes,nodes)]
            if len(reached_spies) == 0:
                current_active_nodes = [1 if n not in spies_info.spies else -1 for n in nodes]
            else:
                current_active_nodes = spies_info.active_nodes
            print('est time',est_time)
            
            
            estimator = estimation_spies.AdaptiveEstimator(self.who_infected, reached_spies,
                                                           active_nodes=current_active_nodes,
                                                           degrees_rv=self.degrees_rv,
                                                           source=self.source)
            

            if len(reached_spies) > 1:
                start = time.time()
                
                pd_ml = estimator.estimate_source()
                print('ml est', pd_ml)
                end = time.time()
                print('ml elapsed:',end-start)
                # estimator.draw_graph()
                
            else:
                # print('\n Not enough spies!!!!\n')
                print('reached spies',[item.node for item in reached_spies])
                # choose a random node
                # ml_estimate = random.randint(0,self.num_infected - 1)
                if len(reached_spies) == 0:
                    pd_ml = 1.0 / self.num_infected
                elif len(reached_spies) == 1 and reached_spies[0].up_node:
                    pd_ml = 1.0 / ((self.degrees_rv.mean() - 1)**(reached_spies[0].level-1))
                else:
                    pd_ml = 1.0 / (self.num_infected - 1)
                
            # hop_distance = nx.shortest_path_length(estimator.graph, self.source, ml_estimate)
            hop_distance = 0
            print('True source: ', self.source, ' ML adaptive estimate pd: ', pd_ml)

            # ml_correct.append(ml_estimate == self.source)
            ml_correct.append(pd_ml)
            hop_distances.append(hop_distance)
            
        results = (ml_correct, hop_distances)
        return results
            