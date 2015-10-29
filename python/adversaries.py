import estimation_spies
import estimation_snapshot
import time
import networkx as nx
import random
import numpy as np

class Adversary(object):
    def __init__(self, source, who_infected):
        self.source = source
        self.who_infected = who_infected
        self.num_infected = len(who_infected)


''' ------------------SNAPSHOT ADVERSARIES------------------------'''
class SnapshotAdversary(Adversary):
    def __init__(self, source, who_infected = []):
        super(SnapshotAdversary, self).__init__(source, who_infected)
        self.ml_correct = []

    def update_data(self, who_infected, degrees):
        self.who_infected = who_infected
        self.degrees = degrees

class RandTreeSnapshotAdversary(SnapshotAdversary):
    def __init__(self, source, degrees_rv, who_infected=[]):
        super(RandTreeSnapshotAdversary, self).__init__(source, who_infected)
        self.degrees_rv = degrees_rv

class ADSnapshotAdversary(RandTreeSnapshotAdversary):
    def __init__(self, source, degrees_rv, who_infected=[], d_o = 1000):
        super(ADSnapshotAdversary, self).__init__(source, degrees_rv, who_infected)
        self.d_o = d_o

    def get_estimates(self, virtual_source):

        estimator = estimation_snapshot.ADEstimatorRandTree(self.who_infected, self.source, self.degrees, self.degrees_rv, self.d_o)
        estimate = estimator.estimate_source(virtual_source)

        # self.ml_correct.append(estimate == self.source)
        self.ml_correct.append(estimate)

class PAADSnapshotAdversary(RandTreeSnapshotAdversary):
    def __init__(self, source, degrees_rv, who_infected=[]):
        super(PAADSnapshotAdversary, self).__init__(source, degrees_rv, who_infected)

    def update_data(self, who_infected, degrees, src_neighbors):
        super(PAADSnapshotAdversary, self).update_data(who_infected, degrees)
        self.src_neighbors = src_neighbors

    def get_estimates(self, virtual_source):

        estimator = estimation_snapshot.PAADEstimatorRandTree(self.who_infected, self.source, self.degrees, self.degrees_rv, num_hops_pa = 1)
        estimate = estimator.estimate_source(virtual_source, self.src_neighbors)

        # self.ml_correct.append(estimate == self.source)
        self.ml_correct.append(estimate)



''' ------------------SPY-BASED ADVERSARIES------------------------'''
class SpyAdversary(Adversary):
    def __init__(self, source, who_infected, spies_info):
        super(SpyAdversary, self).__init__(source, who_infected)
        self.spies_info = spies_info
        
class DiffusionSpiesAdversary(SpyAdversary):
    
    def __init__(self, source, spies_info, who_infected):
        super(DiffusionSpiesAdversary, self).__init__(source, who_infected, spies_info)
    
    def get_estimates(self, max_time, est_times=None):
        '''Estimate the source  '''
        if est_times is None:
            est_times = [max_time]
            
        ml_correct, spy_correct, hop_distances = [],[],[]
        
        spies_info = self.spies_info
        nodes = range(self.num_infected)

        for est_time in est_times:
            
            spies_timestamps = [spy.timestamp for spy in spies_info.active_spies if spy.timestamp <= est_time]
            reached_spies = [spy for spy in spies_info.active_spies if spy.timestamp <= est_time]
            
            if len(reached_spies) == 0:
                current_active_nodes = [1 if n not in spies_info.spies else -1 for n in nodes]
            else:
                current_active_nodes = spies_info.active_nodes
            print('est time',est_time)

            malicious_nodes, timestamps, infectors = [],[],[]
            for spy in reached_spies:
                malicious_nodes.append(spy.node)
                timestamps.append(spy.timestamp)
                infectors.append(spy.infector)
            
            epfl_estimator = estimation_spies.EpflEstimator(self.who_infected, malicious_nodes,
                                                            timestamps, infectors, spies_info.active_nodes)

            spy_estimator = estimation_spies.FirstSpyEstimator(self.who_infected, malicious_nodes,
                                                               timestamps, infectors, spies_info.active_nodes)
                                                            
            if len(reached_spies) > 1:
                # start = time.time()
                
                # pd_ml = epfl_estimator.estimate_source(pd=True)
                pd_ml = 0.0
                spy_est = spy_estimator.estimate_source()
                print('ml est', pd_ml, 'true source', self.source)
                pd_spy = spy_est == self.source
                print('\nspy est', spy_est, 'true source', self.source)
                # end = time.time()
                # print('ml elapsed:',end-start)
                # estimator.draw_graph()
                
            else:
                print('reached spies',[item.node for item in reached_spies])
                # choose a random node
                pd_ml = 1.0 / self.num_infected
                pd_spy = pd_ml
                
            hop_distance = 0
            print(' ML adaptive estimate pd: ', pd_ml)

            # ml_correct.append(ml_estimate == self.source)
            ml_correct.append(pd_ml)
            spy_correct.append(pd_spy)
            hop_distances.append(hop_distance)
            
            
        results = (ml_correct, spy_correct, hop_distances)
        return results
        
class UpDownAdversary(SpyAdversary):
    
    def __init__(self, source, spies_info, who_infected, degrees_rv):
        super(UpDownAdversary, self).__init__(source, who_infected, spies_info)
        
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
            
class DatasetUpDownAdversary(SpyAdversary):
    
    def __init__(self, source, spies_info, who_infected, adjacency, max_infection = 3):
        super(DatasetUpDownAdversary, self).__init__(source, who_infected, spies_info)
        
        self.adjacency = adjacency
        self.max_infection = max_infection
        
    def get_estimates(self, max_time, est_times=None):
        '''Estimate the source  '''
        if est_times is None:
            est_times = [max_time]
            
        ml_correct, hop_distances = [],[]
        
        nodes = range(self.num_infected)
        spies_info = self.spies_info

        for est_time in est_times:
            
            spies_timestamps = [spy.timestamp for spy in spies_info.active_spies if spy.timestamp <= est_time]
            
            reached_spies = [spy for spy in spies_info.active_spies if spy.timestamp <= est_time]
            if len(reached_spies) == 0:
                current_active_nodes = [1 if n not in spies_info.spies else -1 for n in nodes]
            else:
                current_active_nodes = spies_info.active_nodes
            print('est time',est_time)
            
            
            estimator = estimation_spies.DatasetAdaptiveEstimator(self.adjacency, self.who_infected, reached_spies,
                                                           active_nodes=current_active_nodes,
                                                           source=self.source,
                                                           max_infection=self.max_infection)

            if len(reached_spies) > 1:
                # start = time.time()
                
                pd_ml = estimator.estimate_source()
                print('ml est', pd_ml)
                # end = time.time()
                # print('ml elapsed:',end-start)
                # estimator.draw_graph()
                
            else:
                print('reached spies',[item.node for item in reached_spies])
                # choose a random node
                if len(reached_spies) == 0:
                    pd_ml = 1.0 / self.num_infected
                elif len(reached_spies) == 1: # and reached_spies[0].up_node:
                    spy = reached_spies[0]
                    if spy.level == 1:
                        pd_ml = 1.0
                    else:
                        all_paths = nx.single_source_shortest_path(estimator.graph,spy.infector,spy.level-1).items()
                        candidates = 0
                        for path in all_paths:
                            if spy.node not in path:
                                candidates += 1
                        pd_ml = 1.0 / candidates
                else:
                    pd_ml = 1.0 / (self.num_infected - 1)
                
            hop_distance = 0
            print(' ML adaptive estimate pd: ', pd_ml)

            # ml_correct.append(ml_estimate == self.source)
            ml_correct.append(pd_ml)
            hop_distances.append(hop_distance)
            
        results = (ml_correct, hop_distances)
        return results