from spies import *

class Infector(object):
    def __init__(self):
        
        self.who_infected = None
        self.num_infected = 0
        self.source = None
        
class UpDownInfector(Infector):
    ''' Up-down algorithms over irregular trees (this
    is the version we used in our simulations for spies)'''
    
    def __init__(self, spy_probability, degrees_rv):
        super(UpDownInfector, self).__init__()
        
        self.spy_probability = spy_probability
        self.degrees_rv = degrees_rv
        
        self.boundary = []
        self.new_boundary = []
        self.levels = [None]
        
        self.spies_info = None
        
        

    def infect(self, source, max_time):
        
        '''infect_nodes_adaptive_irregular_tree runs the spreading model over an irregular
        tree.
        
        Arguments:
            source: The ID of the source node (usually 0)
            max_time: The maximum amount of timesteps to run the algorithm
            max_infection: Defines the (nominal degree-1) that is assumed in order to choose alpha
            degrees_rv: A random variable that describes the graph degree
            
        Returns:
            infection_details: A list of the characteristics of the infection:
                - tot_num_infected: Total number of infected nodes at each timestep
                - infection_pattern: Which nodes got infected in which order
                - who_infected: Adjacency matrix of infected subgraph
            ml_correct: The ML estimate of the true source'''
    
    
        # ML estimate
        tot_num_infected = [0 for i in range(max_time)]
        
        self.source = source
        
        self.who_infected = [[]]
        self.num_infected = 1
        self.spies_info = SpiesInformation([0],[],[source],[],[])
        self.degrees = self.degrees_rv.rvs(size=1).tolist()
        self.levels = [None]
        # infect twice in one direction, always
        new_neighbors = self.infect_nodes_up_down(source, 1)
        self.spies_info.add_nodes(source, new_neighbors, self.levels, self.spy_probability)
        timesteps = 1
        
        while timesteps < max_time:
        
                
            while self.boundary:
                # print(self.boundary)
                node = self.boundary.pop(0)
                # print('node',node)
                # branch
                new_neighbors = self.infect_nodes_up_down(node, self.degrees[node]-1)
                self.spies_info.add_nodes(node, new_neighbors, self.levels, self.spy_probability)
                
            # Swap the boundary arrays
            self.boundary, self.new_boundary = self.new_boundary, self.boundary
            
            # The estimation for spies comes at the end, whereas we can do snapshot estimation on the fly
            
            tot_num_infected[timesteps] = self.num_infected
            timesteps += 1
            # print('num infected', num_infected)
        
        
        infection_details = (self.who_infected, self.degrees, tot_num_infected)
            
        
        return infection_details, self.spies_info
        
        
    def infect_nodes_up_down(self, node, num_to_infect):
        '''
        infect_nodes infects the nodes listed in children from 'node' over a random tree
        Arguments:
            node:               source of the infection
            children:           array of child ids
            degrees:            array of infected nodes' degrees
            degrees_rv:         random variable describing the degree distribution
            who_infected:       adjacency relations for the infected subgraph
        
        Outputs:
            infection_pattern   (updated)
            who_infected        (updated)
        '''
        children = [self.num_infected + i for i in range(num_to_infect)]
        
        #update the adjacency matrix
        self.who_infected[node] += children
        for child in children:
            self.who_infected.append([node])
        self.num_infected = len(self.who_infected)
        
        #update the degrees
        self.degrees += self.degrees_rv.rvs(size=len(children)).tolist()
        
        #update the metadata
        if self.levels[node] is None:
            self.levels += [1 for child in children]
        elif self.levels[node] > 0:
            up_node = random.choice(children)
            self.levels += [-(self.levels[node]-1) for child in children]
            self.levels[up_node] = -self.levels[up_node] + 2
        elif self.levels[node] < 0:
            self.levels += [self.levels[node]+1 for child in children]
        self.new_boundary += [child for child in children if self.levels[child] != 0]
        return children