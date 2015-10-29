from spies import *
from infectionUtils import *

# class Infector(object):
#     def __init__(self):
        
#         self.who_infected = None
#         self.num_infected = 0
#         self.source = None

class DatasetDiffusionInfector(Infector):
    def __init__(self, adjacency, spy_probability, max_infection, q=0.5):
        super(DatasetDiffusionInfector, self).__init__()
        
        self.spy_probability = spy_probability
        self.adjacency=adjacency
        self.spies_info = None
        self.num_nodes = len(adjacency)
        self.max_infection = max_infection
        self.q = q
        
    def infect(self, source, max_time):
        # ML estimate
        
        tot_num_infected = [0 for i in range(max_time)]
        self.spies_info = DatasetDiffusionSpiesInformation(self.num_nodes,source)
        self.source = source
        
        self.who_infected = [[] for i in range(self.num_nodes)]
        self.infection_pattern = [0 for i in range(self.num_nodes)]
        self.infection_pattern[source] = 1
        self.num_infected = 1
        
        timesteps = 0
        self.boundary = [source]
        self.new_boundary = []
        
        while timesteps < max_time:
            while self.boundary:
                node = self.boundary.pop(0)
                
                current_neighbors = [k for k in self.adjacency[node] if self.infection_pattern[k]==0]
                # Each neighbor spread with probability q 
                to_infect = utilities.update_spies_diffusion(current_neighbors,self.q)
                if len(to_infect) < len(current_neighbors):
                    self.new_boundary.append(node)
                self.infect_nodes_diffusion(source, to_infect)
                
            self.boundary, self.new_boundary = self.new_boundary, self.boundary
            tot_num_infected[timesteps] = self.num_infected
            timesteps += 1

        print('Reached ', self.num_infected,' nodes out of ',self.num_nodes,'. Adjacency has ',len(self.adjacency))
        print('Reached ', len(self.spies_info.spies),' spies, ', len(self.spies_info.active_spies), ' of which are active.')
        # if len(self.spies_info.active_spies)>0:
            # print('Source:',source, 'First spy info:',self.spies_info.active_spies[0].level,self.spies_info.active_spies[0].up_node)
        infection_details = (self.who_infected, tot_num_infected)
            
        return infection_details
        
    def infect_nodes_diffusion(self, node, children):
        '''Infects the children nodes.'''
        #update the adjacency matrix
        self.who_infected[node] += children
        for child in children:
            self.who_infected[child].append(node)
            self.infection_pattern[child] = 1
        self.num_infected = sum(self.infection_pattern)    
        
        self.new_boundary += children
        
                
        #Update the spies??
        self.spies_info.add_nodes(node, children, self.spy_probability)
        

        
class DatasetUpDownInfector(Infector):
    def __init__(self, adjacency, spy_probability, max_infection):
        super(DatasetUpDownInfector, self).__init__()
        
        self.spy_probability = spy_probability
        self.adjacency=adjacency
        self.spies_info = None
        self.num_nodes = len(adjacency)
        self.levels = [-1 for i in range(self.num_nodes)]
        self.directions = [None for i in range(self.num_nodes)]
        self.max_infection = max_infection
        
    def infect(self, source, max_time):
        # ML estimate
        
        tot_num_infected = [0 for i in range(max_time)]
        
        self.source = source
        
        self.who_infected = [[] for i in range(self.num_nodes)]
        self.infection_pattern = [0 for i in range(self.num_nodes)]
        self.infection_pattern[source] = 1
        self.num_infected = 1
        self.spies_info = DatasetSpiesInformation(self.num_nodes,source)
        self.levels[source] = 0
        self.directions[source] = True
        
        virtual_source = source
        timesteps = 0
        
        while timesteps < max_time:
            # print('time ',timesteps)
            # print('infection ',self.who_infected)
            if timesteps == 0:
                current_neighbors = [k for k in self.adjacency[source]]
                virtual_source_candidate, current_neighbors, source_likelihood = pick_random_elements(current_neighbors,1)
                previous_vs = virtual_source
                
                # infect twice in one direction, always
                self.infect_nodes_up_down(source, virtual_source_candidate)
                self.directions[source] = False
                virtual_source_candidate = virtual_source_candidate[0]
                self.pass_branch_message(source, virtual_source_candidate)
                virtual_source = virtual_source_candidate
            else:
                current_neighbors = [k for k in self.who_infected[virtual_source]]
                if (len(current_neighbors) < 2):
                    # if there is nowhere for the virtual source to move, keep it where it is
                    if len(current_neighbors) < 1:
                        blocked = True
                        print('Blocked! Exiting.')
                        break
                        
                    # branch once in every direction
                    for neighbor in current_neighbors:
                        self.pass_branch_message(virtual_source, neighbor)
                    # take care of getting stuck
                    if len(current_neighbors) == 1:
                        self.pass_branch_message(virtual_source, current_neighbors[0])
                        previous_vs = virtual_source
                        virtual_source = current_neighbors[0]
                
                else:           # spread asymmetrically
                    # find a direction to move
                    virtual_source_candidate = [previous_vs]
                    cands = [neighbor for neighbor in current_neighbors if (self.directions[neighbor] == True) and (neighbor != previous_vs)]
                    if not cands:
                        print('We\'re blocked! Exiting.')
                        break
                    virtual_source_candidate = cands[0]
                    previous_vs = virtual_source
                    # the virtual source moves one more hop away from the true source
                
                    # branch twice in one direction
                    self.pass_branch_message(virtual_source, virtual_source_candidate)
                    self.pass_branch_message(virtual_source, virtual_source_candidate)
                    
                    virtual_source = virtual_source_candidate
            tot_num_infected[timesteps] = self.num_infected
            timesteps += 1

        print('Reached ', self.num_infected,' nodes out of ',self.num_nodes,'. Adjacency has ',len(self.adjacency))
        print('Reached ', len(self.spies_info.spies),' spies, ', len(self.spies_info.active_spies), ' of which are active.')
        if len(self.spies_info.active_spies)>0:
            print('Source:',source, 'First spy info:',self.spies_info.active_spies[0].level,self.spies_info.active_spies[0].up_node)
        infection_details = (self.who_infected, tot_num_infected)
            
        return infection_details
    
    def pass_branch_message(self, source, recipient):
        # pass an instruction to branch from the source to the leaves
        # Inputs
        #       source:             source of the infection
        #       recipient:          array of child ids
        #       infection_pattern:  binary array describing whether each node is already infected or not
        #       adjacency:          adjacency relations for the underlying network
        #       max_infection:      maximum number of people who can be infected by a single node     
        #       who_infected:       adjacency relations for the infected subgraph
        #
        # Outputs
        #       infection_pattern   (updated)
        #       who_infected        (updated)
        
        leaf = True
        
        # pass to the neighbors who are your neighbors in the true infection graph
        neighbors = [k for k in self.who_infected[recipient] if (not k == source)]
        neighbors.sort() # this should be commented for a randomized distribution!
        
        for neighbor in neighbors:
            leaf = False
            self.pass_branch_message(recipient, neighbor)
                
        if leaf:
            neighbors = [k for k in self.adjacency[recipient] if self.infection_pattern[k]==0]
            if len(neighbors) > self.max_infection:
                neighbors, remains, leaf_likelihood = pick_random_elements(neighbors,self.max_infection)
            if neighbors:
                # There are neighbors to infect 
                self.infect_nodes_up_down(recipient, neighbors)
     
    def infect_nodes_up_down(self, node, children):
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
        
        if self.levels[node] == 0 and self.directions[node] == False:
            return
            
        #update the adjacency matrix
        self.who_infected[node] += children
        for child in children:
            self.who_infected[child].append(node)
            self.infection_pattern[child] = 1
        self.num_infected = sum(self.infection_pattern)        
        
        #update the metadata
        if self.levels[node] == 0 and self.directions[node] == True:   # start node
            for child in children:
                self.levels[child] = 1
                self.directions[child] = True
        elif self.directions[node] == True:     # up chain
            up_node = random.choice(children)
            for child in children:
                if child == up_node:
                    self.levels[child] = self.levels[node] + 1
                    self.directions[child] = True
                else:
                    self.levels[child] = self.levels[node]-1
                    self.directions[child] = False
        elif self.directions[node] == False:     # down chain
            for child in children:
                self.levels[child] = self.levels[node] - 1
                self.directions[child] = False
                
        #Update the spies
        self.spies_info.add_nodes(node, children, self.levels, self.directions, self.spy_probability)
        
        
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
        self.degrees = self.degrees_rv.draw_values(1)
        self.levels = [None]
        # infect the first neighbor, always
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
        
        
        infection_details = (self.who_infected, tot_num_infected)
            
        
        return infection_details
        
        
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
        # self.degrees += self.degrees_rv.rvs(size=len(children)).tolist()
        self.degrees += self.degrees_rv.draw_values(len(children))
        
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