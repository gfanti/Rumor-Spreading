from infectionUtils import *
from multinomial import *
import networkx as nx

class RandTreeInfector(Infector):

	def __init__(self, source, max_infection, max_time, degrees_rv):
		super(RandTreeInfector, self).__init__()

		self.source = source
		self.max_infection = max_infection
		self.max_time = max_time
		self.degrees_rv = degrees_rv

		self.timesteps = 0

		self.who_infected = [[]]
		self.tot_num_infected = []
		self.degrees = degrees_rv.draw_values(1)
		self.num_infected = 1

	def pass_branch_messages(self, source, recipient):
		# pass an instruction to branch from the source to the leaves
		# Inputs
		#       source:             source of the infection
		#       recipient:          array of child ids

		leaf = True

		# pass to the neighbors who are your neighbors in the true infection graph
		neighbors = [k for k in self.who_infected[recipient] if (not k == source)]
		neighbors.sort() # this should be commented for a randomized distribution!

		for neighbor in neighbors:
			leaf = False
			self.pass_branch_messages(recipient, neighbor)

		if leaf:
			to_infect = self.degrees[recipient] - 1
			neighbors = ([k+len(self.degrees) for k in range(to_infect)])
			self.infect_nodes(recipient, neighbors)

#Adaptive diffusion infector over random trees
class RandTreeADInfector(RandTreeInfector):
	def __init__(self, source, max_infection, max_time, degrees_rv):
		super(RandTreeADInfector, self).__init__(source, max_infection, max_time, degrees_rv)
		self.d_o = max_infection + 1 # assumed regular degree for choosing alphas


	def infect_one_timestep(self):

		if self.timesteps < self.max_time:

			if self.timesteps == 0:
				self.virtual_source = 1
				self.previous_vs = 0

				# infect twice in one direction, always
				# 1) Infect the new VS
				self.infect_nodes(self.source, [self.virtual_source])
				# 2) Pass the branching messages
				self.pass_branch_messages(self.source, self.virtual_source)
				self.m = 1
			else:
				current_neighbors = [k for k in self.who_infected[self.virtual_source]]
				if random.random() < utilities.compute_alpha(self.m, self.timesteps, self.d_o):
					# branch once in every direction
					for neighbor in current_neighbors:
						self.pass_branch_messages(self.virtual_source, neighbor)
				else:
					# pick a direction and branch multiple times
					current_neighbors = [k for k in self.who_infected[self.virtual_source] if (k != self.previous_vs)]	
					virtual_source_candidate = random.choice(current_neighbors)

					self.previous_vs = self.virtual_source

					# branch twice in one direction
					self.pass_branch_messages(self.virtual_source, virtual_source_candidate)
					self.pass_branch_messages(self.virtual_source, virtual_source_candidate)

					self.virtual_source = virtual_source_candidate
					self.m += 1

		self.num_infected = len(self.who_infected)
		self.tot_num_infected.append(self.num_infected)
		self.timesteps += 1

	def infect_nodes(self, node, children):
		'''Infects the children nodes.'''
		#update the adjacency matrix
		self.who_infected[node] += children
		for child in children:
		    self.who_infected.append([node])

		self.degrees += self.degrees_rv.draw_values(len(children))
		self.num_infected = len(self.who_infected)
	    
        
        

#Preferential attachment adaptive diffusion (PAAD) infector over random trees
class RandTreePAADInfector(RandTreeInfector):
	def __init__(self, source, max_infection, max_time, degrees_rv, num_hops_pa = 1):
		super(RandTreePAADInfector, self).__init__(source, max_infection, max_time, degrees_rv)
		self.num_hops_pa = num_hops_pa
	    

	def infect_one_timestep(self):

		if self.timesteps < self.max_time:

			if self.timesteps == 0:
				self.virtual_source = 1
				self.previous_vs = 0

				# infect twice in one direction, always
				# 1) Infect the new VS
				self.infect_nodes(self.source, [self.virtual_source])
				# 2) Pass the branching messages
				self.pass_branch_messages(self.source, self.virtual_source)
			else:
				current_neighbors = [k for k in self.who_infected[self.virtual_source] if (k != self.previous_vs)]	
				# current_neighbors.remove(self.previous_vs)
				weights = [self.degrees[i]-1 for i in current_neighbors]
				weights = [float(i)/sum(weights) for i in weights]
				virtual_source_rv = Multinomial(current_neighbors, weights)
				virtual_source_candidate = virtual_source_rv.draw_values(1)[0]

				self.previous_vs = self.virtual_source


				# branch twice in one direction
				self.pass_branch_messages(self.virtual_source, virtual_source_candidate)
				self.pass_branch_messages(self.virtual_source, virtual_source_candidate)
				# degrees, who_infected = pass_branch_message_randtree(virtual_source, virtual_source_candidate, degrees, degrees_rv, who_infected)[:2]
				# degrees, who_infected = pass_branch_message_randtree(virtual_source, virtual_source_candidate, degrees, degrees_rv, who_infected)[:2]

				self.virtual_source = virtual_source_candidate

		self.num_infected = len(self.who_infected)
		self.tot_num_infected.append(self.num_infected)
		# tot_num_infected[timesteps] = num_infected
		self.timesteps += 1


	def infect_nodes(self, node, children):
		'''Infects the children nodes.'''
		#update the adjacency matrix
		self.who_infected[node] += children
		for child in children:
		    self.who_infected.append([node])

		if (node == 0): # we treat the zero node differently and keep track of its neighbors
			self.src_neighbors = self.degrees_rv.draw_values(self.degrees[node])
			weights = [float(j-1) / (sum(self.src_neighbors) - self.degrees[node]) for j in self.src_neighbors]
			vs = np.random.choice(self.src_neighbors, p = weights)
			self.src_neighbors.remove(vs)
			self.degrees += [vs]
		else:
			self.degrees += self.degrees_rv.draw_values(len(children))

		self.num_infected = len(self.who_infected)

        
    