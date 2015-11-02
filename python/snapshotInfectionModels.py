from infectionUtils import *
from multinomial import *
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.traversal.depth_first_search import dfs_tree


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
	    
        
     
class LocalGraph(nx.Graph):
	def __init__(self, infector):
		super(LocalGraph, self).__init__()
		self.infector = infector
		self.source = self.infector.source
		self.add_node(self.source)

	def generate_local_graph(self):

		boundary_old = [j for j in range(-1, -(self.infector.degrees[0] + 1), -1)]
		for neighbor in boundary_old:
			self.add_edge(self.source, neighbor)
		depth = 1
		self.node_count = -(self.infector.degrees[0] + 1)
		self.expand_local_graph(depth, boundary_old)
		self.choose_initial_labeling()


	def expand_local_graph(self, depth, boundary_old):
		boundary_new = []
		while True:
			if (depth > self.infector.num_hops_pa):
				break
			while boundary_old:
				n = boundary_old.pop(0)
				if self.degree(n) == 1:
					deg_n = self.infector.degrees_rv.draw_values(1)[0] - 1
					for i in range(deg_n):
						self.node_count -= 1
						self.add_edge(n, self.node_count)
						boundary_new += [self.node_count]
				else:
					boundary_new += self.neighbors(n)
			depth += 1
			boundary_new, boundary_old = [], boundary_new

	def get_hop_degrees(self, node, neighbor, hops_left):
		# print("looking from node", node, " to neighbor", neighbor, "hops_left",hops_left)
		current_neighbors = [k for k in self.neighbors(neighbor) if (k != node)]

		if hops_left == 0:
			return len(current_neighbors)
		else:
			sum_degrees = 0
			for n in current_neighbors:
				sum_degrees += self.get_hop_degrees(neighbor, n, hops_left-1)
			return sum_degrees

	def choose_initial_labeling(self):
		''' Choose a labeling of the local graph that respects the n-hop decision rule'''
		source_neighbors = self.neighbors(self.source)

		# Choose the next virtual source
		weights = []
		for neighbor in source_neighbors:
			subtree = dfs_tree(self, neighbor)
			leaves = [i for i in subtree.nodes() if (subtree.out_degree(i)==0 and \
													 not nx.has_path(subtree,self.source,i))]
			weights.append(len(leaves))
		weights = [float(i)/sum(weights) for i in weights]
		virtual_source = np.random.choice(source_neighbors, p=weights)
		
		# Relabel the nodes so that only the infected nodes have positive values
		label = 1
		mapping = {virtual_source: label}
		# Generate a directed tree emanating from the virtual source away from the source
		directed_self = dfs_tree(self, self.source)
		infected_subgraph = dfs_tree(directed_self, virtual_source)
		bfs_edges = list(nx.bfs_edges(infected_subgraph, virtual_source))
		# In this first timestep, we only need the direct neighbors of V.S. to be infected
		bfs_edges = [item[1] for item in bfs_edges if (item[1] != self.source and item[0] == virtual_source)]
		for edge in bfs_edges:
			label += 1
			mapping[edge] = label
		self = nx.relabel_nodes(self, mapping, copy=False)

		

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
				self.local_neighborhood = LocalGraph(self)
				self.local_neighborhood.generate_local_graph()

				self.infect_nodes(self.source, [self.virtual_source])
				# 2) Pass the branching messages
				self.pass_branch_messages(self.source, self.virtual_source)
				self.m = 1
			else:
				# Compute the next v.s. based on the existing graph and the source's local_neighborhood
				current_neighbors = [k for k in self.who_infected[self.virtual_source] if (k != self.previous_vs)]	
				weights = [self.local_neighborhood.get_hop_degrees(self.virtual_source, i, self.num_hops_pa) for i in current_neighbors]
				weights = [float(i)/sum(weights) for i in weights]
				virtual_source_rv = Multinomial(current_neighbors, weights)
				virtual_source_candidate = virtual_source_rv.draw_values(1)[0]

				self.previous_vs = self.virtual_source


				# branch twice in one direction
				self.pass_branch_messages(self.virtual_source, virtual_source_candidate)
				self.pass_branch_messages(self.virtual_source, virtual_source_candidate)
				
				self.virtual_source = virtual_source_candidate
				self.m += 1
		self.num_infected = len(self.who_infected)
		self.tot_num_infected.append(self.num_infected)
		self.timesteps += 1
		
		# nx.draw_networkx(self.local_neighborhood)
		# plt.show()	


	def infect_nodes(self, node, children, num_hops_pa = 1):
		'''Infects the children nodes.'''
		

		boundary_old = [child for child in children]
		#update the adjacency matrix
		self.who_infected[node] += children
		for child in children:
			self.who_infected.append([node])

		# Rename the child nodes to positive nodes
		mapping = {}
		if any([child not in self.local_neighborhood.neighbors(node) for child in children]):
			for neighbor in self.local_neighborhood.neighbors(node):
				if neighbor < 0:
					mapping[neighbor] = children.pop(0)
			self.local_neighborhood = nx.relabel_nodes(self.local_neighborhood, mapping, copy=False)
		
		# Generate the next num_hops_pa levels
		children = [child for child in boundary_old]
		self.local_neighborhood.expand_local_graph(0, boundary_old)
		
		# Update the degrees vector
		for child in children:
			self.degrees += [self.local_neighborhood.degree(child)]

		self.num_infected = len(self.who_infected)


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

        
    