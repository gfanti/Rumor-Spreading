import random


# build a barabasi-albert graph
def buildBarabasiAlbertGraph(num_nodes):
    adjacency = [[] for i in range(num_nodes)]
    total_degree = 0;
    
    for i in range(num_nodes):
        for j in range(i+1,num_nodes):
            probability = random.random() # draw a random number
            
            attractiveness = (len(adjacency[i]))/(total_degree + 0.001) + 0.1;
            
            if (probability < attractiveness) or (total_degree == 0):
                adjacency[i].append(j)
                adjacency[j].append(i)
                
                total_degree += 1
                
    return adjacency
    
# build a small-world graph
def buildSmallWorldGraph(num_nodes, d):
    
    # Inputs:
    #       num_nodes -- number of nodes in the graph
    #       d         -- number of neighbors, divided by 2
    #
    # Outputs:
    #       adjacency -- list of adjacency relations for all nodes
    
    adjacency = [[] for i in range(num_nodes)]
    
    # probability of rewiring
    beta = 0.2
    
    # adjacency set
    S = [i for i in range(-d,0)] + [i for i in range(1,d+1)]
    
    # build up the circulant graph
    for i in range(num_nodes):
        adjacency[i] = [(i+j)%num_nodes for j in S]
    # print('circulant: ',adjacency)
        
    for i in range(num_nodes):
        for j in range(len(adjacency[i])):
            if random.random() < beta:
                oldnode = adjacency[i][j]
                newnode = oldnode
                # if a node is already friends with everyone, you can't rewire
                if len(adjacency[i]) == num_nodes-1:
                    continue
                # otherwise, choose a new edge
                while (newnode == oldnode) or any([k for k in adjacency[i] if k==newnode]) or (newnode == i):
                    newnode = random.randint(0,num_nodes-1)
                # print('newnode: ',newnode)
                # print('i: ',i)
                # print('adjacency[j]')
                adjacency[i].pop(j)
                adjacency[i].append(newnode)
                adjacency[oldnode].remove(i)
                adjacency[newnode].append(i)
                
    return adjacency
    
    
def buildDatasetGraph(filename, min_degree):
    num_nodes = 4941    # hard-coded number of nodes in this dataset
    adjacency = [[] for i in range(num_nodes)]
    
    # open the datafile
    f = open(filename,'rb')
    # get rid of the 2 opening lines
    f.readline()
    f.readline()
    
    edges = f.readlines()
    
    # add all the edges
    for edge in edges:
        edge = edge.split()
        source = int(edge[0]) - 1
        destination = int(edge[1]) - 1
        if (destination < num_nodes):
            adjacency[source].append(destination)
            adjacency[destination].append(source)
    
    # zero out the people with fewer than min_degree friends
    while True:
        loopflag = True
        for i in range(len(adjacency)):
            if len(adjacency[i]) < min_degree and len(adjacency[i]) > 0:
                loopflag = False
                for node in adjacency[i]:
                    adjacency[node].remove(i)
                adjacency[i] = []
        if loopflag:
            break
    
    return adjacency