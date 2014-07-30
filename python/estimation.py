import random 

def rumor_centrality_up(calling_node, called_node):
    if called_node == calling_node:
        for i in who_infected[called_node]:
                rumor_centrality_up(called_node, i)
        up_messages[called_node][0] += 1
        up_messages[called_node][1] = up_messages[called_node][0] * up_messages[called_node][1]
        up_messages[calling_node][0] += up_messages[called_node][0]
        up_messages[calling_node][1] = up_messages[calling_node][1] * up_messages[called_node][1]
    elif len(who_infected[called_node]) == 1:   # leaf node
        up_messages[calling_node][0] += 1 # check
    else:
        for i in who_infected[called_node]:
            if i != calling_node:
                rumor_centrality_up(called_node, i)
        up_messages[called_node][0] += 1
        up_messages[called_node][1] = up_messages[called_node][0] * up_messages[called_node][1]
        up_messages[calling_node][0] += up_messages[called_node][0]
        up_messages[calling_node][1] = up_messages[calling_node][1] * up_messages[called_node][1]
        

def rumor_centrality_down(calling_node, called_node):
    if called_node == calling_node:
        for i in who_infected[called_node]:
            down_messages[called_node] = down_messages[called_node] * up_messages[i][1]
        for i in who_infected[called_node]:
            down_messages[i] = down_messages[called_node]*(up_messages[i][0]/(K-up_messages[i][0]))
            rumor_centrality_down(called_node, i)
    elif len(who_infected[called_node]) == 1:   # leaf node
        down_messages[calling_node] = down_messages[called_node]*(up_messages[calling_node][0]/(K-up_messages[calling_node][0]))
    else:
        for i in who_infected[called_node]:
            if i != calling_node:
                down_messages[i] = down_messages[called_node]*(up_messages[i][0]/(K-up_messages[i][0]))
                rumor_centrality_down(called_node, i)


def rumor_centrality():

    initial_node = 6       # can use arbitrary initial index
    rumor_centrality_up(initial_node,initial_node)
    rumor_centrality_down(initial_node,initial_node)
    max_down = max(down_messages)
    max_down_ind = [i for i, j in enumerate(down_messages) if j == max_down]
    return max_down_ind[random.randrange(0,len(max_down_ind),1)]

global who_infected
global up_messages
global down_messages
global K

who_infected = [ [] for i in range(7)]
who_infected[0] = [1,2]
who_infected[1] = [0,3,4]
who_infected[2] = [0]
who_infected[3] = [1,5]
who_infected[4] = [1,6]
who_infected[5] = [3]
who_infected[6] = [4]

K = len(who_infected)
up_messages = [ [0,1] for i in range(K)]
down_messages = [1]*K

print rumor_centrality()

