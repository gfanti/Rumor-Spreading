import scipy.io as io
import numpy as np
from scipy import stats
import utilities
import runExperiments

'''Run the irregular tree algorithm'''
if __name__ == "__main__":

    trials = 100
    max_time = 2
    max_infection = 99
    
    # Irregular infinite graph
    xk = np.arange(3,5)
    pk = (0.5,0.5)
    dd = 0.05
    # ds = np.arange(1.0,4.05,dd)
    ds = np.array([2])
    
    # Regular infinite graph
    # xk = np.arange(3,4)
    # pk = (1)
    # ds = np.array([2])
    
    # want to find the best max_infection (i.e., d-1) to minimize pd
    
    num_infected_all = []
    pd_ml_all = []
    pd_rand_leaf_all = []
    
    degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
    
    total_pd_ml = 0.0
    total_pd_rand_leaf = 0.0
    pd_ml_progression = []
    pd_rand_leaf_progression = []
    
    for t in range(trials):
        # build up the degrees vector: degrees with which to infect new nodes
        # 1st infected node
        degrees = degrees_rv.rvs(size=1).tolist()  # [x]
        # T = 0.5
        degrees += degrees_rv.rvs(size=1).tolist()  # [x]
        # T = 1
        degrees += [3] + degrees_rv.rvs(size=2).tolist()  # [3,x,x]
        # T = 1.5
        if degrees[1] == 4:
            degrees += [3,3]  # depends on T = 0.5, should sum to 10
        else:
            degrees += degrees_rv.rvs(size=1).tolist()  # [x]
            degrees += [10 - degrees[1] - degrees[-1]]  # depends on T = 0.5, should sum to 10
        # T = 2
        remaining = sum(degrees[-2:]) - 2
        degrees += degrees_rv.rvs(size=remaining).tolist()  # [x,x,x,x]
        
        num_infected, pd_ml, pd_rand_leaf = runExperiments.run_randtree(1, max_time, max_infection, degrees_rv, degrees)
        total_pd_ml += pd_ml[-1]
        total_pd_rand_leaf += pd_rand_leaf[-1]
        # print('pd', pd_ml)
        
        num_infected_all.append(num_infected)
        pd_ml_all.append(pd_ml)
        pd_rand_leaf_all.append(pd_rand_leaf)
        
        pd_ml_progression.append(total_pd_ml / (float(t)+1))
        pd_rand_leaf_progression.append(total_pd_rand_leaf / (float(t)+1))
    
    total_pd_ml /= float(trials)
    total_pd_rand_leaf /= float(trials)
    print('Overall pd_ml is ',total_pd_ml)
    print('1/m is ',1/7)
    
    io.savemat('results/irregular_tree_fixed_m_10',{'pd_ml':np.array(total_pd_ml), 'pd_progression_ml':np.array(pd_ml_progression), 'pd_progression_rand_leaf':np.array(pd_rand_leaf_progression), 'trials':np.array(trials)})
