import scipy.io as io
import numpy as np
from scipy import stats
import runExperiments
import utilities

'''Run the irregular tree algorithm'''
if __name__ == "__main__":

    trials = 10000
    max_time = 9
    
    # Irregular infinite graph
    xk = np.arange(2,4)
    pk = (0.5,0.5)
    dd = 0.5
    ds = np.arange(1.0,4.05,dd)
    
    # Regular infinite graph
    # xk = np.arange(3,4)
    # pk = (1)
    # ds = np.array([2])
    
    # want to find the best max_infection (i.e., d-1) to minimize pd
    
    # max_infection = 100000000
    # max_infection = min(xk)
    # max_infection = max(xk)
    # max_infection = sum([i*j for i,j in zip(xk,pk)])
    
    num_infected_all = []
    pd_ml_all = []
    
    for max_infection in ds.tolist():

        print('Checking d_o = ',max_infection+1)
        degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
        num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 1)[:2]
        
        num_infected_all.append(num_infected)
        pd_ml_all.append(pd_ml)
        
    io.savemat('results/irregular_tree_alt_results',{'pd_ml':np.array(pd_ml_all), 'num_infected':np.array(num_infected_all), 'd_values':ds})
