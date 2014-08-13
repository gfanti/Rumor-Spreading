import scipy.io as io
import numpy as np
from scipy import stats
import utilities

'''Run the irregular tree algorithm'''
if __name__ == "__main__":

    trials = 6000
    max_time = 9
    max_infection = 3
    
    # Irregular infinite graph
    xk = np.arange(3,5)
    pk = (0.5,0.5)
    dd = 0.05
    ds = np.arange(1.0,4.05,dd)
    
    # Regular infinite graph
    # xk = np.arange(3,4)
    # pk = (1)
    # ds = np.array([2])
    
    # want to find the best max_infection (i.e., d-1) to minimize pd
    
    num_infected_all = []
    pd_ml_all = []
    
    for max_infection in ds.tolist():
        print('Checking d_o = ',max_infection+1)
        degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
        num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv)
        
        num_infected_all.append(num_infected)
        pd_ml_all.append(pd_ml)
    # print('Max likelihood Pd: ',pd_ml_all)
    #print('\nMean number of nodes: ',num_infected_all)
    
    io.savemat('results/irregular_tree_results',{'pd_ml':np.array(pd_ml_all), 'num_infected':np.array(num_infected_all), 'd_values':ds})
