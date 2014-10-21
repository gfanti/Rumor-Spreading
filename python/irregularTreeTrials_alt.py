import scipy.io as io
import numpy as np
from scipy import stats
import runExperiments
import utilities

'''Run the irregular tree algorithm'''
if __name__ == "__main__":

    trials = 500
    max_time = 6
    
    # Irregular infinite graph
    xk = np.arange(2,5)
    pk = (0.25,0.5,0.25)
    dd = 0.5
    ds = np.arange(1.0,4.05,dd)
        
    num_infected_all = []
    pd_ml_all = []
    
    for max_infection in [2]:

        print('Checking d_o = ',max_infection+1)
        degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
        num_infected, pd_ml = runExperiments.run_randtree(trials, max_time, max_infection, degrees_rv, 1)[:2]
        
        num_infected_all.append(num_infected)
        pd_ml_all.append(pd_ml)
        
    xk_str = [str(i) for i in xk]
    filename = 'results/irregular_tree_alt_results_' + "_".join(xk_str) + '.mat'
    io.savemat(filename,{'pd_ml':np.array(pd_ml_all), 'num_infected':np.array(num_infected_all), 'd_values':ds})
