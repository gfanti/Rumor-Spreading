import scipy.io as io
import numpy as np
from scipy import stats
import utilities

'''Run the irregular tree algorithm'''
if __name__ == "__main__":

    trials = 100
    max_time = 10
    max_infection = 3
    xk = np.arange(3,5)
    pk = (0.5,0.5)
    degrees_rv = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
    num_infected, pd_ml = utilities.run_randtree(trials, max_time, max_infection, degrees_rv)
    print('Max likelihood Pd: ',pd_ml)