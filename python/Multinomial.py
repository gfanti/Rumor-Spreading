from scipy import stats

class Multinomial:
    def __init__(self, xk, pk):
        self.xk = xk
        self.pk = pk
        self.variable = stats.rv_discrete(name='rv_discrete', values=(xk, pk))
    
    def mean(self):
        return self.variable.mean()
    
    def draw_values(self,howmany):
        return self.variable.rvs(size=howmany).tolist()