from equations.equation import Equation, Parameter, SOLVE_EVERY_TI
import numpy as np
from scipy.integrate import simps

class Duffing1D(Equation):

    def __init__(self):
        initParams = { 'dt' : 0.05,
                    'dx': 1, 
                    'alpha': 0.4,
                    'mu': 0.1,
                    'gamma' : 2.75,
                    'omega' : 0.7}
        
        initRange = { 'dt' : (0, 0.2),
                    'dx': (0, 1), 
                    'alpha': (0,1),
                    'mu': (0, 0.5),
                    'gamma' : (0, 5),
                    'omega' : (0, 1)}

        Equation.__init__(self, 'Duffing1D', initParams, dim=0, n_fields=2, N=1, fieldNames=['x', 'dx/dt'], initRange=initRange)

    def rhs(self, t, X, dXdt):
        v = self.getCurrentParams()

        dx = v['dx']; alpha = v['alpha']; mu = v['mu']; 
        gamma = v['gamma']; omega = v['omega']

        d2Xdt2 = -X + alpha * X ** 3 - X ** 5 - mu * dXdt \
                        + gamma * np.cos(omega * t)

        return (dXdt, d2Xdt2)