from equations.equation import Equation
import numpy as np

class VanDerPol(Equation):

    def __init__(self):
        initParams = { 'dt' : 0.01,
                    'mu': 0.1}
        
        initRange = { 'dt' : (0, 0.2),
                    'mu': (-1, 5)}

        Equation.__init__(self, 'VanDerPol', initParams, dim=0, n_fields=2, N=1, fieldNames=['x', 'dx/dt'], initRange=initRange)

    def rhs(self, t, X, dXdt):
        v = self.getCurrentParams()

        mu = v['mu']

        d2Xdt2 = mu * (1 - X ** 2) * dXdt - X

        return (dXdt, d2Xdt2)