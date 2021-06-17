from equations.equation import Equation
import numpy as np

class Lorenz(Equation):

    def __init__(self):
        initParams = { 'dt' : 0.05,
                    'sigma': 10,
                    'rho' : 28,
                    'beta' : 8/3}
        
        initRange = { 'dt' : (0, 0.2),
                    'sigma': (0, 20),
                    'rho' : (0, 40),
                    'beta' : (0, 5)}

        Equation.__init__(self, 'Lorenz', initParams, dim=0, n_fields=3, N=1, fieldNames=['x', 'y', 'z'], initRange=initRange)

    def rhs(self, t, X, Y, Z):
        v = self.getCurrentParams()

        sigma = v['sigma']; rho = v['rho']; beta = v['beta']

        dXdt = sigma * (Y - X)
        dYdt = X * (rho - Z) - Y
        dZdt = X * Y - beta * Z

        return (dXdt, dYdt, dZdt)