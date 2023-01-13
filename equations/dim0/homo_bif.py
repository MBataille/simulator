from equations.equation import Equation
import numpy as np

class HomoclinicBifurcation(Equation):

    def __init__(self):
        initParams = { 'dt' : 0.01,
                    'epsilon': 0.5,
                    'mu' : -0.1}
        
        initRange = { 'dt' : (0, 0.2),
                    'epsilon': (-1.0, 1.0),
                    'mu' : (-1.0, 1)}

        Equation.__init__(self, 'Lorenz', initParams, dim=0, n_fields=2, N=1, fieldNames=['x', 'dxdt'], initRange=initRange)

    def rhs(self, t, X, Y):
        v = self.getCurrentParams()

        epsilon = v['epsilon']; mu = v['mu']

        dXdt = Y
        dYdt = epsilon * X - X ** 2 + mu * Y + X * Y

        return (dXdt, dYdt)