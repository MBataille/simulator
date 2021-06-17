from equations.equation import Equation
import numpy as np

class Pendulum(Equation):

    def __init__(self):
        initParams = { 'dt' : 0.05,
                    'mu': 0.1,
                    'gamma' : 2.75,
                    'omega0' : 1,
                    'omega' : 0.7}
        
        initRange = { 'dt' : (0, 0.2),
                    'mu': (0, 0.5),
                    'gamma' : (0, 5),
                    'omega0' : (0, 2),
                    'omega' : (0, 2)}

        Equation.__init__(self, 'Pendulum', initParams, dim=0, n_fields=2, N=1, fieldNames=['theta', 'dtheta/dt'], initRange=initRange)

    def rhs(self, t, X, dXdt):
        v = self.getCurrentParams()

        mu = v['mu']; gamma = v['gamma']; omega = v['omega']
        omega0 = v['omega0']

        d2Xdt2 = - omega0 ** 2 * np.sin(X) - mu * dXdt \
                        + gamma * np.cos(omega * t)

        return (dXdt, d2Xdt2)