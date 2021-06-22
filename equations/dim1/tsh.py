from equations.equation import Equation, Parameter
import numpy as np

class TuringSwiftHohenberg(Equation):

    def __init__(self):
        initParams = {
            'dt' : 0.1,
            'dx' : 0.5,
            'eps': -0.22,
            'nu': 1
        }
        Equation.__init__(self, 'Turing-Swift-Hohenberg', initParams, fieldNames=['u'])

    def rhs(self, t, u):
        v = self.getCurrentParams()
        eps = v['eps']; nu = v['nu']; dx = v['dx']

        r = np.roll
        u_nn = r(u, -1) + r(u, 1)

        # This whay it should be *slightly* more efficient
        laplace = (u_nn - 2 * u)/(dx**2)
        bilaplace = (r(u, 2) + r(u, -2) - 4 * u_nn + 6 * u) / (dx**4)

        if self.boundary_condition == 'specular':
            laplace[0] = 2*(u[1] - u[0])/(dx**2)
            laplace[-1] = 2*(u[-2] - u[-1])/(dx**2)
            bilaplace[0] = 2*(u[2] - 4*u[1] + 3*u[0])/(dx**4)
            bilaplace[1] = (u[1] + u[3] - 4*u_nn[1] + 6*u[1])/(dx**4)
            bilaplace[-1] = 2*(u[-3] - 4*u[-2] + 3*u[-1])/(dx**4)
            bilaplace[-2] = (u[-2] + u[-4] - 4*u_nn[-2] + 6*u[-2])/(dx**4)

        #return eps*u - u**3 - nu*self.Laplace1D(u) - self.BiLaplace1D(u)
        return eps*u - u ** 3 - nu * laplace - bilaplace

    def getAllBoundaryConditions(self):
        return ['periodic', 'specular']