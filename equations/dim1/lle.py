from equations.equation import Equation, Parameter, SOLVE_EVERY_TI
from threadpoolctl import threadpool_limits
from scipy.integrate import solve_ivp

import scipy.sparse as sp
import numpy as np
import findiff as fd

T0_UPPER_BOUND = 1e6

def get_finite_difference_coefficients(order, acc):
    d = fd.coefficients(order, acc)['center']
    return d['coefficients'], d['offsets']

def derivative(u:np.ndarray, dt:float, order=1, acc=2):
    coefs, offsets = get_finite_difference_coefficients(order, acc)
    du = np.zeros_like(u)
    for i, coef in enumerate(coefs):
        offset = i - (len(coefs) - 1) // 2
        du += coef * np.roll(u, -offset)
    return du / dt ** order

def better_laplace(X, dx):
    return derivative(X, dx, order=2, acc=8)

def diagonal_matrix(coef, n_t, offset):
    return np.diag(coef * np.ones(n_t-abs(offset)), k=offset)

def derivative_matrix(n_t, dt, order=1, acc=2, sparse=False, kind='center'):
    coefs_dict = fd.coefficients(order, acc)[kind]
    coefs = coefs_dict['coefficients']
    offsets = coefs_dict['offsets']

    D = np.zeros((n_t, n_t))
    for coef, offset in zip(coefs, offsets):
        D += diagonal_matrix(coef, n_t, offset)

        # periodic boundary conditions
        aoff = abs(offset)
        for i in range(aoff):
            j = n_t - aoff + i
            row, col = (i, j) if offset < 0 else (j, i)
            D[row, col] = coef

    D = D / dt ** order
    if sparse:
        return sp.csc_matrix(D)
    return D

class LugiatoLefeverEquation(Equation):

    def __init__(self):
        initParams = { 'dt' : 0.01,
                    'dx': 0.25, 
                    'Delta': 6.25,
                    'S': 2.8,
                    'beta2': -1}

        Equation.__init__(self, 'LugiatoLefever', initParams, dim=1, n_fields=2, N=200, fieldNames=['Re A', 'Im A'], auxFieldNames=['mod A'])

    def rhs(self, t, Re, Im):
        v = self.getCurrentParams()

        A = Re + 1j * Im
        Delta = v['Delta']; S = v['S']; beta2 = v['beta2']

        dA = S - (1 + 1j * Delta) * A - 1j * beta2 * better_laplace(A, v['dx']) + 1j * (Re ** 2 + Im ** 2) * A
        return dA.real, dA.imag
    
    def setNi(self, Ni):
        super().setNi(Ni)

        dx = self.getCurrentParams()['dx']
        D = derivative_matrix(self.Ni, dx, order=2, acc=8, sparse=True)
        self.Dxx = sp.kron(np.array([[0, 1], [-1, 0]]), D, format='csc')
    
    def jac(self, t, Y):
        u, v = self.extractFields(Y)

        p = self.getCurrentParams()
        delta, beta2, dx = p['Delta'], p['beta2'], p['dx']

        principal_diag = np.append(-2 * u * v - 1, 2 * u * v - 1)
        upper_diag = delta - u * u - 3 * v * v
        lower_diag = -delta + 3 * u * u + v * v

        jac_homo = sp.diags([principal_diag, lower_diag, upper_diag],
                    offsets=[0, -self.Ni, self.Ni], format='csc')
        return jac_homo + beta2 * self.Dxx

    def getAuxFields(self, Re, Im):
        return [Re * Re + Im * Im]
    
    def solve(self, t=None):
        if self.t0 > T0_UPPER_BOUND:
            self.t0 = 0

        dt = self.getParam('dt')
        if t is None:
            t = self.t0 + np.linspace(0, SOLVE_EVERY_TI*dt, SOLVE_EVERY_TI + 1) # change this
        else:
            t += self.t0
        # print(t, dt, self.initCond)
        with threadpool_limits(limits=1):
            self.sol = solve_ivp(self.wrhs, (t[0], t[-1]), self.initCond, method='Radau', t_eval=t, jac=self.jac).y
