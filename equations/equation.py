import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
from tkinter import StringVar
import os

from solvers import INTEGRATION_METHODS

DATAFOLDER = 'data/'
ALL_BOUNDARY_CONDITIONS = 'neumann periodic'.split(' ')

T0_UPPER_BOUND = 1e6

class Parameter:
    def __init__(self, name, var):
        self.name = name
        self.var = var
        self.val = float(var.get())

    def getVal(self):
        try:
            val = float(self.var.get())
            if val == 0 and cantBeZero(self.name):
                return self.val
            self.val = val
            return self.val
        except ValueError:
            return self.val

def cantBeZero(name):
    if name == 'dt' or name == 'dx':
        return True
    return False

class Equation:
    def __init__(self, name, initParams, dim=1, N=200, isComplex=False, n_fields=1, fieldNames=['u']):
        self.name = name
        self.dim = dim
        self.initParams = initParams
        self.parameters = self.createParamsDict(initParams)
        self.isComplex = isComplex
        self.fieldNames = fieldNames
        self.t0 = 0

        self.boundary_condition = 'neumann'

        # initialize w/ some solver
        self.solver = INTEGRATION_METHODS[0]

        self.n_fields = n_fields

        # each field will have N points
        self.setNi(N)

    def setNi(self, Ni):
        self.Ni = Ni
        self.N = self.n_fields * Ni

    def setSolver(self, solver):
        self.solver = solver

    def createParamsDict(self, initParams):
        params = {}
        for p in initParams:
            pvar = StringVar() # create var
            pvar.set(str(initParams[p])) # assign initial value
            params[p] = Parameter(p, pvar) # create parameter
        return params

    def getParam(self, p):
        return self.parameters[p].getVal()

    def getCurrentParams(self):
        cparams = {}
        for p in self.parameters:
            cparams[p] = self.parameters[p].getVal()
        return cparams

    def setInitialCondition(self, vals):
        if type(vals) == tuple:
            self.setInitCondFields(vals)
        else:
            self.initCond = vals

    def setInitCondIncoherent(self):
        self.setInitialCondition(0.99 * np.exp(1j * np.random.rand(self.N, self.N) * 2 * np.pi - np.pi))

    def setInitialConditionZero(self):
        dtype = 'complex128' if self.isComplex else 'float64'
        if self.dim == 1:
            self.initCond = np.zeros(self.N, dtype=dtype)
        if self.dim == 2:
            self.initCond = np.zeros((self.N, self.N), dtype=dtype)

    def solve(self):
        if self.t0 > T0_UPPER_BOUND:
            self.t0 = 0

        dt = self.getParam('dt')
        t = self.t0 + np.linspace(0, 60*dt, 61) # change this
        # print(t, dt, self.initCond)
        with threadpool_limits(limits=1):
            if self.dim == 1:
                self.sol = self.solver.solve(self.wrhs, (t[0], t[-1]), self.initCond, t)
            elif self.dim == 2:
                self.sol = self.solver.solve(self.wrhs, (t[0], t[-1]), self.initCond.reshape(self.N * self.N), t)
        self.t0 = t[-1]

    def wrhs(self, t, u):
        """Wrapper for RHS"""
        if self.n_fields == 1:
            return self.rhs(t, u)
        else:
            fields = self.extractFields(u)
            d_fields_dt = self.rhs(t, *fields)
            return self.assembleFields(d_fields_dt)
        

    def updateX(self):
        dx = self.parameters['dx'].getVal()
        Ni = round(self.N / self.n_fields)
        X = (Ni-1)*dx/2
        self.x = np.linspace(-X, X, Ni)

    def getX(self):
        self.updateX()
        return self.x

    def getState(self, k):
        if self.n_fields > 1:
            return self.getFields(k)
        if self.dim == 1:
            return self.sol[:, k]
        elif self.dim == 2:
            return np.reshape(self.sol[:, k], (self.N, self.N)) 

    def isState(self, filename):
        folder = self.getDataFolder()
        return os.path.exists(folder + filename + '.npz')

    def saveState(self, k, filename):
        folder = self.getDataFolder()
        _vals = self.sol[:, k]
        _pnames, _pvals = self.paramsToArray()
        np.savez_compressed(folder + filename + '.npz', vals=_vals, pnames=_pnames, pvals=_pvals)

    def loadState(self, filename):
        folder = self.getDataFolder()
        data = np.load(folder + filename + '.npz')

        self.setInitialCondition(data['vals'])
        self.Ni = int(data['pvals'][0])
        self.N = int(self.Ni * self.n_fields)
        self.initParams = self.arraytoInitParams(data['pnames'][1:], data['pvals'][1:])
        self.parameters = self.createParamsDict(self.initParams)

    def arraytoInitParams(self, pnames, pvals):
        initParams = {}
        for i in range(len(pnames)):
            initParams[pnames[i]] = pvals[i]
        return initParams

    def paramsToArray(self):
        v = self.getCurrentParams()
        pnames = ['Ni']
        pvals = [self.Ni]
        for p in v:
            pnames.append(p)
            pvals.append(v[p])
        return np.array(pnames), np.array(pvals)

    def getFields(self, k):
        return self.extractFields(self.sol[:, k])

    def getField(self, i, k):
        if self.n_fields == 1:
            return self.getState(k)
        Ni = self.Ni
        return self.sol[i * Ni:(i + 1) * Ni, k]

    def extractFields(self, X):
        Ni = self.Ni
        return [X[i * Ni : (i + 1) * Ni] for i in range(self.n_fields)]

    def assembleFields(self, fields):
        assembled = np.zeros(self.N)
        for i in range(len(fields)):
            assembled[i * self.Ni: (i+1) * self.Ni] = fields[i]
        return assembled

    def getInitCondFields(self):
        Ni = self.Ni
        print(Ni, self.N, self.n_fields)
        return [self.initCond[i * Ni : (i + 1) * Ni] for i in range(self.n_fields)]

    def setInitCondFields(self, fields):
        Ni = self.Ni
        self.initCond = self.assembleFields(fields)


    def getInitialConditions(self):
        """ returns a list of strings with the name of each initial condition method
        implemented for the current equation
        """
        methods = [method_name.replace('setInitialCondition', '') for method_name in dir(self)
                  if callable(getattr(self, method_name)) and 'setInitialCondition' in method_name
                   and method_name != 'setInitialCondition' and method_name != 'setInitialConditionZero']
        
        return ['Zero'] + methods

    def getDataFolder(self):
        return DATAFOLDER + self.name + '/'

    def getSavedStatesNames(self):
        # get list of files in data/equation/
        folder = self.getDataFolder()
        _, _, filenames = next(os.walk(folder))

        return [fname.replace('.npz', '') for fname in filenames]

    def Laplace1D(self, x):
        dx = self.getParam('dx')
        laplace = (np.roll(x, 1) + np.roll(x, -1) - 2*x)/(dx*dx)
        if self.boundary_condition == 'neumann':
            laplace[0] = 2 * (x[1] - x[0]) / (dx*dx)
            laplace[-1] = 2 * (x[-2] - x[-1]) / (dx*dx)
        return laplace

    def Derivx(self, x):
        dx = self.getParam('dx')
        derivx = (np.roll(x, -1) - np.roll(x, 1)) / (2*dx)
        if self.boundary_condition =='neumann':
            derivx[0] = 0
            derivx[-1] = 0
        return derivx

    def getBoundaryCondition(self):
        return self.boundary_condition

    def getAllBoundaryConditions(self):
        return ALL_BOUNDARY_CONDITIONS

    def setBoundaryCondition(self, bc):
        self.boundary_condition = bc