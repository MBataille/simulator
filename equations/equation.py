import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
from tkinter import StringVar
import os

from solvers import INTEGRATION_METHODS

DATAFOLDER = 'data/'

class Parameter:
    def __init__(self, name, var):
        self.name = name
        self.var = var
        self.val = float(var.get())

    def getVal(self):
        try:
            val = float(self.var.get())
            self.val = val
            return self.val
        except ValueError:
            return self.val

def Laplace1D(x, dx, neumann=False):
    laplace = (np.roll(x, 1) + np.roll(x, -1) - 2*x)/(dx*dx)
    if neumann:
        laplace[0] = 2 * (x[1] - x[0]) / (dx*dx)
        laplace[-1] = 2 * (x[-2] - x[-1]) / (dx*dx)
    return laplace

def Derivx(x, dx, neumann=False):
    derivx = (np.roll(x, -1) - np.roll(x, 1)) / (2*dx)
    if neumann:
        derivx[0] = 0
        derivx[-1] = 0
    return derivx

class Equation:
    def __init__(self, name, initParams, dim=1, N=200, isComplex=False, img=None, n_fields=1, fieldNames=['u']):
        self.name = name
        self.dim = dim
        self.img = img
        self.initParams = initParams
        self.parameters = self.createParamsDict(initParams)
        self.isComplex = isComplex
        self.fieldNames = fieldNames

        # initialize w/ some solver
        self.solver = INTEGRATION_METHODS[0]

        self.n_fields = n_fields

        # each field will have N points
        self.Ni = N
        self.N = n_fields * N

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
        dt = self.getParam('dt')
        t = np.linspace(0, 60*dt, 61) # change this
        # print(t, dt, self.initCond)
        with threadpool_limits(limits=1):
            if self.dim == 1:
                self.sol = self.solver.solve(self.rhs, (0, 60*dt), self.initCond, t)
            elif self.dim == 2:
                self.sol = self.solver.solve(self.rhs, (0, dt), self.initCond.reshape(self.N * self.N), t)

    def updateX(self):
        dx = self.parameters['dx'].getVal()
        Ni = round(self.N / self.n_fields)
        X = (Ni-1)*dx/2
        self.x = np.linspace(-X, X, Ni)

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
        self.initParams = self.arraytoInitParams(data['pnames'], data['pvals'])
        self.parameters = self.createParamsDict(self.initParams)

    def arraytoInitParams(self, pnames, pvals):
        initParams = {}
        for i in range(len(pnames)):
            initParams[pnames[i]] = pvals[i]
        return initParams

    def paramsToArray(self):
        v = self.getCurrentParams()
        pnames = []
        pvals = []
        for p in v:
            pnames.append(p)
            pvals.append(v[p])
        return np.array(pnames), np.array(pvals)

    def getFields(self, k):
        Ni = round(self.N / self.n_fields)
        return [self.sol[i * Ni : (i + 1) * Ni, k] for i in range(self.n_fields)]

    def getField(self, i, k):
        if self.n_fields == 1:
            return self.getState(k)
        Ni = self.Ni
        return self.sol[i * Ni:(i + 1) * Ni, k]

    def getInitCondFields(self):
        Ni = self.Ni
        return [self.initCond[i * Ni : (i + 1) * Ni] for i in range(self.n_fields)]

    def setInitCondFields(self, fields):
        Ni = self.Ni
        for i in range(len(fields)):
            self.initCond[i * Ni: (i+1) * Ni] = fields[i]

    def getInitialConditions(self):
        """ returns a list of strings with the name of each initial condition method
        implemented for the current equation
        """
        methods = [method_name.replace('setInitialCondition', '') for method_name in dir(self)
                  if callable(getattr(self, method_name)) and 'setInitialCondition' in method_name
                   and 'setInitialCondition' != method_name]
        
        return methods

    def getDataFolder(self):
        return DATAFOLDER + self.name + '/'

    def getSavedStatesNames(self):
        # get list of files in data/equation/
        folder = self.getDataFolder()
        _, _, filenames = next(os.walk(folder))

        return [fname.replace('.npz', '') for fname in filenames]

