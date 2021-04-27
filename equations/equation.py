import numpy as np
from matplotlib import animation
# import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
from tkinter import StringVar
import os

from solvers import INTEGRATION_METHODS

DATAFOLDER = 'data/'
ALL_BOUNDARY_CONDITIONS = 'neumann periodic'.split(' ')
RECORD_FILENAME = 'state'
T0_UPPER_BOUND = 1e6
SOLVE_EVERY_TI = 60 # solve every 60 time steps

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
        self.initParams['noise'] = 0 # add noise to params
        self.parameters = self.createParamsDict(self.initParams)
        self.isComplex = isComplex
        self.fieldNames = fieldNames
  
        # init tick // Perhaps this should go outside?
        self.recording = False
        self.init_tick()

        self.boundary_condition = 'neumann'

        # initialize w/ some solver
        self.solver = INTEGRATION_METHODS[0]

        self.n_fields = n_fields

        # each field will have N points
        self.setNi(N)

    def init_tick(self):
        self.k_sol = 0
        self.t0 = 0
        self.last_params = self.getCurrentParams()

    def tick(self, newInitCondFields=None): #  make 1 time step
        new_params = self.getCurrentParams()
        if self.k_sol % SOLVE_EVERY_TI == 0 or new_params != self.last_params or newInitCondFields is not None:
            # just in case
            self.updateX()
            if newInitCondFields is not None:
                self.setInitCondFields(newInitCondFields)
            self.solve_cycle()
        
        self.k_sol += 1
        
        if self.recording:
            self.saveRecord()
        
        self.last_params = new_params

        # may not be necessary
        # return self.getFields(self.k_sol)

    def getCurrentFields(self): # return Fields at k = k_sol
        return self.getFields(self.k_sol)

    def solve_cycle(self):
        self.solve()
        self.initCond = self.sol[:, -1]
        self.k_sol = 0

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
        #print(self.getCurrentParams())
        if self.n_fields == 1:
            if self.getParam('noise') == 0:
                return self.rhs(t, u)
            else:
                return self.rhs(t, u) + self.GaussianWhiteNoise()
        else:
            fields = self.extractFields(u)
            d_fields_dt = self.rhs(t, *fields)
            if self.getParam('noise') != 0:
                for dfield_dt in d_fields_dt:
                    dfield_dt += self.GaussianWhiteNoise()
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

    def isFolder(self, foldername):
        folder = self.getDataFolder()
        return os.path.exists(folder + foldername)

    def startRecording(self, foldername):
        # assuming this folder doesn't exists yet
        path = self.getDataFolder() + foldername
        if not self.isFolder(foldername):
            os.mkdir(path)
        self.k_recording = 0
        self.recording  = True
        self.foldername_recording = foldername + '/'

    def getRecordStateName(self, k):
        return RECORD_FILENAME + f'_{k}' 

    def saveRecord(self):
        self.k_recording += 1
        filename = self.foldername_recording + self.getRecordStateName(self.k_recording)
        self.saveState(self.k_sol, filename)

    def stopRecording(self):
        self.k_recording = None
        self.path_recording = None
        self.recording = False

    def saveState(self, k, filename):
        folder = self.getDataFolder()
        _vals = self.sol[:, k]
        _pnames, _pvals = self.paramsToArray()
        np.savez_compressed(folder + filename + '.npz', vals=_vals, pnames=_pnames, pvals=_pvals)

    def loadState(self, filename):
        folder = self.getDataFolder()

        if os.path.isdir(folder + filename):
            self.loadState(filename + '/' + self.getRecordStateName(1))
            return

        data = np.load(folder + filename + '.npz')

        self.setInitialCondition(data['vals'])
        self.Ni = int(data['pvals'][0])
        self.N = int(self.Ni * self.n_fields)
        self.initParams = self.arraytoInitParams(data['pnames'][1:], data['pvals'][1:])
        self.initParams['noise'] = 0
        self.parameters = self.createParamsDict(self.initParams)

    def arraytoInitParams(self, pnames, pvals):
        initParams = {}
        for i in range(len(pnames)):
            if pnames[i] in self.initParams:
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
        """ Sets intial condition from a list of fields = [np.array]
        """
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
        """ Returns list of saved (files) and recorded (folders) states in data/equation_name/
        """
        folder = self.getDataFolder()
        _, foldernames, filenames = next(os.walk(folder))

        return [fname.replace('.npz', '') for fname in filenames] + foldernames

    def get_xlims(self):
        """ Returns min and max of the x axis.
        """
        return self.x[0], self.x[-1]

    def get_field_lims(self, field_indx):
        """ Returns min and max of selected field at k = k_sol (current state)
        """
        field = self.getField(field_indx, self.k_sol)
        return field.min(), field.max()
    
    def get_fields_lims(self, indices=None):
        if indices is None: indices = range(self.n_fields)
        fields = self.getFields(self.k_sol)
        fmin, fmax = np.inf, -np.inf
        for i in indices:
            fm, fM = fields[i].min(), fields[i].max()
            if fm < fmin: fmin = fm
            if fM > fmax: fmax = fM
        return fmin, fmax


##### Operators

    def GaussianWhiteNoise(self):
        eta = self.getParam('noise')
        return np.sqrt(eta) * np.random.normal(size=self.getN())

    def BiLaplace1D(self, x): # only with priodic B.C.
        dx = self.getParam('dx')
        r = np.roll
        bilaplace = (r(x, -2) - 4 * r(x, -1) + 6 * x - 4 * r(x, 1) + r(x, 2))/(dx**4)
        return bilaplace

    def Laplace1D(self, x):
        """ Laplacian operator in 1D. Supports periodic and neumann BCs
            Params: x (np array)s
        """
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

    def getN(self):
        return self.Ni