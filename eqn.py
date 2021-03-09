import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
from tkinter import * 

from solvers import INTEGRATION_METHODS

class Parameter:
	def __init__(self, name, val):
		self.name = name
		self.val = val

def Laplace1D(x, dx, neumann=False):
	laplace = (np.roll(x, 1) + np.roll(x, -1) - 2*x)/(dx*dx)
	if neumann:
		laplace[0] = 2 * (x[1] - x[0]) / (dx*dx)
		laplace[-1] = 2 * (x[-2] - x[-1]) / (dx*dx)
	return laplace

class Equation:
	def __init__(self, name, dim, parameters, N=200, isComplex=False, img=None, n_fields=1, fieldNames=['u']):
		self.name = name
		self.dim = dim
		self.img = img
		self.parameters = parameters
		self.N = N
		self.isComplex = isComplex
		self.fieldNames = fieldNames

		# initialize w/ some solver
		self.solver = INTEGRATION_METHODS[0]

		self.n_fields = n_fields

		self.Ni = round(N / n_fields)
		# parameters[0] = dt

	def setSolver(self, solver):
		self.solver = solver

	def createParamsDict(self, initParams):
		params = {}
		for p in initParams:
			pvar = DoubleVar() # create var
			pvar.set(initParams[p]) # assign initial value
			params[p] = Parameter(p, pvar) # create parameter
		return params

	def loadState(self, filename):
		data = np.load(filename, allow_pickle=True)
		return data['vals']

	def saveState(self, state, filename):
		np.savez_compressed(filename, vals=state)

	def setInitialCondition(self, vals):
		self.initCond = vals

	def setInitialConditionFromFile(self, filename):
		self.setInitialCondition(self.loadState(filename))

	def setInitialConditionIncoherent(self):
		self.setInitialCondition(0.99 * np.exp(1j * np.random.rand(self.N, self.N) * 2 * np.pi - np.pi))

	def setInitialConditionZero(self):
		dtype = 'complex128' if self.isComplex else 'float64'
		if self.dim == 1:
			self.initCond = np.zeros(self.N, dtype=dtype)
		if self.dim == 2:
			self.initCond = np.zeros((self.N, self.N), dtype=dtype)

	# def solve(self):
		# dt = self.getParam('dt')
		# t = np.linspace(0, 60*dt, 61)
		# # print(t, dt, self.initCond)
		# with threadpool_limits(limits=1):
		# 	if self.dim == 1:
		# 		self.sol = solve_ivp(self.rhs, (0, 60*dt), self.initCond, method='RK45', t_eval=t).y
		# 	elif self.dim == 2:
		# 		self.sol = solve_ivp(self.rhs, (0, dt), self.initCond.reshape(self.N * self.N), method='RK45', t_eval=t).y

	def solve(self):
		dt = self.getParam('dt')
		t = np.linspace(0, 60*dt, 61)
		# print(t, dt, self.initCond)
		with threadpool_limits(limits=1):
			if self.dim == 1:
				self.sol = self.solver.solve(self.rhs, (0, 60*dt), self.initCond, t)
			elif self.dim == 2:
				self.sol = self.solver.solve(self.rhs, (0, dt), self.initCond.reshape(self.N * self.N), t)

	def updateX(self):
		dx = self.parameters['dx'].val.get()
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

	def getFields(self, k):
		Ni = round(self.N / self.n_fields)
		return [self.sol[i * Ni : (i + 1) * Ni, k] for i in range(self.n_fields)]

	def getField(self, i, k):
		if self.n_fields == 1:
			return self.getState(k)
		Ni = round(self.N / self.n_fields)
		return self.sol[i * Ni:(i + 1) * Ni, k]

	def getInitCondFields(self):
		Ni = round(self.N / self.n_fields)
		return [self.initCond[i * Ni : (i + 1) * Ni] for i in range(self.n_fields)]

	def setInitCondFields(self, fields):
		Ni = round(self.N / self.n_fields)
		for i in range(len(fields)):
			self.initCond[i * Ni: (i+1) * Ni] = fields[i]

	def getParam(self, p):
		return self.parameters[p].val.get()