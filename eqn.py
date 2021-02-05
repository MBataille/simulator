import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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
	def __init__(self, name, dim, parameters, N=101, isComplex=False, img=None):
		self.name = name
		self.dim = dim
		self.img = img
		self.parameters = parameters
		self.N = N
		self.isComplex = isComplex

		# parameters[0] = dt

	def setSolver(self, solver):
		self.solver = solver

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
		dtype = 'complex128' if isComplex else 'float64'
		if self.dim == 1:
			self.initCond = np.zeros(N, dtype=dtype)
		if self.dim == 2:
			self.initCond = np.zeros((N, N), dtype=dtype)

	def solve(self):
		dt = self.getParam('dt')
		t = np.linspace(0, 10*dt, 11)
		# print(t, dt, self.initCond)
		if self.dim == 1:
			self.sol = solve_ivp(self.rhs, (0, 10*dt), self.initCond, method='RK45', t_eval=t).y
		elif self.dim == 2:
			self.sol = solve_ivp(self.rhs, (0, dt), self.initCond.reshape(self.N * self.N), method='RK45', t_eval=t).y

	def updateX(self):
		dx = self.parameters['dx'].val.get()
		X = (self.N-1)*dx/2
		self.x = np.linspace(-X, X, self.N + 2)

	def getState(self, k):
		if self.dim == 1:
			return self.sol[:, k]
		elif self.dim == 2:
			return np.reshape(self.sol[:, k], (self.N, self.N)) 

	def getParam(self, p):
		return self.parameters[p].val.get()