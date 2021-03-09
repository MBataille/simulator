from tkinter import *
from eqn import Equation, Laplace1D, Parameter
import numpy as np

class EqKink(Equation):
	'''params: dt, dx, alfa'''
	def __init__(self):
		self.initParams = { 'dt' : 0.1,
							'dx': 0.5, 
							'eps': 1, 
							'alpha': 0 }
		
		self.params = {}

		for p in self.initParams:
			pvar = DoubleVar() # create var
			pvar.set(self.initParams[p]) # assign initial value
			self.params[p] = Parameter(p, pvar) # create parameter

		Equation.__init__(self, 'Test', 1, self.params)
		

	def rhs(self, t, x): # CB periodica
		dx = self.parameters['dx'].val.get()
		eps = self.parameters['eps'].val.get()
		alpha = self.parameters['alpha'].val.get()

		y = alpha + eps*x - x**3  + Laplace1D(x, dx, neumann=True)
		return y

	def setInitialConditionKink(self):
		self.updateX()
		eps = self.parameters['eps'].val.get()
		self.initCond = np.sqrt(eps) * np.tanh(self.x * np.sqrt(eps/2))