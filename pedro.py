from tkinter import *
from eqn import Equation, Laplace1D, Parameter
import numpy as np

class Pedro(Equation):

	def __init__(self):
		self.initParams = { 'dt' : 0.1,
							'dx': 0.5, 
							'eps': 1, 
							'eta': 0,
							'gamma': 0.5,
							'k': 0.07,
							'omega': 2}

		Equation.__init__(self, 'KMOEP', 1, self.createParamsDict(self.initParams))
		
	def rhs(self, t, x): # CB libre
		dx = self.parameters['dx'].val.get()
		eps = self.parameters['eps'].val.get()
		eta = self.parameters['eta'].val.get()
		gamma = self.parameters['gamma'].val.get()
		k = self.parameters['k'].val.get()
		omega = self.parameters['omega'].val.get()

		y = eta + eps*x - x**3  + Laplace1D(x, dx, neumann=True) + gamma * np.cos(k*x) * np.sin(omega*t)
		return y

	def setInitialConditionKink(self):
		self.updateX()
		eps = self.parameters['eps'].val.get()
		self.initCond = np.sqrt(eps) * np.tanh(self.x * np.sqrt(eps/2))