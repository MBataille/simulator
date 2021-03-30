from tkinter import *
from .equation import Equation, Laplace1D, Parameter
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
		dx = self.getParam('dx')
		eps = self.getParam('eps')
		eta = self.getParam('eta')
		gamma = self.getParam('gamma')
		k = self.getParam('k')
		omega = self.getParam('omega')

		y = eta + eps*x - x**3  + Laplace1D(x, dx, neumann=True) + gamma * np.cos(k*self.x) * np.sin(omega*t)
		return y

	def setInitialConditionKink(self):
		self.updateX()
		eps = self.getParam('eps')
		self.initCond = np.sqrt(eps) * np.tanh(self.x * np.sqrt(eps/2))