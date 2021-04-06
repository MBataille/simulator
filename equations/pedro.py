from tkinter import *
from .equation import Equation, Laplace1D, Parameter
import numpy as np

class Pedro(Equation):

	def __init__(self):
		initParams = { 'dt' : 0.1,
							'dx': 0.5, 
							'eps': 1, 
							'eta': 0,
							'gamma': 0.5,
							'k': 0.07,
							'omega': 2}

		Equation.__init__(self, 'Pedro', initParams, dim=1)
		
	def rhs(self, t, x): 
		v = self.getCurrentParams()
		dx = v['dx']; eps = v['eps']; eta = v['eta']; gamma = v['gamma']; k = v['k']; omega = v['omega']

		y = eta + eps*x - x**3  + Laplace1D(x, dx, neumann=True) + gamma * np.cos(k*self.x) * np.sin(omega*t)
		return y

	def setInitialConditionKink(self):
		self.updateX()
		eps = self.getParam('eps')
		self.initCond = np.sqrt(eps) * np.tanh(self.x * np.sqrt(eps/2))