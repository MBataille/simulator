from tkinter import *
from .equation import Equation, Laplace1D, Parameter
import numpy as np

class EqKink(Equation):
	'''params: dt, dx, alfa'''
	def __init__(self):
		initParams = { 'dt' : 0.1,
							'dx': 0.5, 
							'eps': 1, 
							'alpha': 0 }

		Equation.__init__(self, 'Eqntest', initParams, dim=1)
		

	def rhs(self, t, x): # CB periodica
		v = self.getCurrentParams()
		dx = v['dx']; eps = v['eps']; alpha = v['alpha']

		y = alpha + eps*x - x**3  + Laplace1D(x, dx, neumann=True)
		return y

	def setInitialConditionKink(self):
		self.updateX()
		eps = self.getParam('eps')
		self.initCond = np.sqrt(eps) * np.tanh(self.x * np.sqrt(eps/2))