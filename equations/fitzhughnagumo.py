from tkinter import *
from .equation import Equation, Laplace1D, Parameter
import numpy as np

class FitzHughNagumo(Equation):

	def __init__(self):
		self.initParams = { 'dt' : 0.1,
					'dx': 0.5, 
					'I': 1, 
					'a': 1,
					'b': 1,
					'tau': 1}

		Equation.__init__(self, 'FHN', 1, self.createParamsDict(self.initParams), n_fields=2, N=400, fieldNames=['v', 'w'])

	def rhs(self, t, x):
		v = self.getCurrentParams()
		dx = v['dx']; I = v['I']; a = v['a']; b = v['b']; tau = v['tau']

		V = x[:self.Ni]
		Omega = x[self.Ni:]

		dV = V - V**3 / 3 - Omega + I + Laplace1D(V, dx, neumann=True)
		dOmega = 1/tau * (V + a - b * Omega + Laplace1D(Omega, dx, neumann=True))

		return np.append(dV, dOmega)

	
