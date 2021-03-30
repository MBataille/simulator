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
		dx = self.parameters['dx'].val.get()
		I = self.parameters['I'].val.get()
		a = self.parameters['a'].val.get()
		b = self.parameters['b'].val.get()
		tau = self.parameters['tau'].val.get()
		#B = self.parameters['B'].val.get()

		# assuming X and Y have N/2 points each
		V = x[:round(self.N/2)]
		Omega = x[round(self.N/2):]

		dV = V - V**3 / 3 - Omega + I + Laplace1D(V, dx, neumann=True)
		dOmega = 1/tau * (V + a - b * Omega + Laplace1D(Omega, dx, neumann=True))

		return np.append(dV, dOmega)

	
