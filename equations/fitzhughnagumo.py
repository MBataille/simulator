from tkinter import *
from .equation import Equation, Parameter
import numpy as np

class FitzHughNagumo(Equation):

	def __init__(self):
		initParams = { 'dt' : 0.1,
					'dx': 0.5, 
					'I': 1, 
					'a': 1,
					'b': 1,
					'tau': 1}

		Equation.__init__(self, 'FitzHugh-Nagumo', initParams, dim=1, n_fields=2, N=200, fieldNames=['v', 'w'])

	def rhs(self, t, V, W):
		v = self.getCurrentParams()
		dx = v['dx']; I = v['I']; a = v['a']; b = v['b']; tau = v['tau']

		dV = V - V**3 / 3 - W + I + self.Laplace1D(V)
		dW = 1/tau * (V + a - b * W + self.Laplace1D(W))

		return (dV, dW)

	
