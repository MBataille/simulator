from tkinter import *
from .equation import Equation, Laplace1D, Parameter
import numpy as np

class Brusselator(Equation):

	def __init__(self):
		initParams = { 'dt' : 0.1,
					'dx': 0.5, 
					'A': 1, 
					'B': 3}

		Equation.__init__(self, 'Brusselator', initParams, dim=1, n_fields=2, N=200, fieldNames=['x', 'y'])

	def rhs(self, t, x):
		v = self.getCurrentParams()

		dx = v['dx']; A = v['A']; B = v['B']

		# assuming X and Y have N/2 points each
		X = x[:round(self.N/2)]
		Y = x[round(self.N/2):]

		dX = A + X*(X*Y - B - 1) + Laplace1D(X, dx, neumann=True)
		dY = X*(B -X*Y) + Laplace1D(Y, dx, neumann=True)

		return np.append(dX, dY)

	
