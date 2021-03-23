from tkinter import *
from .equation import Equation, Laplace1D, Parameter
import numpy as np

class Brusselator(Equation):

	def __init__(self):
		self.initParams = { 'dt' : 0.1,
					'dx': 0.5, 
					'A': 1, 
					'B': 3}

		Equation.__init__(self, 'Brusselator', 1, self.createParamsDict(self.initParams), n_fields=2, N=400, fieldNames=['x', 'y'])

	def rhs(self, t, x):
		dx = self.parameters['dx'].val.get()
		A = self.parameters['A'].val.get()
		B = self.parameters['B'].val.get()

		# assuming X and Y have N/2 points each
		X = x[:round(self.N/2)]
		Y = x[round(self.N/2):]

		dX = A + X*(X*Y - B - 1) + Laplace1D(X, dx, neumann=True)
		dY = X*(B -X*Y) + Laplace1D(Y, dx, neumann=True)

		return np.append(dX, dY)

	
