from .equation import Equation, Parameter
import numpy as np

class ComplexGinzburgLandau(Equation):

	def __init__(self):
		initParams = { 'dt' : 0.1,
					'dx': 1, 
					'alpha': 0,
					'beta': -3,
					'eta': 0}

		Equation.__init__(self, 'ComplexGinzburgLandau', initParams, dim=1, n_fields=2, N=200, fieldNames=['Re A', 'Im A'])

	def rhs(self, t, Re, Im):
		v = self.getCurrentParams()

		A = Re + 1j * Im
		dx = v['dx']; alpha = v['alpha']; beta = v['beta']

		dA = A + (1 + 1j*alpha) * self.Laplace1D(A) - (1 +  1j * beta) * np.conj(A) * A * A
		return dA.real + self.GaussianWhiteNoise(), dA.imag + self.GaussianWhiteNoise()

