from .equation import Equation, Parameter
import numpy as np

class Duffing(Equation):

	def __init__(self):
		initParams = { 'dt' : 0.1,
					'dx': 1, 
					'alpha': 0.4,
					'mu': 0.1,
					'gamma' : 2.75,
					'omega' : 0.7,
					'kappa' : 0.42, 
					'delta': 0}

		Equation.__init__(self, 'Duffing', initParams, dim=1, n_fields=2, N=200, fieldNames=['x', 'dx/dt'])

	def rhs(self, t, X, dXdt):
		v = self.getCurrentParams()

		dx = v['dx']; alpha = v['alpha']; kappa = v['kappa']; delta = v['delta']; mu = v['mu']; 
		gamma = v['gamma']; omega = v['omega']

		X_nm = np.roll(X, 1) # x_n-1
		X_np = np.roll(X, -1) # x_n+1
		d2Xdt2 = X * (-1 - 2*kappa) + alpha * X ** 3 - X ** 5 - mu * dXdt \
						+ gamma * np.cos(omega * t) + (kappa + delta) * X_np + (kappa - delta) * X_nm

		return (dXdt, d2Xdt2)

	def setInitialConditionIncoherent(self):
		X = 0.3 * np.random.normal(size=self.Ni)
		dXdt = np.zeros(self.Ni)
		self.setInitialCondition((X, dXdt))

	def getAllBoundaryConditions(self):
		return ['periodic']