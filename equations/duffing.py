from .equation import Equation, Parameter, SOLVE_EVERY_TI
from scipy.signal import hilbert
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
		
		initRange = { 'dt' : (0, 0.5),
					'dx': (0, 1), 
					'alpha': (0,1),
					'mu': (0, 0.5),
					'gamma' : (0, 5),
					'omega' : (0, 1),
					'kappa' : (0, 1), 
					'delta': (0, 1)}

		Equation.__init__(self, 'Duffing', initParams, dim=1, n_fields=2, N=200, fieldNames=['x', 'dx/dt'], auxFieldNames=['theta', 'abs z'], initRange=initRange)

	def initAuxFields(self):
		abszs = np.zeros((SOLVE_EVERY_TI, self.getN()))
		for k in range(SOLVE_EVERY_TI):
			theta, absz = self.getAuxFields(*self.getFields(k), calc_kura=True)
			abszs[k, :] = absz
		self.abszs = np.abs(hilbert(abszs, axis=0))

	def rhs(self, t, X, dXdt):
		v = self.getCurrentParams()

		dx = v['dx']; alpha = v['alpha']; kappa = v['kappa']; delta = v['delta']; mu = v['mu']; 
		gamma = v['gamma']; omega = v['omega']

		X_nm = np.roll(X, 1) # x_n-1
		X_np = np.roll(X, -1) # x_n+1

		if self.boundary_condition == 'neumann':
			X_nm[0] = X[1]
			X_np[-1] = X[-2]

		d2Xdt2 = X * (-1 - 2*kappa) + alpha * X ** 3 - X ** 5 - mu * dXdt \
						+ gamma * np.cos(omega * t) + (kappa + delta) * X_np + (kappa - delta) * X_nm

		return (dXdt, d2Xdt2)

	def kuramoto_local_order_param(self, theta):
		N = self.getN()
		absz = np.zeros(N)
		R = int(np.sqrt(N)/4)
		for i in range(N):
			thetas = theta[i-R:i+R]
			if self.boundary_condition == 'periodic':
				if i < R:
					thetas = np.append(theta[N-1+i-R:],theta[:i+R]) 
			if self.boundary_condition == 'neumann':
				if i < R:
					thetas = theta[:i+R] 
			if len(thetas) == 0: print(i, R)
			absz[i] = np.abs(np.exp(1j * thetas).sum())/(len(thetas))
		return absz

	def getAuxFields(self, X, dXdt, calc_kura=False):
		theta = np.arctan2(dXdt, X)
		if calc_kura:
			absz = self.kuramoto_local_order_param(theta)
		else:
			try:
				absz = self.abszs[self.k_sol-1]
			except AttributeError:
				absz = self.kuramoto_local_order_param(theta)
		return theta, absz

	def setInitialConditionIncoherent(self):
		X = 0.3 * np.random.normal(size=self.Ni)
		dXdt = np.zeros(self.Ni)
		self.setInitialCondition((X, dXdt))

	def setInitialConditionGaussian(self):
		dXdt = np.zeros(self.Ni)
		x = self.getX()
		A = 4
		sigma = 15
		X = A * np.exp(-x**2/(2*sigma**2))
		self.setInitialCondition((X, dXdt))

	def getMarkers(self, X, dXdt, theta, absz):
		pass