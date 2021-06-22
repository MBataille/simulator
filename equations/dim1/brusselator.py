from equations.equation import Equation, Parameter
import numpy as np

### 
# In the case of many fields you will have to:
#	1. specify n_fieds and fieldNames
#
#	2. in the rhs function, recieve each field
#	as an argument and return them as a tuple
#	(separated by commas).
#
#	3. (optional) if you want to define special
#	initial conditions, then 



class Brusselator(Equation):

	def __init__(self):
		initParams = { 'dt' : 0.1,
					'dx': 0.5, 
					'A': 1, 
					'B': 3}

		Equation.__init__(self, 'Brusselator', initParams, dim=1, n_fields=2, N=200, fieldNames=['X', 'Y'])

	def rhs(self, t, X, Y):
		v = self.getCurrentParams()
		dx = v['dx']; A = v['A']; B = v['B']

		dXdt = A + X*(X*Y - B - 1) + self.Laplace1D(X) 
		dYdt = X*(B -X*Y) + self.Laplace1D(Y) 

		return dXdt, dYdt

	def setInitialConditionGaussians(self):
		x = self.getX()
		N = len(x)

		# X is a Gaussian centered in 0
		# Y is just an array of zeros
		X = np.exp(-x**2/2)
		Y = np.zeros(N)

		fields = (X, Y)
		self.setInitialCondition(fields)
