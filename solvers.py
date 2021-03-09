import numpy as np
from scipy.integrate import solve_ivp

class IntMethod:
	def __init__(self, code, name, description, fromScipy):
		self.code = code
		self.name = name
		self.description = description
		self.fromScipy = fromScipy

	def solve(self, rhs, t_span, initCond, t_eval):
		if self.fromScipy:
			return solve_ivp(rhs, t_span, initCond, method=self.code, t_eval=t_eval).y
		else:
			print('Integration Method not implemented yet')
			### Implement a new integration method.

#### Create available integration methods:

## 1. Descriptions: (the tedious and long part)

descriptions = {'RK45': ' Explicit Runge-Kutta method of order 5(4).',
				'RK23': 'Explicit Runge-Kutta method of order 3(2). ',
				'DOP853': 'Explicit Runge-Kutta method of order 8.',
				'Radau': 'Implicit Runge-Kutta method of the Radau IIA family of order 5.\nCan\'t be applied on complex domain',
				'BDF': 'Implicit multi-step variable-order (1 to 5) method\n based on a backwarddifferentiation formula for the derivative approximation. ',
				'LSODA': 'Adams/BDF method with automatic stiffness detection and switching.\nCan\'t be applied on complex domain'}


## Rk45
rk45 = IntMethod('RK45', 'Exp Runge-Kutta 5th Order', descriptions['RK45'],True)
rk23 = IntMethod('RK23', 'Exp Runge-Kutta 3th Order', descriptions['RK23'],True)
dop853 = IntMethod('DOP853', 'Exp Runge-Kutta 8th Order', descriptions['DOP853'],True)
radau = IntMethod('Radau', 'Imp Runge-Kutta 5th Order', descriptions['Radau'],True)
bdf = IntMethod('BDF', 'Imp variable-order', descriptions['BDF'],True)
lsoda = IntMethod('LSODA', 'Adams/BDF method', descriptions['LSODA'],True)
#euler

INTEGRATION_METHODS = [rk45, rk23, dop853, radau, bdf, lsoda]