import numpy as np
from scipy.integrate import solve_ivp

class IntMethod:
    def __init__(self, code, name, description, func=None):
        self.code = code
        self.name = name
        self.description = description
        self.func = func

    def solve(self, rhs, t_span, initCond, t_eval):
        if self.func is None: # then calls corresponding method from Scipy
            return solve_ivp(rhs, t_span, initCond, method=self.code, t_eval=t_eval).y
        else:
            return self.func(rhs, t_span, initCond, t_eval=t_eval)

#### Custom integration methods

# Runge-Kutta 4 (for perf?)

def rk4(rhs, t_span, initCond, t_eval=None):
    if t_eval is None:
        dt = 0.1
    else:
        dt = t_eval[1] - t_eval[0]
    ts = np.linspace(t_span[0], t_span[1], int((t_span[1] - t_span[0]) / (dt) + 1))
    xs = np.zeros((len(initCond), len(ts)), dtype=type(initCond[0]))

    x = initCond
    xs[:, 0] = initCond
    for i, t in enumerate(ts[1:]):
        t = ts[i+1]
        k1 = dt * rhs(t, x)
        k2 = dt * rhs(t + 0.5 * dt, x + 0.5 * k1)
        k3 = dt * rhs(t + 0.5 * dt, x + 0.5 * k2)
        k4 = dt * rhs(t + dt, x + k3)

        xs[:, i+1] = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x = xs[:, i+1]
    return xs



#### Create available integration methods:

## 1. Descriptions: (the tedious and long part)

descriptions = {'RK45': ' Explicit Runge-Kutta method of order 5(4).',
                'RK23': 'Explicit Runge-Kutta method of order 3(2). ',
                'DOP853': 'Explicit Runge-Kutta method of order 8.',
                'Radau': 'Implicit Runge-Kutta method of the Radau IIA family of order 5.\nCan\'t be applied on complex domain',
                'BDF': 'Implicit multi-step variable-order (1 to 5) method\n based on a backwarddifferentiation formula for the derivative approximation. ',
                'LSODA': 'Adams/BDF method with automatic stiffness detection and switching.\nCan\'t be applied on complex domain',
                'm-RK4': 'Explicit Runge-Kutta method of order 5(4) implemented by Martin.'}


## Rk45
rk45 = IntMethod('RK45', 'Exp Runge-Kutta 5th Order', descriptions['RK45'])
rk23 = IntMethod('RK23', 'Exp Runge-Kutta 3th Order', descriptions['RK23'])
dop853 = IntMethod('DOP853', 'Exp Runge-Kutta 8th Order', descriptions['DOP853'])
radau = IntMethod('Radau', 'Imp Runge-Kutta 5th Order', descriptions['Radau'])
bdf = IntMethod('BDF', 'Imp variable-order', descriptions['BDF'])
lsoda = IntMethod('LSODA', 'Adams/BDF method', descriptions['LSODA'])
#euler

customrk4 = IntMethod('m-RK4', 'Exp custom made Kunge-Kutta 4th Order', descriptions['m-RK4'], rk4)

INTEGRATION_METHODS = [rk45, rk23, dop853, radau, bdf, lsoda, customrk4]