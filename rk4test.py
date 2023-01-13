import numpy as np
import matplotlib.pyplot as plt

# exp decay
alpha = -0.1
def rhs(t, X):
    return np.array([X[1], -X[0] + 0.4 * X[0] ** 3 - X[0] ** 5 - 0.1 * X[1] + 2.75 * np.cos(0.7 * t)])

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

ts = np.linspace(0, 10, 1001)
xs = rk4(rhs, (0, 10), np.array([0., 0.]), t_eval=ts)

plt.plot(ts, xs[0, :])
plt.show()