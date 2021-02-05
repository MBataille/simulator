import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class GinzburgLandau:
	def __init__(self, N, L, b, c):
		self.N = N
		self.L = L
		self.delta = L/(N-1)
		self.b = b
		self.c = c

	def rhs(self, t, A):
		a = A.reshape(self.N, self.N)
		laplace = 1/(self.delta**2) * (np.roll(a, 1, axis=1) + np.roll(a,-1,axis=1) \
			+ np.roll(a, 1, axis=0) + np.roll(a,-1,axis=0) - 4 * a)
		# no flux
		# laplace[0, :] = 2 / (self.delta**2) * (a[1, :] - a[0, :])
		# laplace[self.N - 1, :] = 2 / (self.delta**2) * (a[self.N - 2, :] - a[self.N - 1, :])
		# laplace[:, 0] = 2 / (self.delta**2) * (a[:, 1] - a[:, 0])
		# laplace[:, self.N - 1] = 2 / (self.delta**2) * (a[:, self.N - 2] - a[:, self.N - 1])

		return (a + (1 + 1j * self.b) * laplace - (1 + 1j * self.c) * a * a * np.conj(a)).reshape(self.N * self.N)

	def loadState(self, filename):
		data = np.load(filename, allow_pickle=True)
		return data['vals']

	def saveState(self, state, filename):
		np.savez_compressed(filename, vals=state)

	def setInitialCondition(self, vals):
		self.a0 = vals

	def setInitialConditionFromFile(self, filename, roll=0):
		#self.setInitialCondition(cut(halve(np.roll(self.loadState(filename), roll, axis=1))))
		self.setInitialCondition(np.roll(halve(self.loadState(filename)), roll, axis=1))

	def setInitialConditionIncoherent(self):
		self.setInitialCondition(0.99 * np.exp(1j * np.random.rand(self.N, self.N) * 2 * np.pi - np.pi))

	def solve(self, t, dt=1000):
		self.sol = solve_ivp(self.rhs, (0, t[-1]), self.a0.reshape(self.N * self.N), method='RK45', t_eval=t, max_step=dt).y

	def getState(self, k):
		return np.reshape(self.sol[:, k], (self.N, self.N)) 

def cut(data):
	N = len(data)
	#data2 = np.zeros((int(N/2), int(N/2)))
	d2 = data[:, :int(N/2)]
	return np.append(d2, d2[:,::-1], axis=1)

def halve(data):
	data2 = np.zeros((int(len(data)/2), int(len(data)/2)))
	data2 = data[::2, ::2]
	data2 = data[::2,1::2]
	data2 = data[1::2, ::2]
	data2 = data[1::2,1::2]
	return data2

# Nicolau, b = c1, c = -c3
# Aranson kramer, b = 0, c = 1.5
# Nicolau, b = 1.5, c = -0.77
t1 = GinzburgLandau(128, 1, 2, -0.7)
#t1.setInitialConditionFromFile('gl_test.npz', roll=60)

#t1.setInitialConditionFromFile('ginzburglandau/s001a092.npz')
t1.setInitialConditionIncoherent()

plt.imshow(np.angle(t1.a0), cmap='hsv')
plt.show()
t = np.linspace(0, .1, 101)
t1.solve(t)
#plt.imshow(np.angle(t1.getState(-1)), cmap='hsv')
#plt.show()

fig = plt.figure()
#ax = plt.axes(xlim=(0, 256), ylim=(0, 256))

im = plt.imshow(np.angle(t1.getState(0)), cmap='hsv')

def init():
	return [im]

def animate(i):
	vals = t1.getState(i)
	im.set_data(np.angle(vals))
#	im.set_title('t = {}'.format(i))
	return [im]

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=30, blit=True)

plt.show()