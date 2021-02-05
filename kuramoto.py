import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class Kuramoto: ## 1D Kuramoto w global coupling and phase lag
    def __init__(self, N, K, alpha, gamma):
        self.N = N
        self.K = K
        self.alpha = alpha
        self.gamma = gamma

        # initialize frequencies
        #self.lorentzian()
        self.omega = np.zeros(self.N)

        # initialize fourier transform of coupling kernel
        self.initG()

    def initG(self):
        g = np.zeros(self.N) + self.K / (2 * np.pi)
        self.g_ft = np.fft.fft(g)

    def lorentzian(self): #  heterogeneous distribution for the frequencies
        xi = np.random.rand(self.N) * 2 - 1
        y = self.gamma * np.tan(np.pi * xi / 2)
        y[abs(y) > 100] = 0
        self.omega = y

    def rhs(self, t, theta): # Right Hand Side of eq
        e_ft = np.fft.fft(np.exp(theta * 1j))
        zf = np.fft.ifft(self.g_ft * e_ft) * (2 * np.pi / self.N)
        dtheta = self.omega - (np.conj(zf) * np.exp(1j * (theta + self.alpha))).imag

        return dtheta


    def loadState(self, filename):
        data = np.load(filename, allow_pickle=True)
        return data['vals']

    def saveState(self, state, filename):
        np.savez_compressed(filename, vals=state)

    def setInitialCondition(self, vals):
        self.theta0 = vals

    def setInitialConditionFromFile(self, filename, roll=0):
        #self.setInitialCondition(cut(halve(np.roll(self.loadState(filename), roll, axis=1))))
        self.setInitialCondition(self.loadState(filename))

    def setInitialConditionIncoherent(self):
        self.setInitialCondition(np.random.rand(self.N) * 2 * np.pi - np.pi)

    def solve(self, t, dt=1000):
        self.sol = solve_ivp(self.rhs, (0, t[-1]), self.theta0, method='RK45', t_eval=t, max_step=dt).y

    def getState(self, k):
        return (self.sol[:, k] + np.pi) % (2*np.pi) - np.pi

def main():

    #### constantes ####

    N = 256
    K = .01
    alpha = 0
    gamma = 0.01


    #### Inicializar y condicion inicial #####

    k1 = Kuramoto(N, K, alpha, gamma)
    k1.setInitialConditionIncoherent()


    #### Resolver ####

    T = 101
    t = np.linspace(0, 10, T)
    k1.solve(t)



    #### Animar ####

    fig, ax = plt.subplots()

    line, = ax.plot(np.arange(N), k1.theta0, 'o')

    def init():
        line.set_ydata(k1.theta0)
        return [line]

    def animate(i):
        vals = k1.getState(i)
        line.set_ydata(vals)
        return [line]

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=T-1, interval=30, blit=True)

    plt.xlabel(r'Osciladores')
    plt.ylabel(r'$\theta$')
    plt.show()

    #anim.save('Kuramoto_desync.mp4')

if __name__ == '__main__':
    main()