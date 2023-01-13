from equations.equation import Equation, Parameter
import numpy as np

class Rene(Equation):

    def __init__(self):
        initParams = { 'dt' : 0.1,
                            'dx': 0.5, 
                            'eps': 1, 
                            'eta': 0,
                            'D': 1,
                            'delta': 0.5
                        }

        Equation.__init__(self, 'Rene', initParams, dim=1)
        
    def rhs(self, t, u): 
        v = self.getCurrentParams()
        eta = v['eta']; eps = v['eps']; D = v['D']; dx = v['dx']; delta = v['delta']
        N = u.shape[0]

        u_fourier = np.fft.fft(u)
        k = np.fft.fftfreq(N) * 2 * np.pi / dx

        u_tilda = np.fft.ifft(np.exp(-1j * k * delta) * u_fourier).real
        u_xx = np.fft.ifft(- k * k * u_fourier).real

        y = eta + eps * u_tilda - u_tilda ** 3 + D * u_xx
        return y

    def setInitialConditionKink(self):
        self.updateX()
        eps = self.getParam('eps')
        self.initCond = np.sqrt(eps) * np.tanh(self.x * np.sqrt(eps/2))