import numpy as np

def fft(dt, y):
    fa = 1.0/dt  # scan frequency
    Y = np.fft.fft(y)
    N = int(len(Y)/2+0.5)
    X = np.linspace(0, fa/2, N, endpoint=True)
    return np.array([X, 2.0*np.abs(Y[:N])/N])
