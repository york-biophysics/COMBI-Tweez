import json
# from os import getcwd, listdir, makedirs, path

import numpy as np
from lmfit import Minimizer, Parameters
from numpy.fft import rfft, rfftfreq
from scipy.signal import periodogram

from brownian import k_brown
from functions import anti_aliased_lorentzian, lorentzian
from settings import D, Fs, axl_correction, gamma, lat_correction

use_moving_average = True
anti_aliasing = False


def moving_average(data, window=10):
    extended_data = np.hstack([[data[0]] * (window - 1), data])
    weightings = np.repeat(1.0, window) / window
    return np.convolve(extended_data, weightings)[window-1:-(window-1)]


def psd(x, Fs=Fs):
    # Calculate the 1D Fourier Transform
    # BUT THIS is just the +ve part of the two-sided
    # f, p = rfftfreq(x.size, d=1./Fs), np.absolute(rfft(x))**2
    f, p = periodogram(x, fs=Fs, return_onesided=True, scaling='density')
    return np.array([f, p])

def fit(f, p, axis, cutoffs=None, antialiasing=None):
    f = np.array(f)
    p = np.array(p)
    # Cutoff low frequencies
    if not cutoffs:
        if axis == 'z':
            # z stiffness lower so cutoff needs to be lower
            cutoffs = [50, 9000]
        else:
            cutoffs = [130, 9000]
        cutoff_condition = (f > cutoffs[0]) & (f < cutoffs[1])
        p = p[cutoff_condition]
        f = f[cutoff_condition]

    # Set initial parameters for fitting
    fc = 1000 if axis != "z" else 200
    A = 1e3  # if axis != "z" else 500

    # if axis == "z":
    #     p, f = moving_average(p), moving_average(f)
    #     plt.plot(f, p, c='g')

    params = Parameters()
    params.add('Dv', value=A)
    params.add('fc', value=fc, min=-0, max=Fs/2)
    if anti_aliasing:
        params.add('fNyq', value=Fs/2)
        minner = Minimizer(anti_aliased_lorentzian, params, fcn_args=(f, p))
    else:
        minner = Minimizer(lorentzian, params, fcn_args=(f, p))
    result = minner.minimize()
    p_fit = p + result.residual

    Dv, fc = result.params['Dv'].value, result.params['fc'].value
    Dv_err, fc_err = result.params['Dv'].stderr, result.params['fc'].stderr
    b = np.sqrt(D/Dv)  # in m/V
    # Taking derivative to get error (df^2 = (df/dx dx)^2 + ...)
    if Dv_err:
        bErr = np.sqrt(D / (4 * Dv**3)) * Dv_err
    else:
        bErr = 0
    # k = 2πγfc
    k = 2*np.pi*gamma*fc*1e6  # in pN/μm
    if fc_err:
        kErr = 2*np.pi*gamma*fc_err*1e6  # in pN/μm
    else:
        kErr = 0

    b, bErr = b * 1e6, bErr * 1e6  # in μm/V
    h_correction = axl_correction if axis == 'z'else lat_correction
    k, kErr = k/h_correction, kErr/h_correction
    f = np.append([b], f)
    p_fit = np.append([k], p_fit)
    return np.array([f, p_fit])