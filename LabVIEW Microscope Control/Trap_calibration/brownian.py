import matplotlib.pyplot as plt
import numpy as np
import statistics
from scipy.optimize import curve_fit
from settings import kT


def k_brown(sig, b, axis, filename, plot=False):
    """
    Take signal and get brownian result
    b is sensitivity in μm/V
    """
    # Calibrated position + stiffness
    x = sig*b*1e3  # bead position [nm]
    x = x - statistics.mean(x)
    k_brown = 4.114 / np.std(x)**2
    k_brown = k_brown * 1e3  # convert to pN/um

    return k_brown

# if plot:
#     plt.clf()
#     hist, bins = np.histogram(x, bins=500, density=True)
#     centres = (bins[:-1] + bins[1:]) / 2
#     popt, pcov = curve_fit(gaussian, centres, hist, p0=[0.0, 0.0, np.std(x)])
#     plt.scatter(centres, hist, marker='x')
#     plt.plot(centres, gaussian(centres, *popt), 'r--')
#     plt.xlim(xmin=-3*popt[2], xmax=3*popt[2])
#     save_file = 'plots/hist_' + axis + '_' + filename[0:-4]
#     plt.savefig(save_file, bbox_inches='tight')
#     k_brown = (1e6*kT) / (popt[2]**2)
