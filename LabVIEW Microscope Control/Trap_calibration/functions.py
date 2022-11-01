from lmfit import Minimizer, Parameters
from numpy import pi


def anti_aliased_lorentzian(params, f, data, sides=1):
    N = 10
    N_values = [i-N for i in range(2*N+1)]
    model = 0
    for n in N_values:
        model += (params['Dv']/(sides*pi**2)) / \
            ((f + n*params['fNyq']/2)**2 + params['fc']**2)
    return model - data


def lorentzian(params, f, data, sides=1):
    model = (params['Dv']/(sides*pi**2))/(f**2 + params['fc']**2)
    return model - data
