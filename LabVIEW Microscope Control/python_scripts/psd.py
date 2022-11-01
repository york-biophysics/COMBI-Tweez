import numpy as np
def buildArray(f0, df, y):
	x = np.arange(0, len(y))*df + f0
	return x