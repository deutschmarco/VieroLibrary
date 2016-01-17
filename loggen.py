import numpy as np

def loggen(minval, maxval, npoints, linear = None):
	points = np.arange(npoints)/(npoints - 1)
	if (linear != None):
		return (maxval - minval)*points + minval
	else:
		return 10.0 ** ( (np.log10(maxval/minval)) * points + np.log10(minval) )
