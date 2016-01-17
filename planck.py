import numpy as np

h = 6.626e-34
c = 3.0e+8
k = 1.38e-23

def planck(wav, T): #where wav is in m?
	#nuvector = c * 1.e6 / lambdavector # Hz from microns??
    a = 2.0 * h * c**2
    b = h * c / (wav * k * T)
    intensity = a / ( (wav**5) * (np.exp(b) - 1.0) )
    return intensity