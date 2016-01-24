import pdb
import numpy as np
from VieroLibrary.loggen import loggen
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
pi=3.141592653589793

def comoving_distance(z,h=0.6774,OmM=0.3089,OmL=0.6911,Omk=0.0,dz=0.001,inverse_h=None):
    H0 = 100. * h #km / s / Mpc
    D_hubble = 3000. / h # h^{-1} Mpc = 9.26e25 / h; (meters)
    cosmo = FlatLambdaCDM(H0 = H0 * u.km / u.s / u.Mpc, Om0 = OmM)
    n_z = z/dz 
    i_z = np.arange(n_z)*dz

    D_c = 0.0
    for i in i_z:
        E = np.sqrt(OmM*(1.+i)**3. + OmL)
        D_c += D_hubble * dz / E

    return D_c

def comoving_volume(zed1, zed2, mpc=None):
    if zed1 < zed2:
        z1 = zed1
        z2 = zed2
    else:
        z1 = zed2
        z2 = zed1
    comovo=1e-9* 4./3.* pi * (comoving_distance(z2)**3. - comoving_distance(z1)**3.)
    if mpc != None:
        comovo *= 1e3**3.0

    return comovo

def comoving_number_density(number, area, z1, z2, ff=1.0, mpc=None, verbose=None):
    #if z2 != None: zin2 = 0.0
    vol = comoving_volume(z2,z1,mpc=1)
    num = (number/(area*ff)) * (180.0/pi)**2.0 * 4.0 * pi
    comovnumden=num/vol

    return comovnumden
