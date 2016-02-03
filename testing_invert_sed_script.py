import numpy as np
import pdb
from VieroLibrary.invert_sed import sedint
from VieroLibrary.invert_sed import simple_flux_from_greybody
from VieroLibrary.invert_sed import single_simple_flux_from_greybody
from VieroLibrary.black import black
from VieroLibrary.loggen import loggen
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from lmfit import Parameters, minimize, fit_report

zin = np.array([1.0])
Lin = np.array([1e12])#,2e12,1e13])
lambdavector = np.array([250,350,500])
nwv = len(lambdavector)

L_sun = 3.839e26 # W
c = 299792458.0 # m/s
nuvector = c * 1.e6 / lambdavector # Hz

nsed = 1e4
lambda_mod = loggen(1e3, 8.0, nsed) # microns
nu_mod = c * 1.e6/lambda_mod # Hz
cosmo = FlatLambdaCDM(H0 = 70.5 * u.km / u.s / u.Mpc, Om0 = 0.273)
conversion = 4.0 * np.pi *(1.0E-13 * cosmo.luminosity_distance(zin) * 3.08568025E22)**2.0 / L_sun # 4 * pi * D_L^2    units are L_sun/(Jy x Hz)
int_sed = Lin / conversion # Jy x Hz


Ain = np.array([1.e-36])#,3.e-36,1.2e-36]) #good starting parameter
Tin = np.array([23.0])#,24.0,28.0])
betain = np.array([2.0])#,2.0,2.0])
alphain =np.array([2.0])#,2.0,2.0])
#Ain = np.array([1.e-36]) #good starting parameter
#Tin = np.array([23.0])
#betain = np.array([2.0])
#alphain =np.array([2.0])

test1 = black(nu_mod,Tin)

#pdb.set_trace()


#Ain = (2.e-36,3.e-36,1.2e-36)
#Ain = (2.e-3,3.e-3,1.2e-3)

ngal = len(Ain)
print ngal
fit_params = Parameters()
fit_params.add('Ain', value= Ain[0])
#pdb.set_trace()
#fit_params.add('Tin', value= Tin, vary = False)
#fit_params.add('betain', value= betain, vary = False)
#fit_params.add('alphain', value= alphain, vary = False)
#test = sedint(fit_params, nu_mod, int_sed.value, ngal)
#test = sedint(fit_params, nu_mod, int_sed.value, ngal, Tin, betain, alphain)

#pdb.set_trace()

test2 = single_simple_flux_from_greybody(lambdavector, Trf = Tin, b = betain, Lrf = Lin, zin = zin)
print 1e3*test2