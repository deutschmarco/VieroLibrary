import numpy as np
import pdb
from VieroLibrary.invert_sed import sed
from VieroLibrary.invert_sed import sedint
from VieroLibrary.invert_sed import simple_flux_from_greybody
from VieroLibrary.black import black
from VieroLibrary.loggen import loggen
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from lmfit import Parameters, minimize, fit_report
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#pylab.figure()

L_sun = 3.839e26 # W
c = 299792458.0 # m/s

lambdavector = np.array([250,350,500])
nuvector = c * 1.e6 / lambdavector # Hz

nsed = 1e4
lambda_mod = loggen(1e3, 8.0, nsed) # microns
nu_mod = c * 1.e6/lambda_mod # Hz
cosmo = FlatLambdaCDM(H0 = 70.5 * u.km / u.s / u.Mpc, Om0 = 0.273)

nz = 5.#20.0#0
nl = 5.#20.0#0
nt = 5.#20.0#0
ngal=1

z_vector = loggen(0.1,5.,nz,linear=1)
l_vector = loggen(8.0,13.0,nl,linear=1)
t_vector = loggen(10.0,75.0,nt,linear=1)

z_l_t_matrix = np.array([z_vector,l_vector,t_vector])

betain = 2.0
alphain= 2.0

A_lookup_table = np.zeros([nl,nt,nz])
flux_250 = np.zeros([nl,nt,nz])
flux_350 = np.zeros([nl,nt,nz])
flux_500 = np.zeros([nl,nt,nz])

#nuvector = c * 1.e6 / lambdavector # Hz

for iz in np.arange(nz):
	conversion = 4.0 * np.pi *(1.0E-13 * cosmo.luminosity_distance(z_vector[iz]) * 3.08568025E22)**2.0 / L_sun # 4 * pi * D_L^2    units are L_sun/(Jy x Hz)
	for il in np.arange(nl):
		Lir = 10**l_vector[il] / conversion # Jy x Hz
		for it in np.arange(nt):
			Ain = 1.e-36
			fit_params = Parameters()
			fit_params.add('Ain', value= Ain)
			Pfin = minimize(sedint, fit_params, args=(nu_mod,Lir.value,ngal,t_vector[it]/(1.+z_vector[iz]),betain,alphain))
			flux_mJy=sed(Pfin.params,nuvector,ngal,t_vector[it]/(1.+z_vector[iz]),betain,alphain)
			flux_250[il,it,iz] = 1e3*flux_mJy[0][0]
			flux_350[il,it,iz] = 1e3*flux_mJy[0][1]
			flux_500[il,it,iz] = 1e3*flux_mJy[0][2]
			A_lookup_table[il,it,iz] = Pfin.params.valuesdict()['Ain'] 
			#print str(z_vector[iz]) ,str(l_vector[il]) ,str(t_vector[it]) ,str(A_lookup_table[il,it,iz]), str(1e3*flux_mJy) 

