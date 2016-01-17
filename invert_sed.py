import pdb
import numpy as np
from VieroLibrary.black import black
from VieroLibrary.loggen import loggen
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from lmfit import Parameters, minimize, fit_report

L_sun = 3.839e26 # W
c = 299792458.0 # m/s

def find_nearest_index(array_in,value):
	ng = len(value)
	#idx = (np.abs(array_in-value)).argmin()
	idx = (np.abs(array_in-np.reshape(value,(ng,1)))).argmin(axis=1)
	return idx

def sed(p, nu_in, ng, T, betain, alphain): #; m = [A, T, Beta, Alpha] - return SED (array) in Jy
	#P[0] = A
	#P[1] = T
	#P[2] = Beta
	#P[3] = Alpha

	v = p.valuesdict()
	A = v['Ain']
	#T = v['Tin']
	#betain = v['betain']
	#alphain = v['alphain']

	ns = len(nu_in)
	base = 2.0 * (6.626)**(-2.0 - betain - alphain) * (1.38)**(3. + betain + alphain) / (2.99792458)**2.0
	expo = 34.0 * (2.0 + betain + alphain) - 23.0 * (3.0 + betain + alphain) - 16.0 + 26.0
	K = base * 10.0**expo
	w_num = A * K * (T * (3.0 + betain + alphain))**(3.0 + betain + alphain) 
	w_den = (np.exp(3.0 + betain + alphain) - 1.0)
	w_div = w_num/w_den 
	nu_cut = (3.0 + betain + alphain) * 0.208367e11 * T

	graybody = np.reshape(A,(ng,1)) * nu_in**np.reshape(betain,(ng,1)) * black(nu_in, T) / 1000.0 
	powerlaw = np.reshape(w_div,(ng,1)) * nu_in**np.reshape(-1.0 * betain,(ng,1))
	graybody[np.where(nu_in >= np.reshape(nu_cut,(ng,1)))]=powerlaw[np.where(nu_in >= np.reshape(nu_cut,(ng,1)))]

	return graybody

def sed_old(p, nu_in): #; m = [A, T, Beta, Alpha] - return SED (array) in Jy
	#P[0] = A
	#P[1] = T
	#P[2] = Beta
	#P[3] = Alpha

	v = p.valuesdict()
	A = v['Ain']
	T = v['Tin']
	betain = v['betain']
	alphain = v['alphain']

	ns = len(nu_in)
	base = 2.0 * (6.626)**(-2.0 - betain - alphain) * (1.38)**(3. + betain + alphain) / (2.99792458)**2.0
	expo = 34.0 * (2.0 + betain + alphain) - 23.0 * (3.0 + betain + alphain) - 16.0 + 26.0
	K = base * 10.0**expo
	w = A * K * (T * (3.0 + betain + alphain))**(3.0 + betain + alphain) / (np.exp(3.0 + betain + alphain) - 1.0)
	nu_cut = (3.0 + betain + alphain) * 0.208367e11 * T
	
	graybody = np.zeros(ns)

	if (len(nu_in) > 1):
		nu_gt_cut = np.where(nu_in > nu_cut)
		n_zero1 = len(nu_gt_cut)
		nu_le_cut = np.where(nu_in <= nu_cut)
		n_zero2 = len(nu_le_cut)

		if (n_zero1 > 0): graybody[nu_gt_cut] = w * nu_in[nu_gt_cut]**(-1.0 * betain)
		if (n_zero2 > 0): graybody[nu_le_cut] = A * nu_in[nu_le_cut]**betain * black(nu_in[nu_le_cut], T) / 1000.0

	else:
		for i in range(len(nu_in)): 
			if (nu_in[i] > nu_cut):
				graybody[i] = w * nu_in[i]**(-1.0 * betain) 
			else: 
				graybody[i] = A * nu_in[i]**betain * black(nu[i], T)/1000.
	return graybody

def sedint(p, nu_in, Lir, ng, T, betain, alphain): # m = [A, T, Beta, Alpha] - return integrated SED flux (one number) in Jy x Hz
	#P[0] = A
	#P[1] = T
	#P[2] = Beta
	#P[3] = Alpha

	v = p.valuesdict()
	A0 = v['Ain']
	A=np.asarray(A0)
	#T = v['Tin']
	#betain = v['betain']
	#alphain = v['alphain']


	#print 'A is ' + str(A)
	ns = len(nu_in)
	#pdb.set_trace()
	#ng = len(A)
	base = 2.0 * (6.626)**(-2.0 - betain - alphain) * (1.38)**(3. + betain + alphain) / (2.99792458)**2.0
	expo = 34.0 * (2.0 + betain + alphain) - 23.0 * (3.0 + betain + alphain) - 16.0 + 26.0
	K = base * 10.0**expo
	w_num = A * K * (T * (3.0 + betain + alphain))**(3.0 + betain + alphain) 
	w_den = (np.exp(3.0 + betain + alphain) - 1.0)
	w_div = w_num/w_den 
	nu_cut = (3.0 + betain + alphain) * 0.208367e11 * T
	
	#nu_cut_ind = find_nearest_index(nu_in,nu_cut)

	graybody = np.reshape(A,(ng,1)) * nu_in**np.reshape(betain,(ng,1)) * black(nu_in, T) / 1000.0 
	powerlaw = np.reshape(w_div,(ng,1)) * nu_in**np.reshape(-1.0 * betain,(ng,1))
	graybody[np.where(nu_in >= np.reshape(nu_cut,(ng,1)))]=powerlaw[np.where(nu_in >= np.reshape(nu_cut,(ng,1)))]


	dnu = nu_in[1:ns] - nu_in[0:ns-1]
	dnu = np.append(dnu[0],dnu)

	return np.ravel([np.sum(graybody * dnu, axis=1) - Lir]) 

def sedint2(p, nu_in, Lir, ng): # m = [A, T, Beta, Alpha] - return integrated SED flux (one number) in Jy x Hz
	#P[0] = A
	#P[1] = T
	#P[2] = Beta
	#P[3] = Alpha

	v = p.valuesdict()
	A = v['Ain']
	T = v['Tin']
	betain = v['betain']
	alphain = v['alphain']


	print 'A is ' + str(A)
	ns = len(nu_in)
	#ng = len(A)
	base = 2.0 * (6.626)**(-2.0 - betain - alphain) * (1.38)**(3. + betain + alphain) / (2.99792458)**2.0
	expo = 34.0 * (2.0 + betain + alphain) - 23.0 * (3.0 + betain + alphain) - 16.0 + 26.0
	K = base * 10.0**expo
	w_num = A * K * (T * (3.0 + betain + alphain))**(3.0 + betain + alphain) 
	w_den = (np.exp(3.0 + betain + alphain) - 1.0)
	w_div = w_num/w_den 
	nu_cut = (3.0 + betain + alphain) * 0.208367e11 * T
	
	#nu_cut_ind = find_nearest_index(nu_in,nu_cut)

	graybody = np.reshape(A,(ng,1)) * nu_in**np.reshape(betain,(ng,1)) * black(nu_in, T) / 1000.0 
	powerlaw = np.reshape(w_div,(ng,1)) * nu_in**np.reshape(-1.0 * betain,(ng,1))
	graybody[np.where(nu_in >= np.reshape(nu_cut,(ng,1)))]=powerlaw[np.where(nu_in >= np.reshape(nu_cut,(ng,1)))]

	#pdb.set_trace()

	dnu = nu_in[1:ns] - nu_in[0:ns-1]
	dnu = np.append(dnu[0],dnu)

	return np.ravel([np.sum(graybody * dnu, axis=1) - Lir]) 

def sedint_old(p, nu_in, Lir): # m = [A, T, Beta, Alpha] - return integrated SED flux (one number) in Jy x Hz
	#P[0] = A
	#P[1] = T
	#P[2] = Beta
	#P[3] = Alpha

	v = p.valuesdict()
	A = v['Ain']
	T = v['Tin']
	betain = v['betain']
	alphain = v['alphain']


	ns = len(nu_in)
	ng = len(A)
	base = 2.0 * (6.626)**(-2.0 - betain - alphain) * (1.38)**(3. + betain + alphain) / (2.99792458)**2.0
	expo = 34.0 * (2.0 + betain + alphain) - 23.0 * (3.0 + betain + alphain) - 16.0 + 26.0
	K = base * 10.0**expo
	w_num = A * K * (T * (3.0 + betain + alphain))**(3.0 + betain + alphain) 
	w_den = (np.exp(3.0 + betain + alphain) - 1.0)
	w_div = w_num/w_den #A * K * (T * (3.0 + betain + alphain))**(3.0 + betain + alphain) / (np.exp(3.0 + betain + alphain) - 1.0)
	nu_cut = (3.0 + betain + alphain) * 0.208367e11 * T
	
	#Original Function
	graybody = np.zeros(ns)
	nu_gt_cut = np.where(nu_in > nu_cut)
	nu_le_cut = np.where(nu_in <= nu_cut)
	n_zero1 = len(nu_gt_cut)
	n_zero2 = len(nu_le_cut)
	if (n_zero1 > 0): graybody[nu_gt_cut] = w_div * nu_in[nu_gt_cut]**(-1.0 * betain)
	if (n_zero2 > 0): graybody[nu_le_cut] = A * nu_in[nu_le_cut]**betain * black(nu_in[nu_le_cut], T) / 1000.0

	#graybody = np.ravel(np.zeros([ng,ns]))

	#nu_gt_cut = (np.where(nu_in > np.reshape(nu_cut, (ng,1)) ) ) [1]
	#nu_gt_cut = np.reshape(nu_gt_cut,(ng,ns))
	##nu_gt_cut = np.where(np.tile(nu_in, (ng,1)) > np.reshape(nu_cut, (ng,1)) )
	#n_zero1 = len(nu_gt_cut)
	#nu_le_cut = (np.where(nu_in <= np.reshape(nu_cut, (ng,1)) ) )[1] 
	#nu_le_cut = np.reshape(nu_le_cut,(ng,ns))
	##nu_le_cut = np.where(np.tile(nu_in, (ng,1)) <= np.reshape(nu_cut, (ng,1)) )
	#n_zero2 = len(nu_le_cut)
	#pdb.set_trace()

	##if (n_zero1 > 0): graybody[nu_gt_cut] = np.reshape(w_div, (ng,1)) * np.tile(nu_in, (ng,1))[nu_gt_cut]**(-1.0 * np.reshape(betain,(ng,1)) )
	##if (n_zero2 > 0): graybody[nu_le_cut] = np.reshape(A,(ng,1))      * np.tile(nu_in, (ng,1))[nu_le_cut]**np.reshape(betain,(ng,1)) * black(np.tile(nu_in, (ng,1))[nu_le_cut], T) / 1000.0

	#blackbody = black(nu_in,np.reshape(T,(ng,1)))
	#if (n_zero1 > 0): 
	#	graybody[nu_gt_cut] = np.reshape(w_div, (ng,1)) * np.tile(nu_in, ng)[nu_gt_cut]**(-1.0 * np.reshape(betain,(ng,1)) )
	#if (n_zero2 > 0): 
	#	graybody[nu_le_cut] = np.reshape(A,(ng,1))      * np.tile(nu_in, ng)[nu_le_cut]**np.reshape(betain,(ng,1)) * black(np.tile(nu_in, (ng,1))[nu_le_cut], T) / 1000.0

	dnu = nu_in[1:ns] - nu_in[0:ns-1]
	dnu = np.append(dnu[0],dnu)

	return np.ravel([np.sum(graybody * dnu) - Lir]) 

def simple_flux_from_greybody(lambdavector, Trf = None, b = None, Lrf = None, zin = None, ngal = None):
	''' 
	Return flux densities at any wavelength of interest (in the range 1-10000 micron),
	assuming a galaxy (at given redshift) graybody spectral energy distribution (SED),
	with a power law replacing the Wien part of the spectrum to account for the
	variability of dust temperatures within the galaxy. The two different functional
	forms are stitched together by imposing that the two functions and their first
	derivatives coincide. The code contains the nitty-gritty details explicitly.
	Cosmology assumed: H0=70.5, Omega_M=0.274, Omega_L=0.726 (Hinshaw et al. 2009)

	Inputs:
	alphain = spectral index of the power law replacing the Wien part of the spectrum, to account for the variability of dust temperatures within a galaxy [default = 2; see Blain 1999 and Blain et al. 2003]
	betain = spectral index of the emissivity law for the graybody [default = 2; see Hildebrand 1985]
	Trf = rest-frame temperature [in K; default = 20K]
	Lrf = rest-frame FIR bolometric luminosity [in L_sun; default = 10^10]
	zin = galaxy redshift [default = 0.001]
	lambdavector = array of wavelengths of interest [in microns; default = (24, 70, 160, 250, 350, 500)];
	
	AUTHOR:
	Lorenzo Moncelsi [moncelsi@caltech.edu]
	
	HISTORY:
	20June2012: created in IDL
	November2015: converted to Python
	'''

	nwv = len(lambdavector)
	nuvector = c * 1.e6 / lambdavector # Hz

	nsed = 1e4
	lambda_mod = loggen(1e3, 8.0, nsed) # microns
	nu_mod = c * 1.e6/lambda_mod # Hz

	cosmo = FlatLambdaCDM(H0 = 70.5 * u.km / u.s / u.Mpc, Om0 = 0.273)
	conversion = 4.0 * np.pi *(1.0E-13 * cosmo.luminosity_distance(zin) * 3.08568025E22)**2.0 / L_sun # 4 * pi * D_L^2    units are L_sun/(Jy x Hz)

	Lir = Lrf / conversion # Jy x Hz

	Ain = np.zeros(ngal) + 1.0e-36 #good starting parameter
	betain =  np.zeros(ngal) + b 
	alphain=  np.zeros(ngal) + 2.0

	fit_params = Parameters()
	fit_params.add('Ain', value= Ain)
	#fit_params.add('Tin', value= Trf/(1.+zin), vary = False)
	#fit_params.add('betain', value= b, vary = False)
	#fit_params.add('alphain', value= alphain, vary = False)

	#pdb.set_trace()
	#THE LM FIT IS HERE
	#Pfin = minimize(sedint, fit_params, args=(nu_mod,Lir.value,ngal))
	Pfin = minimize(sedint, fit_params, args=(nu_mod,Lir.value,ngal,Trf/(1.+zin),b,alphain))

	#pdb.set_trace()
	flux_mJy=sed(Pfin.params,nuvector,ngal,Trf/(1.+zin),b,alphain)

	return flux_mJy






