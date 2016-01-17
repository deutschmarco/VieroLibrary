import numpy as np
import pdb
from my_first_functions import lambda_to_ghz
from number_galaxies_given_redshifts import *
from viero_flux import *
def viero_cib(znodes,wavelength,aqs=2,NM=20.0,Mlo=8.0,Mhi=13.0): 
  #redshift centers array
  NZ=len(znodes)-1
  z=np.zeros(NZ)
  for i in range(NZ):
    z[i]=(znodes[i]+znodes[i+1])/2.0
  #mass array
  M=np.linspace(Mlo,Mhi,num=NM)

  #Number of galaxies in sky per redshift slice
  #DONT FORGET NEEDS MODIFYING AT HIGH Z
  numbers=np.zeros([NM,NZ])
  for i in range(NZ):
    numbers[:,i]=number_galaxies_given_redshifts(znodes[i],znodes[i+1],aqs=2,NM=NM,Mlo=Mlo,Mhi=Mhi)

  fluxes=viero_flux(M,z,wavelength)

  nuInu=np.zeros([NM,NZ])
  #FIND CIB
  for i in range(NZ):
    #print 1e-3*fluxes[:,i]*(lambda_to_ghz(wavelength)*1e9)*1e-26*1e9*numbers[:,i]/(4.*3.141592654)
    nuInu[:,i]=1e-3*fluxes[:,i]*(lambda_to_ghz(wavelength)*1e9)*1e-26*1e9*numbers[:,i]/(4.*3.141592654)
  #pdb.set_trace()
      
  return nuInu
