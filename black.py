import pdb
import numpy as np

def black(nu_in, T):

  #h = 6.623e-34     ; Joule*s
  #k = 1.38e-23      ; Joule/K
  #c = 3e8           ; m/s
  # (2*h*nu_in^3/c^2)*(1/( exp(h*nu_in/k*T) - 1 )) * 10^29

  a0 = 1.4718e-21   # 2*h*10^29/c^2
  a1 = 4.7993e-11   # h/k

  #nu_cutoff = 1.2 * T / a1 * ( 3.0 * np.log(T / a1) - 46 )

  #n_nu_in = len(nu_in)
  #n_T = len(T)
  ##if (len(nu_in) > len(T)):
  ##if (np.shape(nu_in)[0] > np.shape(T)[0]):
  #if (n_nu_in > n_T):
  #  nonzero = np.where(nu_in <= nu_cutoff)
  #  ret = nu_in * 0.0
  #else: 
  #  nonzero = np.where(nu_cutoff >= nu_in)
  #  ret = nu_in * T

  #nz = nonzero[0]
  #if (nz[0] != -1): 
  #  ret[nonzero] = a0 * nu_in[nonzero]**3.0 / ( np.exp(a1 * nu_in[nonzero] / T) - 1.0 )

  #pdb.set_trace()

  #ret = a0 * nu_in**3.0 / ( np.exp(a1 * nu_in / T) - 1.0 )
  num = a0 * nu_in**3.0
  den = np.exp(a1 * np.outer(1.0/T,nu_in)) - 1.0 
  ret = num / den                 

  
  return ret
  #if (len(ret) == 1): 
  #  return ret[0] 
  #else: 
  #  return ret
