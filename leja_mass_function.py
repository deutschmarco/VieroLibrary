import pdb
import numpy as np
#import dschechter

def dschechter(X,P):
  #Fits a double Schechter function but using the same M*
  # X is alog10(M)
  # P[0]=alpha, P[1]=M*, P[2]=phi*, P[3]=alpha_2, P[4]=M*_2, P[5]=phi*_2

  #pdb.set_trace()
  rsch1 = np.log(10.) * P[2] * (10.**((X-P[1])*(1+P[0]))) * np.exp(-10.**(X-P[1]))
  rsch2 = np.log(10.) * P[5] * (10.**((X-P[4])*(1+P[3]))) * np.exp(-10.**(X-P[4]))

  return rsch1+rsch2


def leja_mass_function(mass,redshift,aqs=0):
  #aqs = 0  -  All
  #aqs = 1  -  Quiescent 
  #aqs = 2  -  Star Forming

  a1=np.array([-0.39,-0.10,-0.97])
  a2=np.array([-1.53,-1.69,-1.58])

  p1a=np.array([-2.46,-2.51,-2.88])
  p1b=np.array([ 0.07,-0.33, 0.11])
  p1c=np.array([-0.28,-0.07,-0.31])
  p2a=np.array([-3.11,-3.54,-3.48])
  p2b=np.array([-0.18,-2.31, 0.07])
  p2c=np.array([-0.03, 0.73,-0.11])
  ma= np.array([10.72,10.70,10.67])
  mb= np.array([-0.13, 0.00,-0.02])
  mc= np.array([ 0.11, 0.00, 0.10])
  
  mf=np.zeros([len(redshift),len(mass)])

  for i in range(len(redshift)):
    aone=a1[aqs]
    atwo=a2[aqs]
    phione=p1a[aqs]+p1b[aqs]*redshift[i]+p1c[aqs]*redshift[i]**2
    phitwo=p2a[aqs]+p2b[aqs]*redshift[i]+p2c[aqs]*redshift[i]**2
    mstar = ma[aqs]+ mb[aqs]*redshift[i]+ mc[aqs]*redshift[i]**2
    mone  = mstar
    mtwo  = mstar

    P0=np.array([aone,mone,10**phione,atwo,mtwo,10**phitwo])
    mf[i,:]=dschechter(mass,P0)

  return mf
