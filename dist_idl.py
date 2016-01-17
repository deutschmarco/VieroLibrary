
def dist_idl(n1,m1=None):
    ''' Copy of IDL's dist.pro
    Create a rectangular array in which each element is 
    proportinal to its frequency'''

    import numpy as np

    #n1 = n[0]
    if m1 == None:
        m1 = n1
    #print m1
    #print n1

    x = np.arange(float(n1))
    for i in range(len(x)): x[i]= min(x[i],(n1 - x[i])) ** 2.

    a = np.zeros([float(n1),float(m1)])

    i2 = m1/2 + 1

    #print i2
    for i in np.arange(i2):
        y = np.sqrt(x + i ** 2.)
        a[:,i] = y
        if i != 0:
            a[:,m1-i]=y 

    return a


