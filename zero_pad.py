if __name__ == "__main__":
  import sys
  fib(int(sys.argv[1]))

import numpy as np
def zero_pad(cmap,l2=0):
  ms=np.shape(cmap)
  if ms[0] <= l2 and ms[1] <=l2:
    zmap=np.zeros([l2,l2])
    zmap[:ms[0],:ms[1]]=cmap 
  else:
    zmap=cmap
  return zmap
