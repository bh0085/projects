import scipy.signal as ss
import numpy as np
from numpy import *
def dist(dims, c = False):
  if not shape(dims): dims = (dims, dims)
  xy = reshape(vstack([mod(arange(dims[0]*dims[0]),dims[1]) 
                       ,arange(dims[0]*dims[1]) / dims[0]]).T
                      ,(dims[0],dims[1],2)) 
  
  xy =  sqrt(sum(xy**2,2))
  xy =np.min(dstack([xy,xy[::-1]]),2)
  xy =np.min(dstack([xy.T,xy.T[::-1]]),2)
  if c:
    xy = roll(roll(xy, dims[0]/2, 0), dims[1]/2,1)
    
  return xy
