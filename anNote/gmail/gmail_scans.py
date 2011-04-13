import attachments

import subprocess, os
import matplotlib.image as mpimg
import numpy as np
from scipy import signal as ss
from numpy import *
from compbio.utils.idl import dist

def get_all():
  attachments.check_att(att_dir )
  

def process(my_smooth = True, order_filt = True):
  img_res = 100
  pdf_files = [f for f in os.listdir(att_dir) if '.pdf' or '.PDF' in f]
  for f in pdf_files:
    inp_f = os.path.join(att_dir,f)
    out_f = os.path.join(att_dir,os.path.splitext(f)[0] + '.png')
    cmd = 'convert -density {2} {0} {1}'.format(inp_f, out_f, img_res)
    out = subprocess.Popen(cmd, shell = True)
    out.communicate()
    #os.remove(inp_f)
  
  pngs = [os.path.join(att_dir,f) for f in os.listdir(att_dir) if '.png' in f]
  for p in pngs:
    img = mpimg.imread(p)
    blue= squeeze(img[:,:,2])
    others = np.sum(img[:,:,0:2],2)
    final = blue - others
    
    final[less(final,0)] = 0.

    winsize = 5
    wid = winsize/2
    funny = False
    circle_only = True
    
    if order_filt:
      order_mask = less(dist(winsize, c = True), wid)
      n = len(nonzero(order_mask)[0])
      lowest =int( n * .3)
      convo = ss.order_filter(final, order_mask, lowest)
      convo = ss.medfilt(convo,7)
      return final, convo, order_mask
    if circle_only:
      gaussian = 1.* less(dist(winsize,c = True),wid)
    elif funny:
      gaussian = exp(-1 * (dist(winsize) /wid )**2)
    else:
      gaussian = exp(-1 * (dist(20, c = True)/20 /wid )**2)
    gaussian /= sum(gaussian)
    
    
    if my_smooth:
      convo = ss.convolve2d(final, gaussian)
    else:
      convo =  blur_image(final, 20)
    return final, convo, gaussian


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = ss.convolve(im,g, mode='valid')
    return(improc)
