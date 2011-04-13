#!/usr/bin/env python
'''
Check mail and process scanned pdfs into full color and blue channel pngs.
'''

#STANDARD IMPORTS
import os, subprocess, itertools as it, numpy as np, Image, re
from matplotlib import image as mpimg
from numpy import *

import scipy.signal as ss
#LOCAL IMPORTS
import attachments
import simplejson

att_dir = os.path.join(os.environ.get("HOME"),'.anNote')
def run():

  for r,d,fs in os.walk(os.path.join(att_dir,'tmp')):
    null = [ os.remove(os.path.join(r,f)) for f in fs]
  attachments.check_att(att_dir)
  pdfs = dict([(os.path.join(att_dir,'tmp/tmp_{0:04d}.pdf'.format(idx)), \
                  os.path.join(att_dir,f)) 
               for idx, f in enumerate(os.listdir(att_dir)) if '.pdf' in f.lower()])
  for k,v in pdfs.iteritems(): os.rename(v,k)

  #Convert PDFs
  res = 300
  cvsub = ['''convert -density {2} {0} {1}; rm {0}'''.\
             format(p , p.replace('.pdf','.png'), res) 
           for p in pdfs.keys()]
  for c in cvsub: 
    print 'calling for ' +  c
    subprocess.call(c,shell = True)

  #Find all files produced by convert
  inp_files=[ os.path.join(att_dir,'tmp',e) for e in 
              it.chain(\
      *[filter( lambda x: os.path.splitext(os.path.basename(key))[0]\
                  in x and True or False,
                os.listdir(os.path.join(att_dir,'tmp')))
        for key in  pdfs.keys()])]
  
  bluechannels = []
  x_inches = .25
  y_inches = .15
  skip = 5

  #open them and get the blue channels
  for i in inp_files:
    full = mpimg.imread(i)
    
    if isflipped(full):
      full = transpose(full,(1,0,2))

    blue= squeeze(full[:,:,2])
    others = np.sum(full[:,:,0:2],2)
    blue = blue - others
    blue[less(blue,0)] = 0.
    blue[greater(blue,0)] = 1.
    xrad = floor(res/skip*x_inches)
    yrad = floor(res/skip*y_inches)


    bluesmall = blue[::skip,::skip]
    bluesmall = ss.order_filter( bluesmall,  ones((xrad,1)), xrad-1)
    bluesmall = ss.order_filter( bluesmall,  ones((1,yrad)), yrad-1)

    #FLIP THE CLUSTERS JUST TO BE CONFUSING
    cls = cluster_img(bluesmall)
    cls = [cl.T * 5 for cl in cls]
      
    root = os.path.join(att_dir,'pages')
    if not os.path.isdir(root): os.mkdir(root)
    num_max = max(array([0] + list(it.chain(*[re.findall(re.compile('[\d]+'),f) for f in os.listdir(root)])),int))
    this_folder = os.path.join(root,'page_{0:05d}'.format(num_max+1))
    os.mkdir(this_folder)
    
    Image.fromarray(transpose(array(blue*255,dtype=np.uint8),(1,0))).save(open(os.path.join(this_folder,'blue.png'),'w'))
    Image.fromarray(transpose(array(full*255,dtype=np.uint8),(1,0,2))).save(open(os.path.join(this_folder,'full.png'),'w'))
    
    outline_folder = os.path.join(this_folder,'highlights')
    for idx,c in enumerate(cls):
      bounds =array([ np.min(c[0]),np.min(c[1]),np.max(c[0]),np.max(c[1])])
      b0 = array(bounds)
      bounds =bounds + 100 *array( [-1,-1,1,1] )
      clip_bounds(bounds, shape(blue.T))
      subimg = full[bounds[0]:bounds[2], bounds[1]:bounds[3]]
      coords = { 'full_bounds':list(bounds),
               'cluster_bounds':list(b0)}
      Image.fromarray(transpose(array(subimg*255,dtype=np.uint8),(1,0,2))).save(open(os.path.join(this_folder,'hl_{0:02d}.png').format(idx),'w'))
      fopen = open('hl_{0:02d}.txt'.format(idx), 'w')
      fopen.write(simplejson.dumps(coords))
      fopen.close()
  

def clip_bounds(rect,dims):
  rect[less(rect,0)] = 0.
  rect = np.min(vstack((rect, [dims[0],dims[1],dims[0],dims[1]])),0)


def isflipped(img):
  '''check whether a (THRESHOLDED) image is flipped. In other words, do the lines run from x=0... or y=0.... Note that an image that I define as "not flipped" for these purposes will by default appear flipped in matplotlib. Note also, that if the image is full of pictures, there is no guarantee that this method will work. What is for sure is that if the image is flipped improperly, later clusterings will most likely fail to identify blocks of text.

input:
  A thresholded image that will be shrunk and checked for line order.

output:
  True if words run in the y-direction.

'''
  imsml = greater(1- np.min(img[::5,::5],2),.5)
  dims = shape(imsml)
  smx = ss.order_filter(imsml, ones((dims[0]/72*2 + 1,1)), dims[0]/36 -1)
  smy = ss.order_filter(imsml, ones((1,dims[1]/72*2 + 1 )), dims[1]/36 -1)
  clx = cluster_img(smx)
  cly = cluster_img(smy)

  if len(clx)> len(cly): return True
  else: return False

def cluster_img(img):
  '''Cluster the points in an input image into disconnected subsets

input:
  img: a black and white [nx, nm] image.

output:
  clusters: a [2,nc] list of cluster elements.
  '''
  clusters = []

  nzs = nonzero(img)
  sx,sy = shape(img)

  
  nzs = set([n[0] + n[1]*sx for n in zip(*nzs)])
  neighbors = lambda pt: set([pt + 1, pt -1, pt + sx , pt - sx])
  
  while nzs:
    queued = set([nzs.pop()])
    added = set()
    while queued:
      p = queued.pop()
      nzs.difference_update([p])
      added.update([p])
      queued.update([n for n in neighbors(p) 
                     if ( not (n in added or n in queued)) 
                     and n in nzs])

    clusters.append(array( list(added) ))
  return [array([floor(c/sx), floor(mod(c, sx))]) for c in clusters]


if __name__ == '__main__':
  run()
