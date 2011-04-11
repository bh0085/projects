#!/usr/bin/env python
'''
Check mail and process scanned pdfs into full color and blue channel pngs.
'''

#STANDARD IMPORTS
import os, subprocess, itertools as it, numpy as np
from matplotlib import image as mpimg
from numpy import *
#LOCAL IMPORTS
import attachments

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
  #open them and get the blue channels
  for i in inp_files:
    img = mpimg.imread(i)
    blue= squeeze(img[:,:,2])
    others = np.sum(img[:,:,0:2],2)
    final = blue - others
    final[less(final,0)] = 0.
    bluechannels.append(final)

  return bluechannels
  

if __name__ == '__main__':
  run()
