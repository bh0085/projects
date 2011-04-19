#!/usr/bin/env python

import subprocess as spc
import compbio.utils.bsub as bsub
import inspect
import os, sys, inspect
import scipy.io as sio
from numpy import *

def make_tests():
    mirnaf = os.path.join(os.path.dirname(inspect.stack()[0][1]), 'miRNA.mat')
    mirna = sio.loadmat(mirnaf)
    expr = mirna['expression']

    e_norms = sum(expr**2,1)

    cluster_dists = e_norms[:,newaxis] + e_norms[newaxis,:] \
        - 2 * dot(expr, expr.T)
    sims = - cluster_dists


    inp_dicts = []
    percentiles = [.01,1.,5.,10.,20.,25.,50.,75.]
    for p in percentiles:
        inp_dicts.append(dict(similarities = sims,
                              self_similarity = percentile(sims.flatten(),p)))
        
    eyeball = bsub.eyeball(os.path.abspath(inspect.stack()[0][1]), 
                           ['test_bsubfun'], 
                           inp_dicts)
    return eyeball
def test_bsubfun(input_dict, run_id):
    #sub = spc.Popen('find', shell =True, stdout = spc.PIPE).\
    #   communicate()[0]
    #out_dict = dict(output = sub)
    
    

    tmpnames = bsub.tmp_fnames(run_id,2)
    sio.savemat(tmpnames[0], inp_dict)
    
    cstr = 'mlab -r ap_frompy({0}, {1})'.\
                        format(tmpnames[0],tmpnames[1])


    #raise Exception()
    sub = spc.Popen(cstr,shell = True, stdout = spc.PIPE).\
        communicate()[0]
    out_dict = sio.loadmat(tmpnames[1])['out_struct']
    return out_dict

def usage():
  print '''
usage: btol.py run_id
1
Run a batch process on BTOL with inputs stored in 
data/batch/inputs/{run_id}.inp in pickle serial.
'''
  exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 3: usage()
    run_id = sys.argv[2]
    run_func = globals()[sys.argv[1]]
    input_dict = bsub.load_inp(run_id)
    output_dict = run_func(input_dict, run_id)
    bsub.save_out( output_dict, run_id)
    
