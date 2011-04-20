#!/usr/bin/env python

import subprocess as spc
import compbio.utils.bsub as bsub
import compbio.utils.bs_macros as bsm
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

def bic_clustering(input_dict, run_id):
    return bsm.runmat('ap_max_bic', input_dict, run_id)

def test_bsubfun(input_dict, run_id):
    '''
A sample function to demonstrate the calling of a matlab script (here, 
ap_frompy) from within python. Taking an input dictionary and a run_id,
this script is designed to be called using the 'eyeball' class from 
utils/bsub.py.

inputs:
  input_dict: {similarities: a similarity matrix for the input points,
               self_similarity: a single value for the self similarity
                                of datapoints. Control cluster size.

outputs:
  outpt_dict: {indexes: cluster exemplar indices.}

'''
    return bsm.runmat('ap_frompy', input_dict, run_id)

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
    
