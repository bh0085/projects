#!/usr/bin/env python

import subprocess as spc
import compbio.utils.bsub as bsub
import compbio.config as config
import compbio.utils.bs_macros as bsm

import os, sys, inspect, pipes
import scipy.io as sio
from numpy import *


#A few batch processes.
def bic_clustering(input_dict, run_id):
    '''
A matlab/bsub process to compute the BIC maximal clustering for an
input dictionary containing a similarity matrix.

inputs:
  input_dict:  {similarities: a similarity matrix}

outputs:
  output_dict: {inds:cluster exemplar indices,      (MAX BIC)
                self_similarity:float, self similarity (MAX BIC)
                
                inds_[#]: (same as above, ALL BIC)
                self_similarity_[#}: (...)
                bic_[#]: (...)
                }
'''
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

#Local launchpoint for the batch scripts
def local_launch():
    bs_cmd = pipes.quote( bsub.cmd(os.path.relpath(inspect.stack()[0][1],config.root),
                       'remote_make_tests',
                       run_id =bsub.get_run_id(0)))
    
    cmd = '''
ssh tin '{0}'
'''.format(bs_cmd)
    return cmd
  
    
#Remote launchpoint for bsub.
def remote_make_tests():
    mirnaf = os.path.join(os.path.dirname(inspect.stack()[0][1]), 'miRNA.mat')
    mirna = sio.loadmat(mirnaf)
    expr = mirna['expression']
    e_norms = sum(expr**2,1)
    cluster_dists = e_norms[:,newaxis] + e_norms[newaxis,:] \
        - 2 * dot(expr, expr.T)
    sims = - cluster_dists
    inp_dicts = []
    percentiles = [.01]#[.01,1.,10.,50.,75.]
    for p in percentiles:
        inp_dicts.append(dict(similarities = sims,
                              self_similarity = percentile(sims.flatten(),p)))
    eyeball = bsub.eyeball(os.path.abspath(inspect.stack()[0][1]), 
                           ['test_bsubfun'], 
                           inp_dicts)
    return eyeball



def usage():
  print '''
usage: btol.py run_id
1
Run a batch process on BTOL with inputs stored in 
data/batch/inputs/{run_id}.inp in pickle serial.
'''
  exit(1)

if __name__ == '__main__':
    run_id = sys.argv[2] if len(sys.argv) > 2 else 0
    run_func = globals()[sys.argv[1]]
    input_dict = bsub.load_inp(run_id)
    output_dict = run_func(input_dict, run_id)
    bsub.save_out( output_dict, run_id)
    
