#!/usr/bin/env python
'''
bsub_clusters.py

So far, mostly testing routines for my python-matlab interface over bsub.

local...  Spawns a single bsub process from a local machine
          Returns a launcher whose outputs and statuses can be queried.

remote... Is called by local_ remotely. Spawns a bunch of threads and
          exports the output to the local machine. Returns a dictionary
          containing the filepath of the exported data once subprocesses 
          have all run to completion.
'''
import subprocess as spc
import compbio.utils.bsub as bsub
import compbio.utils.bsub_utils as butils
import compbio.config as config
import compbio.utils.bs_macros as bsm


import os, sys, inspect, pipes
import scipy.io as sio
from numpy import *


#A few batch processes.
def bic_clustering(run_id):
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
    input_dict = butils.load_data(run_id, 'input')
    return bsm.runmat('ap_max_bic', input_dict, run_id)

def test_bsubfun(run_id):
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
    input_dict = butils.load_data(run_id,'input')
    return bsm.runmat('ap_frompy', input_dict, run_id)


#Local launchpoint for the batch scripts
def launcher(input_dicts, host = None, name = None):
    scriptfile = os.path.abspath(inspect.stack()[0][1])
    scriptroot = 'prog'
    func = 'remote_make_tests'
    if name == None:
        run_id = 'bcl_launch0'
    else:
        run_id = name
    launcher = bsub.local_launcher(scriptfile,
                                   scriptroot,
                                   func = func,
                                   input_dicts = input_dicts,
                                   run_id = run_id,
                                   host = host)
    return launcher
                                 
#Remote launchpoint for bsub.
def remote_make_tests(run_id):
    '''
the idea is that this function will queue up the batch jobs
and submit them with bsub. Using eyeball, it will then wait
until all jobs are done and when they are, export output back to
gliese.

inputs:
  run_id

output:
  the datapath (same for local and remote) of data output from 
  threads.
'''
    test = False
    if test:
        mirnaf = os.path.join(os.path.dirname(inspect.stack()[0][1]), 'miRNA.mat')
        mirna = sio.loadmat(mirnaf)
        expr = mirna['expression']
        e_norms = sum(expr**2,1)
        cluster_dists = e_norms[:,newaxis] + e_norms[newaxis,:] \
            - 2 * dot(expr, expr.T)
        sims = - cluster_dists
        inp_dicts = []
        percentiles = logspace(-2,1.99,3)
        for p in percentiles:
            inp_dicts.append(dict(similarities = sims,
                                  self_similarity = percentile(sims.flatten(),p)))
    else:
        inp_dicts = butils.load_data(run_id, 'input')
    
    eyeball = bsub.eyeball(run_id,
                           os.path.abspath(inspect.stack()[0][1]), inp_dicts,
                           func = 'test_bsubfun',
                           name = run_id+'_test_', 
                           mem = 2)
                          
  
    eyeball.launch()
    eyeball.await()
    eyeball.package()
    eyeball.complete()


def usage():
  print '''
usage: bsub_clusters function run_id

Call function with run_id.
'''
  exit(1)

if __name__ == '__main__':
    run_id = sys.argv[2]
    run_func = globals()[sys.argv[1]]
    output_dict = run_func(run_id)
    if output_dict == None:
        output_dict = {'blank':'Nothing output in call to {0}'.\
                           format(sys.argv[1])}
    butils.save_data( output_dict, run_id, 'output')
    
