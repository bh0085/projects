#!/usr/bin/env python

import subprocess as spc
import compbio.utils.bsub as bsub
import inspect
import os, sys

def make_tests():
    eyeball = bsub.eyeball(os.path.abspath(inspect.stack()[0][1]), 
                           ['test_bsubfun'], 
                           [{} for i in range(5)])
    return eyeball
def test_bsubfun(input_dict, run_id):
    sub = spc.Popen('find', shell =True, stdout = spc.PIPE).\
        communicate()[0]
    out_dict = dict(output = sub)
    
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
    
