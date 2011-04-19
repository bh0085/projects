import subprocess as spc

def make_tests():
    

def test_bsubfun():
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

if __name__ = '__main__':
    if len(sys.arv < 3): usage()
    run_id = sys.argv[2]
    run_func = globals()[sys.argv[1]]
    input_dict = bsub.load_inp(run_id)
    output_dict = run_func(input_dict, run_id)
    bsub.save_out( output_dict, run_id)
    
