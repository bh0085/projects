#!/usr/bin/env python
'''
Remote usage examples:

bsub -q compbio-week -o ${HOME}/bsub_tree_rfam.log ${PROGRAMMING_PATH}/projects/RNAfoldz/tree_scripts.py 'bsub_all' 'bsub_tree_all'
bsub -q compbio-week -o ${HOME}/bsub_tree_ribos.log ${PROGRAMMING_PATH}/projects/RNAfoldz/tree_scripts.py 'bsub_riboswitches' 'bsub_tree_riboswitches'

'''
import sys, os, inspect
import compbio.utils.bsub as bsub
import compbio.utils.bsub_utils as bsu
import tree_utils as tutils

def bsub_all(run_id):
        inp_dicts = []
        for r in  range(0,1493):
                rfid = 'RF{0:05}'.format(r)
                inp_dicts.append(dict(family = rfid,
                                      inp_run_id = 'FA_{0}'.format(rfid),
                                      run_id = 'FA_tree_{0}'.format(rfid)))
	eyeball = bsub.eyeball(run_id, 
			       os.path.abspath(inspect.stack()[0][1]),
			       inp_dicts,
			       func = 'run_treebuilder',
			       mem = 3)
	eyeball.launch()
	return dict(cmds=eyeball.cmds)

def bsub_riboswitches(run_id):
        inp_dicts = []
        for k,v in rutils.switch_dicts().iteritems():
                rfid = 'RF{0:05}'.format(v)
                inp_dicts.append(dict(family = rfid,
                                      inp_run_id = 'RS_{0}'.format(rfid),
                                      run_id = 'RS_tree_{0}'.format(rfid)))
	eyeball = bsub.eyeball(run_id, 
			       os.path.abspath(inspect.stack()[0][1]),
			       inp_dicts,
			       func = 'run_treebuilder',
			       mem = 3)
	eyeball.launch()
	return dict(cmds=eyeball.cmds)

def run_treebuilder(run_id):
  '''
So that output files can be located by family, 
run_ids should have a consistent format: 

RS_{familyname}

In any event, the family name is stored in the
input dict so there will be no problems running the
script is the run_id does not conform to this format.
'''
  inp_dict = bsub.load_data(run_id, 'input')
  inp_run_id = inp_dict['inp_run_id']
  rfid = inp_dict['family']

  output = tutils.run(rfid,run_id, inp_run_id, reset = True)
  return output



def usage():
  print '''
usage: rna_ali2d function run_id

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
    bsu.save_data( output_dict, run_id, 'output')
    
