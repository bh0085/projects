#!/usr/bin/env python
import utils as rutils
import compbio.utils.bsub as bsub
import compbio.utils.bsub_utils as bsu
import os, inspect
import sys
def make_ribostructs():    
    for k,v in rutils.switch_dicts().iteritems():
        rfid = 'RF{0:05}'.format(v)
        savename = 'Riboswitch_list_{1}_{0}'.format(rfid,k)
        structs, rutils.family_clustered_suboptimals(rfid =rfid,savename = savename, draw = True)

def bsub_riboswitches(run_id):
        inp_dicts = []
        for k,v in rutils.switch_dicts().iteritems():
                rfid = 'RF{0:05}'.format(v)
                savename = 'Riboswitch_list_{1}_{0}'.format(rfid,k)
                inp_dicts.append(dict(family = rfid,
                                      run_id = 'RS_{0}'.format(rfid),
                                      savename = savename))
	eyeball = bsub.eyeball(run_id, 
			       os.path.abspath(inspect.stack()[0][1]),
			       inp_dicts,
			       func = 'run_structmaker',
			       mem = 3)
	eyeball.launch()
	return dict(cmds=eyeball.cmds)

def run_structmaker(run_id):
  '''
So that output files can be located by family, 
run_ids should have a consistent format: 

RS_{familyname}

In any event, the family name is stored in the
input dict so there will be no problems running the
script is the run_id does not conform to this format.
'''
  inp_dict = bsub.load_data(run_id, 'input')
  fam = inp_dict['family']
  savename = inp_dict['savename']

  structs, energies, seq = rutils.family_clustered_suboptimals(\
          fam,
          savename = savename,
          draw = False)

  output = {'family':fam,
            'structs':structs,
            'energies':energies,
            'seq':seq}

  return output

def runswitches():
	all_outs = []
	for k,ofs in switch_dicts().iteritems():
		print k
		out = get_consensus(reset = True,
				    run_id = 'riboswitch_{0}_RF{1}'.format(k,ofs))
		all_outs.append(out)
	for out in all_outs:
		rplots.show_output(out, show = 'conservation')
		rplots.show_output(out, show = 'embeddings')
	return all_outs
def runmany(run_id):
	print 'TESTING WITH A LIMITED RANGE OF FAMILIES'
	inp_dicts = [dict([('ofs',r)]) for r in range(0,1493)]
	eyeball = bsub.eyeball(run_id, 
			       os.path.abspath(inspect.stack()[0][1]),
			       inp_dicts,
			       func = 'run',
			       name = 'ra2_runs_',
			       mem = 3)
	eyeball.launch()
	return dict(cmds=eyeball.cmds)
def run(run_id):
	data = bsu.load_data(run_id, 'input')
	ofs = data['ofs']
	outputs = get_consensus(ofs, 
			       run_id = run_id,
			       reset = True)
	return(outputs)

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
    
