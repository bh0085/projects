#!/usr/bin/env python
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
    
