default_name = 'proj0'
last_run = ''
default_window = 200

default_genname = 'avida.cfg'
default_evename = 'events.cfg'
default_envname = 'environment.cfg'
default_orgname = 'default-heads.org'
default_ananame = 'analyze.cfg'
default_insname = 'instset-heads.cfg'

default_task_grid_frequency = 2
default_printrate = 'g 0:1'

default_framerate =10

import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import sys, time, os, re
import cfg_avida as ca
import avida_plots as ap
import avida_utils as au
import avida_presets as presets



#make config now incorporates the params file so that we can know
#in advance what kind of data to print out.
def make_cfg(plot_params,
             name =default_name,
             geo = 'square',
             mut = .001,
             delta = .0005,
             dim = 20):

    if not os.path.isdir(name):
        os.mkdir(name)
    default_file = open(os.path.join(name,'avida.cfg'))
    outfile = open(os.path.join(name,'avida.cfg.auto'),'w')
    outfile.write('#avida.cfg generated by python.')
    lines = default_file.readlines()
    

    auto_dir =  os.path.join(name,os.path.join('autogen'))
    if not os.path.isdir(auto_dir):
        os.mkdir(auto_dir)
    eve_fname = ca.make_eve(name,default_evename)
    org_fname = ca.make_org(name,default_orgname)
    ana_fname = ca.make_ana(name,default_ananame)
    env_fname = ca.make_env(name,default_envname)
    ins_fname = ca.make_ins(name,default_insname)

    computed = presets.differential_mut(name,
                                        mut = mut,
                                        geo_name = geo, 
                                        delta = delta,
                                        dim = dim,
                                        eve_fname = eve_fname,
                                        org_fname = org_fname,
                                        ana_fname = ana_fname,
                                        env_fname = env_fname,
                                        ins_fname = ins_fname,
                                        gen_name = default_genname)
    for pp in plot_params:
        pp['computed'] = computed

    print
    print 'setting print params: '
    print

    for p in plot_params:
        if p['type'] == 'ts':
            colstr = ''
            for col in p.get('command_params',[]):
                colstr += col + ' '
            printstr = p['update']+' '+p['command'] +' '+p['fname'] + ' '+colstr
            ca.alter_eve(eve_fname, printstr)
        elif p['type'] == 'grid':
            printstr = p['update'] + ' ' + p['command'] + ' ' + p['fname']
            ca.alter_eve(eve_fname, printstr)
    
    f = open(eve_fname).read()
    print f
    
def runblank(name = default_name):
    proc = spawn(name, silent = True)
    com = proc.communicate()[0]
    return 

def runview(plot_params,
            kill_at_end = True,
            name = default_name,
            silent = True):



    #Now, we use a name just to locate the avida.cfg file.
    #Soon I'll probably eliminate this and just make send
    #all of the autogen prefs to a single folder.
    proc = spawn(name, silent = silent)
    time.sleep(.25)
    max_loops = 1000
    loops = 0

    while 1:
        try:
            for p in plot_params:
                ap.make_plot(p)

                
            time.sleep(1./default_framerate)
            if proc.poll() != None:
                print proc.stdout.read()
                print proc.poll()
                print
                print 'Avida finished running'
                print
                break

            loops+=1
            if loops>max_loops:
                print
                print 'Visual loop maxed out frame count.'
                print '...killing avida'
                print
                proc.kill()
                print 'avida killed'
                print
                print '...exiting'
                print 
                break
        except KeyboardInterrupt, e:
            print 'process interrupted:'
            print 'killing avida'
            proc.kill()
            break
        
        #except Exception, e:
        #    if e.args[0] == 'kill':
        #        print 'process interrupted by internal exception:'
        #        print 'killing avida'
        #        proc.kill()
        #        break
        #    else: raise e



    #cols, data = au.parse_printed( filename = filename)
    #return cols, data


def spawn(name=None, silent = True):
    import subprocess

    genesis = os.path.join(os.path.join(name,'autogen'),'avida.cfg.auto')

    silent = True
    if silent:
        proc = subprocess.Popen('avida -c '+genesis + '', shell = True,  \
                                    stdout=open('/dev/null','w'))
    else:
        proc = subprocess.Popen('avida -c '+genesis + '', shell = True)
    return proc




def run_batch():
    mut = .002
    delta = .0005
    geometries = ['star','square']
    dim =  20
    iters = 50
    import compbio.config as cfg
    
    root = cfg.dataPath('avida_runs')
    if not os.path.isdir(root):
        os.mkdir(root)
    os.chdir( root)
    exec_subdirs = ['{0}_{1}'.format(i,j) 
                    for i in geometries 
                    for j in range(iters)]
    for e in exec_subdirs:
        if not os.path.isdir(os.path.join(root,e)):
            os.mkdir(os.path.join(root,e))
        os.chdir(os.path.join(root,e))
        p = get_params(psets['lots'])
        make_cfg(geo = e.split('_')[0],
                 mut = mut,
                 delta = delta,
                 dim = dim)
        
        
    
    
    

def match_names(cols, regex):
    cols_out = []
    for c in cols:
        if re.search(regex, c):
            cols_out.append(c)
    return cols_out

#unused options include: lin2D, res, all,'tasks'
psets = {'lots':['res','all','lin2D','tasks','tasks2D'],
         '2grids':['lin2D','tasks2D']}
def get_params(names = psets['2grids']):    
    all_plot_params = []
    
    
    for i in range(len(names)):
        plot_params= {}
        name = names[i]
        plot_params['name'] = name


        if '2D' in name:
            plot_params['type'] = 'grid'
        else:
            plot_params['type'] = 'ts'

        if name =='cData2D':
            filename = 'celldata_grid.data'
            plot_params['command'] = 'DumpCellDataGrid'
            plot_params['update'] = default_printrate

        if name =='lin2D':
            filename = 'lineage_grid.data'
            plot_params['command'] = 'DumpLineageGrid'
            plot_params['update'] = default_printrate
            
        if name =='res':
            filename = 'resource_printed.dat'
            plot_params['command'] = 'PrintResourceData'
            plot_params['update']= default_printrate

        if name =='tasks2D':
            filename = 'tasks_grid.data'
            plot_params['command'] = 'DumpTaskGrid'
            plot_params['update'] = 'g 0:'+str(default_task_grid_frequency)
                       
        if name == 'all':
            filename = 'custom_data.dat'
            plot_params['command'] = 'PrintData'
            plot_params['update'] = default_printrate
            cols =  ['ave_fitness']
            plot_params['command_params'] = [','.join(cols)]
            widlambda = lambda x : 1
            plot_params['widths'] =  [ widlambda(cols[j]) for j in range(len(cols))]

        if name =='cdata':
            filename = 'celldata.dat'
            plot_params['command'] = 'PrintCellData'
            plot_params['update'] = default_printrate

        if name == 'tasks':
            filename = 'tasks.dat'                                   
            plot_params['command'] = 'PrintTasksData'
            plot_params['update'] = default_printrate
        if name == 'tasks_deme':
            filename = 'tasks_deme.dat'                                   
            plot_params['command'] = 'PrintDemeTasksData'
            plot_params['update'] = default_printrate

        plot_params['window'] = default_window
        plot_params['figure'] = i
        plot_params['fname'] = filename
        all_plot_params.append(plot_params)

    return all_plot_params

