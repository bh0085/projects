import os, re
import cb.utils.plots as myplots
import cb.utils.colors as mycolors
from numpy import *
import numpy as np

def data_dirs(root = '/data/avida_runs'):
    ddirs = []
    for r,dirs,files in os.walk(root):
        for d in dirs:
            if d == 'data':
                ddirs.append(os.path.join(r,d))
    
    return ddirs

def dir_infos(dirs):
    infos = []
    for d in dirs:
        run_name = os.path.split(os.path.split(d)[0])[1]
        splits = run_name.split('_')
        if len(splits) != 5:
            continue
        chars   = {'name':splits[1],
                   'mut':float(splits[2]),
                   'dim':int(splits[3]),
                   'iter_num':int(splits[4]),
                   'dir':d}
        infos.append(chars)
    return infos


#choose a binning scheme and the fashion in which to
#deal with extra variables (fix them or integrate over them)
#
#then return binned directories for extraction of binned data.
def info_groups_1d(infos,
                   bins1d = 'name',
                   **kwargs
                   ):
    #default values
    params = dict(
        dim = 40,
        iter_num = 'any',
        name = 'square',
        mut = .002
        )
    #kwargs may only contain params.
    assert len(set(kwargs.keys()).difference(params.keys())) == 0 
    params.update(kwargs)
    params[bins1d] = 'any'

    binvals = sorted(set([i[bins1d] for i in infos]))
    out = dict([(k,[]) for k in binvals])

    for i in infos:
        skip = False
        for k, v in params.iteritems():
            if v == 'any': continue
            if i[k] != v : skip = True
        if not skip:
            out[i[bins1d]].append(i)
    
    return out
        
def batch_dimensions(igroups):
    lens = []
    gl = len(igroups)
    for i,e in enumerate(igroups.iteritems()):
        k,v = e
        for j in v:
            d0 = j['dir']
            fname = os.path.join(d0,'tasks.dat')
            if not os.path.isfile(fname): continue
            tfile = open(fname)
            tdat = tfile.readlines()
            tfile.close()
            
            infolines = '\n'.join(tdat[0:4])
            l = l0 = 4
            while 1:
                if tdat[l][0] != '#': break
                l+=1

            col_lines = tdat[l0:l]
            task_count = len(col_lines)
            ctuples =[re.compile('(\d+): (\S+)').search(line).groups() 
                      for line in col_lines]
            cmaps =dict([(int(e[0]) -1,e[1]) for e in ctuples])
            data_array = array([[float(e) for e in line.strip().split(' ')] 
                                for line in tdat[l+1:]])
            lens.append(len(data_array))
    return max(lens), task_count
            


def binned_data_1d(igroups, dimensions):
    '''
takes igroups as computed from info_groups_1d.

Also, matrix dimensionality information that comes
from batch_stats.
'''
    max_gens, n_tasks =dimensions 

    lens = []
    gl = len(igroups)
    datum = []
    descriptors = []
    for i,e in enumerate(igroups.iteritems()):
        k,v = e
        descriptors.append((k,v))
        datum.append([])
        for j, info in enumerate(v):
            d0 = info['dir']
            fname = os.path.join(d0,'tasks.dat')
            if not os.path.isfile(fname): continue
            
            tfile = open(fname)
            tdat = tfile.readlines()
            tfile.close()
            
            infolines = '\n'.join(tdat[0:4])
            l = l0 = 4
            while 1:
                if tdat[l][0] != '#': break
                l+=1

            col_lines = tdat[l0:l]
            ctuples =[re.compile('(\d+): (\S+)').search(line).groups() 
                      for line in col_lines]
            cmaps =dict([(int(e[0]) -1,e[1]) for e in ctuples])


            data_array = array([[float(e) for e in line.strip().split(' ')] 
                                for line in tdat[l+1:]])

            if not len(data_array) == max_gens: 
                continue

            datum[i].append(data_array)
            
    return descriptors, datum, cmaps

def show_binned_data_1d(descriptors, datum, cmaps):
    '''
Grab metadata and parsed grids from binned_data_1d and
plot them.
'''
    f = myplots.fignum(1, (6,8))
    gl = len(descriptors)
    dshapes = [shape(d) for d in datum]
    ntasks = dshapes[0][-1]
    task_colors = mycolors.getct(ntasks)
    
    
    for i,e in enumerate(datum):
        ax = f.add_subplot('{0}1{1}'\
                               .format(gl,i + 1))
        ax.set_title(descriptors[i][0])
        
        
        
        sums =np.sum(e, 0)
        for i,task in enumerate(sums.T):
            p = ax.plot(log(2 + task[::20]), 
                    color =task_colors[i],
                    label = cmaps[i])
        
    
          
def fixation_props(descriptors, datum, cmaps):
    '''
Compute first appearance, 10% fixations, 50% fixations for 
each task.
'''
    nt = len(cmaps)
    tasks = [i for i in range(1,nt)]
    dprops = []

    keys=['1','5','10','10%','25%','50%']
    for i,d in enumerate(datum):
        dprops.append({})
        for t in tasks:
            dprops[i][t] = dict([(k, []) for k in keys]) 
    
    for i,d in enumerate(datum):
        #max number (not including updates)
        max_pop = max([ np.max(e[:,1:]) for e in d])
        for t in tasks:
            tprops = dprops[i][t]
            for k, run in enumerate(d):
                run_task = run[:,t]
                nz = nonzero(greater(run_task,0))[0]
                tprops['1'].append((nz[0] if len(nz) > 0 else -1
                                    ,run[nz[0] if len(nz) > 0 else -1
                                         ,0]))
                tprops['5'].append((nz[4] if len(nz) > 4 else -1
                                    ,run[nz[4] if len(nz) > 4 else -1
                                         ,0]))
                tprops['10'].append((nz[9] if len(nz) > 9 else -1
                                    ,run[nz[9] if len(nz) > 9 else -1
                                         ,0]))
                
                
                perc10=nonzero(greater(run_task,max_pop * .1))[0]
                tprops['10%'].append((perc10[0],run[perc10[0],0])
                                     if len(perc10) > 0 else (-1,-1))
                perc25=nonzero(greater(run_task,max_pop * .25))[0]
                tprops['25%'].append((perc25[0],run[perc25[0],0]) 
                                     if len(perc25) > 0 else (-1,-1))
                perc50=nonzero(greater(run_task,max_pop * .5))[0]
                tprops['50%'].append((perc50[0],run[perc50[0],0])
                                     if len(perc50) > 0 else (-1,-1))
    
    
    return dprops

def mutation_fixations(infos):
    mut_vals =  [0.0005, 0.001, 0.002, 0.003, 0.004]
    fixations = {}
    for m in mut_vals:
        print 'computing groups -- m = {0}\n'.format(m)
        igroups = info_groups_1d(infos,
                                 bins1d = 'name',
                                 mut = m,
                                 dim = 40)
        print 'dimensions,descriptors...\n'
        print 'dims'
        dimensions = batch_dimensions(igroups)
        print 'descr'
        descriptors, datum, cmaps = binned_data_1d(igroups, dimensions)
        print 'fixation'
        fixations[m] = fixation_props(descriptors, datum, cmaps)
        
    return fixations
                                 
def fixation_stats(all_fixations):
    out = {}
    ref_mut = .002
    ref_name_idx= 1
    
    times_10p = {}
    for t,v in all_fixations[ref_mut][ref_name_idx].iteritems():
        mean_10 = mean([e[1] for e in v['10%'] if e[1] != -1])
        times_10p[t] = mean_10
    out['fix_ages'] = times_10p

    ages = out['fix_ages']
    age_tuples = [(k,v) for k, v in ages.iteritems()]
    ages_ranked = sorted(age_tuples, key = lambda x: x[1])
    
    out['ranked_ages'] = ages_ranked

    return out


def show_fixations(all_fixations,cmaps):

    fstats = fixation_stats(all_fixations)
    ranked_tasks = fstats['ranked_ages']
    
    names = ['s_star','torus','torus_repl']

    #argsort(ages, key = lambda)

    mut_vals =  [0.0005, 0.001, 0.002, 0.003, 0.004]
   
    xs = []
    ys = []
    cs = []

    name_colors = mycolors.getct(3)
    task_colors = mycolors.getct(9)

    all_deltas = zeros((len(mut_vals),
                        len(all_fixations.values()[0]),
                        8))
    all_fracs =  zeros((len(mut_vals),
                        len(all_fixations.values()[0])))
    all_counts =  zeros((len(mut_vals),
                        len(all_fixations.values()[0])))
    mut_spool = []
    
    for i, m in enumerate(mut_vals):
        fix_map = all_fixations[m]
        for j, e in enumerate(fix_map):
            name_c= name_colors[j]
            
            task_10ptime = array([ [item[1] for item in e[rt[0]]['1']] 
                                   for rt in ranked_tasks])
            idxs_allowed = nonzero(greater(np.min(task_10ptime, 0),0))[0]
            frac_allowed =float( len(idxs_allowed))/ shape(task_10ptime)[1]
            
            '''roll deltas for all sxsful completions'''
            if len(idxs_allowed )!= 0:
                nrml = task_10ptime[:,idxs_allowed]
                #nrml = 1
                deltas = np.mean((roll(task_10ptime[:,idxs_allowed],-1) \
                                       - task_10ptime[:,idxs_allowed]) / \
                                    nrml ,1)
                all_deltas[i,j,:] = deltas[:-1]
                all_counts[i,j] = len(idxs_allowed)
                all_fracs[i,j] = frac_allowed
                mut_spool.append({'tuple':(i,j),
                                  'mut':m,
                                  'name': j})

            for k, e2 in enumerate(e.iteritems()):
                t,v = e2
                task_c = task_colors[k]
                p10_times = [item[1] 
                             for item in v['5'] if item[0] != -1]
                n= len(p10_times)
                these_x = (zeros(n) + i ) + random.uniform(size = n)/3
                xs.extend(these_x)
                ys.extend(p10_times)
                cs.extend([task_c] * n)
    f = myplots.fignum(2,(6,6))
    ax = f.add_subplot(111)
    ax.set_title('fixation times for all tasks')
    ax.set_xlabel('mutation rate')
    ax.set_ylabel('fixation time')
    #ax.scatter(xs, ys, 20, color = cs,alpha = .4)
                

    f2 = myplots.fignum(3, (8,8))
    ax = f2.add_axes([.3,.3,.6,.6])
    ax.set_title('fixation time (fold change over previous tasks)')
    ax.set_xlabel('task')
    ax.set_ylabel('condition')
    
    xlabels = [(cmaps[e[0][0]],cmaps[e[1][0]]) 
               for e in zip(ranked_tasks, roll(ranked_tasks,-1,0))][0:-1]

    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(['{0} -> {1}'.format(*xl) for xl in xlabels],
                       rotation = -90,
                       va = 'top',
                       ha = 'left')


    rows = []
    labels = []
    for ms in sorted(mut_spool, key = lambda x:x['name']):
        tup = ms['tuple']
        rows.append(all_deltas[tup[0],tup[1],:])
        ct  = all_counts[tup]
        frac =all_fracs[tup]
        mut = ms['mut']
        labels.append('{0}: mut rate={1};n={2}'.\
                          format(names[ms['name']],mut,int(ct),frac)
                      )
    

    im = ax.imshow(rows,interpolation = 'nearest')
    f2.colorbar(im)    

    
    ax.set_yticks(range(len(mut_spool)))
    ax.set_yticklabels(labels)

    f2.savefig(myplots.figpath('graph_evolution_acceleration.pdf'))
