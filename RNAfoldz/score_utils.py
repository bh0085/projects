#!/usr/bin/env python
import utils as rutils
import plots as rplots
import compbio.utils.bsub_utils as bsu

from numpy import *
import numpy as np, itertools as it

import matplotlib.pyplot as plt

import compbio.utils.plots as myplots
import compbio.utils.colors as mycolors
import compbio.config as cfg
import pickle

figsize = (8,8)
figtype = 'ps'
figfile = cfg.dataPath('figs/gpm2/pt3_scores/{{0}}.{0}'.format(figtype))

draw_many = True
draw_single = True


def ribo_struct_outfile(rfid):
    si = 'RS_{0}'.format(rfid)
    return bsu.load_data(si, 'output')



def save_muts_structs(out,out_tree):
    ofile = open(cfg.dataPath('RNAfoldz/out.pickle'),'w')
    otfile = open(cfg.dataPath('RNAfoldz/out_tree.pickle'), 'w')
    pickle.dump(out, ofile)
    pickle.dump(out_tree,otfile)
    ofile.close(), otfile.close()

def get_muts_structs():
    ofile = open(cfg.dataPath('RNAfoldz/out.pickle'),'w')
    otfile = open(cfg.dataPath('RNAfoldz/out_tree.pickle'),'w')
    ofile.close; otfile.close()
    return pickle.load(out, ofile), pickle.load(out_tree,otfile)

def get_tree(rfid, fam_type = 'all'):
    if fam_type == 'riboswitch':
        rname = rfid
        si ='FA_{0}'.format(rfid) 
        ti = 'FA_tree_{0}'.format(rfid)

    elif fam_type == 'all':
        rname = rfid
        si = 'FA_{0}'.format(rfid)
        ti = 'FA_tree_{0}'.format(rfid) 


    print 'Loading family for {0}'.format(rname)

    try:
        structs = bsu.load_data(si, 'output')
        trees   = bsu.load_data(ti, 'output')
        print 'Success! Analyzing tree output'
    except Exception, e:
        print 'Failure! Did I make a booboo?'
        return None
    
        

    str_esrt = argsort(structs['energies'])[::-1]
    #SORT STRUCTURES IN DECREASING ORDER OF ENERGY (MATCH TREES)
    structs['structs'] = [structs['structs'][j] for j in str_esrt]
    structs['energies'] = [structs['energies'][j] for j in str_esrt]
    
    mc, tc, sc = {},{},{}
    sc['energies'] = structs['energies']
    sc['structs'] = structs['structs']
    sc['seq'] = structs['seq']
    for j, t in enumerate(trees):
        if t == None:continue
        mc[j] = {}
        tc[j] = {}
        for idx in range(len(t['structs'])):
            t_infos = t['structs']
            t_str = t_infos[idx]['pairs']
            s_str = structs['structs'][idx]
            e = structs['energies'][idx]
            
            t_times =t_infos[idx]['times']
            t_muts  =t_infos[idx]['muts']
            
            frac_resolved = t_times['total']  /\
                (t_times['total_incl_unresolved'])
            frac_paired   = t_times['paired'] /\
                (t_times['unpaired'] + t_times['paired'])

            n_2cons = t_muts['comp']
            n_1cons = t_muts['wob'] 
            n_0cons = t_muts['ucom']
            n_pluscons = t_muts['reco']
            n_nocons =   t_muts['bbad']

            frac_silent = (n_2cons+n_1cons+n_pluscons)/\
                (n_0cons + n_nocons+\
                     n_2cons+n_1cons+n_pluscons)
            
            frac_double = (n_2cons)/\
                (n_2cons+n_1cons)
                                      
            frac_destructive=(n_0cons)/\
                (n_2cons+n_1cons+n_0cons)

            total_muts = (n_0cons + n_nocons+\
                                   n_2cons+n_1cons+n_pluscons)
            total_silent = (n_2cons+n_1cons)
            total_pair_mutants = (n_2cons+n_1cons+n_0cons)

            tc[j][idx] = dict(frac_resolved = frac_resolved,
                         frac_paired = frac_paired,
                         total_time = array([t_times['total_incl_unresolved']]*len(t_times['total'])),
                         total_time_res = t_times['total'])
            mc[j][idx] = dict(frac_silent = frac_silent,
                         frac_double = frac_double,
                         frac_destructive = frac_destructive,
                         total_muts = total_muts,
                         total_silent = total_silent,
                         total_pair_mutants = total_pair_mutants)

            
    print '''Done!
Results:
  {0} subtrees computed.

'''.format(len(mc))
    return mc, tc, sc

def check_trees(fam_type = 'riboswitch'):
    if fam_type == 'riboswitch':
        sdicts = rutils.switch_dicts()
        rfids = [ 'RF{0:05}'.format(n) for n in sdicts.values()]
        names = sdicts.keys()
        struct_ids = ['RS_{0}'.format(rfid) for rfid in rfids]
        tree_ids = ['RS_tree_{0}'.format(rfid) for rfid in rfids]
    elif fam_type == 'all':
        rfids = [ 'RF{0:05}'.format(n) for n in range(0,1493)]
        names = rfids
        struct_ids = ['FA_{0}'.format(rfid) for rfid in rfids]
        tree_ids = ['FA_tree_{0}'.format(rfid) for rfid in rfids]

    switch_muts = {}
    switch_times= {}
    switch_structs= {}

    for i, rname in enumerate(names):
        print 'Loading family for {0}'.format(rname)
        si, ti = zip(*[struct_ids, tree_ids])[i]
        try:
            structs = bsu.load_data(si, 'output')
            trees   = bsu.load_data(ti, 'output')
            print 'Success! Analyzing tree output'
        except Exception, e:
            print 'Failure! Did I make a booboo?'
            continue
        
            

        str_esrt = argsort(structs['energies'])[::-1]
        #SORT STRUCTURES IN DECREASING ORDER OF ENERGY (MATCH TREES)
        structs['structs'] = [structs['structs'][j] for j in str_esrt]
        structs['energies'] = [structs['energies'][j] for j in str_esrt]
        
        mc, tc, sc = {},{},{}
        sc['energies'] = structs['energies']
        sc['structs'] = structs['structs']
        sc['seq'] = structs['seq']
        for j, t in enumerate(trees):
            if t == None:continue
            mc[j] = {}
            tc[j] = {}
            for idx in range(len(t['structs'])):
                t_infos = t['structs']
                t_str = t_infos[idx]['pairs']
                s_str = structs['structs'][idx]
                e = structs['energies'][idx]
                
                t_times =t_infos[idx]['times']
                t_muts  =t_infos[idx]['muts']
                
                frac_resolved = t_times['total']  /\
                    (t_times['total_incl_unresolved'])
                frac_paired   = t_times['paired'] /\
                    (t_times['unpaired'] + t_times['paired'])

                n_2cons = t_muts['comp']
                n_1cons = t_muts['wob'] 
                n_0cons = t_muts['ucom']
                n_pluscons = t_muts['reco']
                n_nocons =   t_muts['bbad']

                frac_silent = (n_2cons+n_1cons+n_pluscons)/\
                    (n_0cons + n_nocons+\
                         n_2cons+n_1cons+n_pluscons)
                
                frac_double = (n_2cons)/\
                    (n_2cons+n_1cons)
                                          
                frac_destructive=(n_0cons)/\
                    (n_2cons+n_1cons+n_0cons)

                total_muts = (n_0cons + n_nocons+\
                                       n_2cons+n_1cons+n_pluscons)
                total_silent = (n_2cons+n_1cons)
                total_pair_mutants = (n_2cons+n_1cons+n_0cons)

                tc[j][idx] = dict(frac_resolved = frac_resolved,
                             frac_paired = frac_paired,
                             total_time = array([t_times['total_incl_unresolved']]*len(t_times['total'])),
                             total_time_res = t_times['total'])
                mc[j][idx] = dict(frac_silent = frac_silent,
                             frac_double = frac_double,
                             frac_destructive = frac_destructive,
                             total_muts = total_muts,
                             total_silent = total_silent,
                             total_pair_mutants = total_pair_mutants)

                
        print '''Done!
Results:
  {0} subtrees computed.

'''.format(len(mc))

        switch_muts[rname] = mc
        switch_times[rname] = tc
        switch_structs[rname] = sc
    return switch_muts, switch_times, switch_structs
    
def show_trees(switch_muts, switch_times, switch_structs):
    for rname,muts, times, structs in zip(switch_times.keys(),
                                          switch_muts.values(), 
                                          switch_times.values(), 
                                          switch_structs.values(),):
        if not rname[0:2] == 'RF':
            rtype = 'switch'
            rfid = 'RF{0:05}'.format(rutils.switch_dicts()[rname])
        else:
            rfid = rname
            rtype = 'all'

        print 'Showing: {0} '.format(rfid)
        #if rname != 'SAM_alpha':continue
        show_paired_v_energy(rname,rfid,muts,times,structs,rtype)
        return

def metric(mvals, tvals,
           mtype = 'frac_silent'):
    '''
Metrics:
  time metrics:
    frac_paired

  mut metrics:
    frac_silent

  stat metrics:
    frac_ug
    frac_good

Returns:
  value: a float or a vector depending on whether the metric
         has been computed for the entire structure or for 
         the average
  idxs:  idxs where the metric has been computed

'''
    npairs = len(tvals['frac_paired'])
    
    if mtype in ['frac_ug']:
        idxs_good = range(npairs)
    elif mtype in ['frac_silent', 'frac_double']:
        idxs_good = nonzero( isfinite(mvals['frac_silent']))[0]
    elif mtype in ['frac_paired','n_comp']:
        idxs_good = nonzero(tvals['total_time_res'])[0]
    idxs_bad = [ i for i in range(npairs) if not i in idxs_good ]


    frac_good = float(len(idxs_good))/\
        (len(idxs_good)+len(idxs_bad))


    gap_ct = len(nonzero(equal(tvals['frac_resolved'],0))[0])
    ug_ct = len(nonzero(greater(tvals['frac_resolved'],0))[0])
    frac_ug = float(ug_ct) / (ug_ct + gap_ct)

    frac_good = float(len(idxs_bad))/\
        (len(idxs_good)+len(idxs_bad))

    pfrac = [tvals['frac_paired'][idx] for idx in idxs_good]
    sfrac = [mvals['frac_silent'][idx] for idx in idxs_good]
    dfrac = [mvals['frac_double'][idx] for idx in idxs_good]
    n_comp = array([mvals['frac_double'][idx]*mvals['total_silent'][idx] 
                    for idx in idxs_good])
    n_comp[nonzero(1-isfinite(n_comp))] = 0
    if mtype == 'frac_ug':
        return frac_ug, idxs_good,frac_good
    elif mtype == 'frac_silent':
        return sfrac, idxs_good,frac_good
    elif mtype == 'frac_double':
        return dfrac, idxs_good,frac_good
    elif mtype == 'frac_paired':
        return pfrac, idxs_good,frac_good
    elif mtype == 'n_comp':
        return n_comp, idxs_good, frac_good
                
    else: raise Exception('mtype {0} not implemented'.format(mtype))
    
def show_paired_v_energy(rname,rfid, all_muts, all_times, structs,rtype):
    if all_times == {}:
        return
    resolved_frac =  [ mean(list(it.chain(*[s_times['frac_resolved'] 
                                       for s_times in t_times.values()])))
                       for t_times in all_times.values()]
    total_lens = [mean(list(it.chain(*[s_times['total_time'] 
                                       for s_times in t_times.values()])))
                  for t_times in all_times.values()]
    total_lens_res = [mean(list(it.chain(*[s_times['total_time_res'] 
                                       for s_times in t_times.values()])))
                  for t_times in all_times.values()]    
    focus_tree = all_times.keys()[argmax(total_lens_res)]
    
    muts = all_muts[focus_tree]
    times = all_times[focus_tree]


    ns = len(muts.keys())
    s2 = dict(structs)
    s2['energies'] = s2['energies'][:ns]
    s2['structs'] = s2['structs'][:ns]
    structs = s2

    energies = structs['energies']


    
    f = myplots.fignum(3,figsize)

    xvals, yvals , pfracs, ugfracs = [], [], [], []
    for i, vals in enumerate(zip(muts.values(),times.values())):
        mvals, tvals = vals
        xvals.append( energies[i])
        
        frac_ug = metric(mvals, tvals, 'frac_ug')[0]
        pfrac,pinds,frac_good = metric(mvals, tvals, 'frac_silent')
        sfrac = metric(mvals, tvals, 'frac_silent')[0]
        

        ugfracs.append(frac_ug)
        pfracs.append(  mean(pfrac))
        yvals.append( mean(sfrac)*frac_good)
        
    colors = array(pfracs)
    colors = (colors - min(colors)) /(max(colors) -  min(colors))
    colors = colors[:,newaxis] * [0,1,0]

    ax = f.add_subplot(111)

    ax.scatter(xvals,yvals,array(ugfracs) * 200,
               color = colors)
    ax.set_ylabel('mutation score')
    ax.set_xlabel('free energy (-kCal)')
    ax.annotate('''Evaluated structures positioned by energy 
and a mutation based evolutionary score.
Color indicates fractional frequency of double mutants.
Radius indicates percentage of ungapped base pairs.''' , [0,1],xycoords = 'axes fraction',
                        xytext = [10,-10], textcoords='offset pixels',
                        va = 'top')
          
    myplots.padded_limits(ax, xvals, yvals, .2)
    f.savefig(figfile.format('{1}_frac_double_{0}'.format(rname,rtype)))
    f.clear()
            
    colors = array(pfracs)
    colors = (colors - min(colors)) /(max(colors) -  min(colors))
    colors = colors[:,newaxis] * [1,0,0]

    f = myplots.fignum(3,figsize)
    


    xvals, yvals , pfracs, ugfracs = [], [], [], []
    for i, vals in enumerate(zip(muts.values(),times.values())):
        mvals, tvals = vals
        xvals.append( energies[i])
        
        frac_ug = metric(mvals, tvals, 'frac_ug')[0]
        pfrac,pinds,frac_good = metric(mvals, tvals, 'frac_paired')
        sfrac = metric(mvals, tvals, 'frac_silent')[0]
        

        ugfracs.append(frac_ug)
        pfracs.append(  mean(pfrac)*frac_good)
        yvals.append( mean(sfrac)*frac_good)
        
    colors = array(pfracs)
    colors = (colors - min(colors)) /(max(colors) -  min(colors))
    colors = colors[:,newaxis] * [1,0,0]

    ax = f.add_subplot(111)

    ax.scatter(xvals,yvals,array(ugfracs) * 200,
               color = colors)
    ax.set_ylabel('mutation score')
    ax.set_xlabel('free energy (-kCal)')
    ax.annotate('''Evaluated structures positioned by energy 
and a mutation based evolutionary score.
Color indicates a second score from paired BL.
Radius indicates percentage of ungapped base pairs.''' , [0,1],
                xycoords = 'axes fraction',
                        xytext = [10,-10], textcoords='offset pixels',
                        va = 'top')
          
    myplots.padded_limits(ax, xvals, yvals, .2)
    f.savefig(figfile.format('{1}_frac_silent_{0}'.format(rname,rtype)))


    seq = structs['seq']
    n, selection_type = [4,'both']
    idxs = get_interesting_inds(xvals, yvals, structs, energies, pfracs,
                                rname, rtype, ns,rfid,figsize, colors,
                                seq,n,selection_type)

    if draw_single:
        show_rna_structs(xvals, yvals, structs, energies, pfracs,
                         rname, rtype, ns,rfid,figsize, colors,
                         seq,n,selection_type, idxs)
    if draw_many:
        for n, selection_type in \
                [[5,'ptime'],[5,'energy'],[ns,'energy']]:
            
            m_idxs = get_interesting_inds(xvals, yvals, structs, energies, pfracs,
                                rname, rtype, ns,rfid,figsize, colors,
                                seq,n,selection_type)
            show_rna_structs(xvals, yvals, structs, energies, pfracs,
                                    rname, rtype, ns,rfid,figsize, colors,
                                    seq,n,selection_type, m_idxs)
        
    return idxs,focus_tree

def get_interesting_inds(xvals, yvals, structs, energies, pfracs,
                         rname, rtype,ns,rfid,figsize,colors,
                         seq, n, selection_type):

      if selection_type == 'energy':
          vert_idxs = argsort(energies)[::-1][:n][::-1]
      elif selection_type == 'ptime':
          vert_idxs  =  argsort(pfracs)[::-1][:n][::-1]
      elif selection_type == 'both':
          vert_idxs1  =  argsort(pfracs)[::-1][:n][::-1]
          vert_idxs2  =  argsort(energies)[::-1][:n][::-1]
          vert_idxs = list(set(vert_idxs1).union(set(vert_idxs2)))
      return vert_idxs
def show_rna_structs(xvals, yvals, structs, energies, pfracs,
                     rname, rtype,ns,rfid,figsize,colors,
                     seq, n, selection_type,vert_idxs):


      verts = rutils.struct_verts([structs['structs'][i] for i in vert_idxs] 
                                  ,seq,rfid)
      
      f = myplots.fignum(3,figsize)
      ax = f.add_subplot(111)
      myplots.padded_limits(ax, xvals, yvals, .2)
      
      for vi, v in enumerate(verts):
            

            i = vert_idxs[vi]
            dims = [30]
            shadow_width = 10
            pkw0 = {'linewidth':shadow_width,
                    'color':'white',
                    'alpha':1,
                    'zorder':1.1}
            rplots.show_rna([xvals[i],yvals[i]], v,
                            dims = dims,
                            pkw = pkw0)
            
            pkw0 = {'linewidth':shadow_width,
                    'color':'white',
                    'alpha':.8,
                    'zorder':vi+2}
            rplots.show_rna([xvals[i],yvals[i]], v,
                            dims = dims,
                            pkw = pkw0)
            
            
            pkw1 = {'linewidth':2,
                    'color':colors[i],
                    'zorder':vi+2}
            rplots.show_rna([xvals[i],yvals[i]], v,
                            dims = dims, pkw = pkw1)
            ax.set_ylabel('mutation score')
            ax.set_xlabel('free energy (-kCal)')
            ax.annotate('''Suboptimal foldings, positioned by energy and
a mutation based evolutionary score.
Color indicates a second score from paired BL.''' , [0,1],xycoords ='axes fraction',
                        xytext = [10,-10], textcoords='offset pixels',
                        va = 'top')
          
      
      f.savefig(figfile.format('{3}_frac_silent_{0}_{1}{2}'.\
                                   format(rname,selection_type,n,rtype)))
      return vert_idxs



'''
NOTES:

SAM_alpha appears to have one --very-- highly conserved structure
that is strongly suboptimal.

'''
