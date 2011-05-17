import utils as rutils
import plots as rplots
import score_utils as sutils
import compbio.utils.bsub_utils as bsu

from numpy import *
import numpy as np, itertools as it

import matplotlib.pyplot as plt

import compbio.utils.plots as myplots
import compbio.utils.colors as mycolors
import compbio.utils.memo as mem
import compbio.config as cfg

import pickle

figsize = (8,8)
figtype = 'ps'
figfile = cfg.dataPath('figs/gpm2/pt3_fana/{{0}}.{0}'.format(figtype))

do_make_figs = True
if do_make_figs:
    do_make_subopts = False

#flist = [50,311,140,143,495,637,1304]
flist = [311,1304]

def setFamData(rfid = None, ftype = None,**kwargs):

    assert rfid; assert ftype;
    fprefix = 'FA' if ftype == 'all' else 'RS'
    sdat = bsu.load_data('{1}_{0}'.format(rfid,fprefix), 'output')
    tdat = bsu.load_data('{1}_tree_{0}'.format(rfid,fprefix), 'output')

    return sdat, tdat

def setTree(rfid = None, ftype = None,**kwargs):
    assert setTree; assert ftype
    return sutils.get_tree(rfid, fam_type = ftype)

def show_conservation(fidx = 0, reset = False):
    fnum = flist[fidx]
    rfid = 'RF{0:05}'.format(fnum)
    print rfid
    if fnum ==50: ftype = 'riboswitch'
    else: ftype = 'all'
    
    
    out = mem.getOrSet(setFamData,
                              **mem.rc({}, reset =reset,
                                       on_fail = 'compute',
                                       hardcopy = False,
                                       register = 'fdat'+rfid,
                                       ftype = ftype,
                                       rfid = rfid))

    
    mvals, tvals, structs = mem.getOrSet(setTree,
                                         **mem.rc({},reset = reset,
                                                  on_fail = 'compute',
                                                  hardcopy = True,
                                                  register = 'st'+rfid,
                                                  rfid = rfid,
                                                  ftype = ftype))
    
    idxs, tidx  = sutils.show_paired_v_energy(rfid,rfid,mvals,tvals,structs,ftype)
    
    all_pairs = structs['structs']
    all_energies = structs['energies']
    
    pints,eints, mints, tints = [structs['structs'][i] for i in idxs],\
        [ structs['energies'][i] for i in idxs],\
        [ mvals[tidx][i] for i in idxs],\
        [ tvals[tidx][i] for i in idxs]
    seq = structs['seq']
    
    if do_make_subopts:
        subopts = rutils.suboptimals(seq, n = 400)
        verts = rutils.struct_verts(subopts, seq, rfid)
        f = myplots.fignum(4,figsize)
        rplots.grid_rnas(verts, dims = [40])
        f.savefig(figfile.format('{0}_grid_rnas'.\
                                     format(rfid)))

    
                



    aff = rutils.struct_affinity_matrix(all_pairs, len(seq))
    pca = rutils.project_structs(all_pairs,
                          ptype ='pca',
                          affinities = aff,
                          n_comp = 3) 

    for metric in ['n_comp']:# ['frac_silent','frac_paired','n_comp']:
      scolors = []
      for i in range(len(tvals[tidx])):
          m_silent, pidxs, frac_good = sutils.metric(
              mvals[tidx][i],tvals[tidx][i],
              mtype = metric)
          
          scolors.append(mean(m_silent))
      scolors = myplots.rescale(scolors, [0.,1.])[:,newaxis] * array([1.,0.,0.])
      
      
      f = myplots.fignum(4,figsize)
      ax = f.add_subplot(111)
      xvals, yvals = pca[:,:2].T
      myplots.padded_limits(ax, xvals, yvals)
      
      ax.scatter(xvals,yvals,300,linewidth = 1,
                 edgecolor = 'black', color = scolors)

      ax.scatter(pca[idxs,0],pca[idxs,1], 2100 ,alpha = 1, 
                 color = 'black')
      ax.scatter(pca[idxs,0],pca[idxs,1], 2000 ,alpha = 1, 
                 color = 'white')
      ax.scatter(pca[idxs,0],pca[idxs,1], 400 ,alpha = 1, 
                 color = scolors[idxs],
                 )


      ax.annotate('''Conservation metric: {0}
Projected onto C=2 Principal Components'''.format(metric),
                  [0,1],xycoords = 'axes fraction', va = 'top',
                  xytext = [10,-10],textcoords='offset points')
      
      f.savefig(figfile.format('{0}_pca_{1}'.\
                                 format(rfid, metric)))



    #verts = rutils.struct_verts(sint, seq, rfid)
    
