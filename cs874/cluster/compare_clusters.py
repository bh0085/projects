import scipy.io as sio
import matplotlib.pyplot as plt
import compbio.utils.colors as mc
from numpy import *
import numpy as np
from scipy.stats import hypergeom

def test_constraints():
  prob2 = sio.loadmat('prob2.mat')
  
  domain_names = ['X', 'XV', 'Y', 'Y8', 'Y12']
  
  domains = [prob2.get(d) for d in domain_names]
  #domain_clusters = [prob2.get('ids_' + d) for d in domain_names]
  tissue_clusters = prob2.get('tissue_category')
  tc_inds = tissue_clusters.flatten() -1
  
  all_clusters = prob2.get('all_ids')
  all_fracs = []
  for cslice in all_clusters.T :
    c_inds = array(cslice).flatten()
    cpairs = set(['{0:d}x{1:d}'.format(ix,iy) 
                for ix, x in enumerate(c_inds) for iy, y in enumerate(c_inds)
                if ix < iy and x == y ])
    tcpairs = set(['{0:d}x{1:d}'.format(ix,iy) 
                for ix, x in enumerate(tc_inds) for iy, y in enumerate(tc_inds)
                if ix < iy and x == y ])
  
    max_pairs =( len(tc_inds) * len(tc_inds)  - len(tc_inds)) / 2
    total_pairs = len(cpairs.union(tcpairs))
    shared_pairs =len(cpairs.intersection(tcpairs))

    print 'found'
    print ' max pairs: {0}'.format(max_pairs)
    print ' total pairs: {0}'.format(total_pairs)
    print ' tissue pairs: {0}'.format(len(tcpairs))
    print ' cluster pairs: {0}'.format(len(cpairs))
    print ' shared pairs: {0}'.format(shared_pairs)

    all_fracs.append(float(shared_pairs)/ len(cpairs))
    

  lens = prob2.get('lens')
  f = plt.figure(5)
  f.clear()
  ax = f.add_subplot(111, title = 'fraction of pairs from expression clusters falling into tissue clusters',
                     xlabel = 'constraints (log base 10)', ylabel = 'frac correct')
  ax.plot(log10(lens),all_fracs )

  best_i = argmax(all_fracs)
  ax.scatter(log10(lens[best_i]), all_fracs[best_i], 300, color = 'black', facecolor = 'none', linewidth = 3)
  ax.annotate('Best clustering by agreement on pairs with tissues.\n nConstraints = {0:3g}'.format(float(lens[best_i][0])),\
                (log10(lens[best_i]),
                 all_fracs[best_i]), xytext = [15,0], textcoords = 'offset pixels', va  = 'top'),
                                                                                                                            
  f.savefig('figs/constraint_performance.tiff',format = 'tiff')
  return all_fracs


def run( domain_name = 'X', projection_name = 'Y8'  ):
  prob2 = sio.loadmat('prob2.mat')
  
  domain_names = ['X', 'XV', 'Y', 'Y8', 'Y12']
  
  domains = [prob2.get(d) for d in domain_names]
  #domain_clusters = [prob2.get('ids_' + d) for d in domain_names]
  tissue_clusters = prob2.get('tissue_category')
  

  clusters = domain_clusters[domain_names.index(domain_name)]
  pdom = domains[domain_names.index(projection_name)]
  cdom = domains[domain_names.index(domain_name)]

  f = plt.figure(1)
  f.clear() 
  random.seed(1)
  ct = array(mc.getct(218))
  
  #px, py = 2, 2
  sstrings = ['21{0:d}'.format(i+1) for i in range(4)]
  
  inds = arange(shape(dom)[1])
  
  c_inds = array(clusters).flatten() -1
  tc_inds = tissue_clusters.flatten() -1

  colors = ct[c_inds,:]

  ax = f.add_subplot(sstrings[0], title = \
                       'Clusters from genespace affinity. Projection to first two elements')  
  ax.scatter(*cdom[inds,0:2].T,s= 100, c = colors)
  ax = f.add_subplot(sstrings[1], title = \
                     'Clusters from genespace affinity. Projection to MVE')  
  ax.scatter(*pdom[inds,0:2].T,s= 100, c = colors)

  cpairs = set(['{0:d}x{1:d}'.format(ix,iy) 
                for ix, x in enumerate(c_inds) for iy, y in enumerate(c_inds)
                if ix < iy and x == y ])
  tcpairs = set(['{0:d}x{1:d}'.format(ix,iy) 
                for ix, x in enumerate(tc_inds) for iy, y in enumerate(tc_inds)
                if ix < iy and x == y ])
  f.savefig('figs/cluster_projectsions.tiff',format = 'tiff')
  
  max_pairs =( len(tc_inds) * len(tc_inds)  - len(tc_inds)) / 2
  total_pairs = len(cpairs.union(tcpairs))
  shared_pairs =len(cpairs.intersection(tcpairs))

  print 'using affinity propagation with affinites over domain {0}'.format(domain_name)
  print 'found'
  print ' max pairs: {0}'.format(max_pairs)
  print ' total pairs: {0}'.format(total_pairs)
  print ' tissue pairs: {0}'.format(len(tcpairs))
  print ' cluster pairs: {0}'.format(len(cpairs))
  print ' shared pairs: {0}'.format(shared_pairs)

  hg =  hypergeom( len(tcpairs), len(cpairs), max_pairs )
  return hg

