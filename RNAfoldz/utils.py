import subprocess as spc, sys, re, os,inspect, itertools as it, pickle
import matplotlib.mlab as mlab
import Bio.AlignIO as aio
import Bio
import compbio.utils.colors as mycolors
import compbio.utils.seismic as seismic
import compbio.utils.plots as myplots
import compbio.config as cfg


import rfam
import hcluster
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import myxml.parser as xparse
import myxml.svg as svg

import plots as rplots


def switch_dicts():
    switch_families = '''RF00050	FMN	Cis-reg; riboswitch
RF00059	TPP	Cis-reg; riboswitch
RF00162	SAM	Cis-reg; riboswitch
RF00167	Purine	Cis-reg; riboswitch
RF00168	Lysine	Cis-reg; riboswitch
RF00174	Cobalamin	Cis-reg; riboswitch
RF00234	glmS	Cis-reg; riboswitch
RF00504	Glycine	Cis-reg; riboswitch
RF00521	SAM_alpha	Cis-reg; riboswitch
RF00522	PreQ1	Cis-reg; riboswitch
RF00634	SAM-IV	Cis-reg; riboswitch
RF01054	preQ1-II	Cis-reg; riboswitch
RF01055	MOCO_RNA_motif	Cis-reg; riboswitch
RF01056	Mg_sensor	Cis-reg; riboswitch
RF01057	SAH_riboswitch	Cis-reg; riboswitch
RF01480	rli52	Cis-reg; riboswitch
RF01481	rli53	Cis-reg; riboswitch
RF01482	rli55	Cis-reg; riboswitch
RF01483	rli56	Cis-reg; riboswitch
RF01485	rli61	Cis-reg; riboswitch
RF01486	rli62	Cis-reg; riboswitch
RF01491	rli54	Cis-reg; riboswitch'''.split('\n')

    switch_dicts = dict([(fam.split('\t')[1],int(fam.split('\t')[0][2:]))
		     for fam in switch_families])
    return switch_dicts
    




def rna_draw(seq, struct, name, out_type = 'svg'):
    lines = '{0}\n{1}\n'.format(seq,struct)
    if out_type == 'png':
        outfile = cfg.dataPath('rnafold/{0}.png'.format(name))
        rprc = spc.Popen('RNAplot -o svg; convert rna.svg {0}'.format(outfile), shell = True,
                         stdin = spc.PIPE, stdout = spc.PIPE)
        
        out = rprc.communicate(input = lines)[0].splitlines()
  
        from matplotlib._png import read_png
        image = read_png(outfile)
    elif out_type== 'svg':
        outfile = cfg.dataPath('rnafold/{0}.svg'.format(name))
        
        rprc = spc.Popen('RNAplot -o svg; mv rna.svg {0}'.format(outfile), shell = True,
                         stdin = spc.PIPE, stdout = spc.PIPE)
        
        out = rprc.communicate(input = lines)[0].splitlines()
        
        
        struct_svg =  open(outfile).read()
        data = xparse.parse(struct_svg)
        arr = svg.get_polys(data)[0]

    else:
        raise Exception()
    
    return arr
def project_lstruct(pstruct, l):
	out = zeros(l)
	for p in pstruct:
		out[p[0]] = -1.
		out[p[1]] = 1.
	return out

def ungapped_seq(seq, seq_id, name = 'unspecified name' ):
    ungapped_ref =  Bio.SeqRecord.SeqRecord( \
        Bio.Seq.Seq(''.join(\
                [let for let in str(seq.seq) 
                 if let.upper() in 'AUCGT']),
                    seq.seq.alphabet)
        ,seq_id, name = name)
    return ungapped_ref
    
def pairs_stk(struct, l):
    stk = ['.'] * l
    for p in struct: stk[p[0]] ,stk[p[1]] = '(',')'
    stk = ''.join(stk)
    return stk

def stk_pairs(struct):
    #FETCH PAIRS IN STRUCT
   pairs = []
   pqueue = []
   for i, cstr in enumerate(struct):
       if cstr in '[<({':
           pqueue.append(i)
       elif cstr in ']>)}':
           pairs.append((pqueue.pop(),i))    
   assert not pqueue   
   return pairs


def stk_parse(fopen):
    ali = aio.parse(fopen, 'stockholm').next()
    fopen.seek(0)
    str_re = re.compile('^#=GC SS_cons\s*(.*)$')
    ref_re = re.compile('^#=GC RF\s*(.*)$')

    struct, ref = '', ''
    for l in fopen.xreadlines():
        match = str_re.search(l)
        if match: 
            struct += match.group(1)

        match = ref_re.match(l)
        if match: 
            ref += match.group(1)
    
    return ali, ref, struct
            
def stk_format(seq, struct):

    return '''# STOCKHOLM 1.0

{0:20} {2}
{1:20} {3}
//
'''.format(seq.id,'#=GC SS_cons', seq.seq, struct)




def stk_mat(stk):
    mat = zeros( (len(stk), 3))
    map_dict = {'(':2, ')':0,'.':1}
    mark = [ ( i, map_dict[elt] ) for i, elt in enumerate(stk)]
    mat[zip(*mark)] = 1
    return mat

def struct_project_l(pairs, l):
    '''
    Project a structure onto a vector of length l.
    Values are -1, 0, 1 depending on whether a pair exists, opens or closes.
    
    Returns projection, an "l-"vector
'''
    return project_lstruct(pairs, l)
def struct_project_l2(pairs,l):
    '''
    Project a structure onto an l-squared-vector.
    Every pair puts two ones in the vector: one at i+j*l, the other at i*l + j
'''
    p = array(list(pairs))
    #import scipy.sparse.lil as lil
    z = zeros(l**2)
    z[sum(p * array([1,l]),1)] = 1 
    z[sum(p * array([l,1]),1)] = 1 
    return z
def struct_project_l3(pairs,l):
    return stk_mat(pairs_stk(pairs,l))


def cluster_2(structs,counts, seq,
              n_comp = 10,
              ptype = 'pca' ):


    

    if len(seq) < 300:
        vecs = project_structs(structs, l = len(seq) ,ptype = 'pairs_full')
    else:
       #COMPUTE PROJECTIONS NECESSARY FOR THE CLUSTERING.
        affinities = struct_affinity_matrix(structs,len(seq), aff_type = 'pairs', ss_multiplier = None)
        vecs = project_structs(structs, ptype = 'pca',
                               affinities = affinities,
                               n_comp = min([len(affinities),50]))
    

    #CLUSTER STRUCT PROJECTIONS
    clusters = cluster_structs(vecs, k =30)
    
    return clusters

   
    
def show_subopts(structs, polys, energies):
    srted = argsort(energies)
    e = array(energies)
    cols = [1.,0.,0.] * ((e - min(e)) / (max(e) - min(e)))[:,newaxis]
    plf2 = myplots.fignum(7,(10,10))

    rplots.grid_rnas(polys[srted], 
                     colors =cols[srted],
                     size = (8,8), dims = [180,50])    
    raise Exception()


def cluster_2_show(clusters, polys): 
    sortorder = argsort(clusters)
    ct_colors = mycolors.getct(len(set(clusters)))
    ct_dict = dict([(cluster, ct_colors[i]) for i, cluster in enumerate(set(clusters))])
    
    plf2 = myplots.fignum(8,(10,10))
    
    rplots.grid_rnas(polys[sortorder], 
                     colors  = [ct_dict[i] for i in clusters[sortorder]],
                     size = (5,5), dims = [180,50])    

def project_structs(structs,
                    ptype ='l',
                    affinities = None,
                    n_comp = None,
                    l = None,
                    vecs = None):
    '''
Project RNA structures in any one of several ways.
Different projections require different inputs.

inputs:
 ptype:  ['pca', 'rnd', 'l', 'full_pairs']
 affinities: aff matrix for pca
 n_comps:    n_comps for pca
 l:          length for l projections and pairs projections
 vecs:       vecs from l/pairs projection for random projections


outputs:
 projections in the form of an [N, X] matrix where X is the 
 size of the projection and N is the number of input structures.

The projections requiring the fewest input variables
are 'full_l' and 'full_pairs' as these only require a list
of structures (specifiied as base pairs) and a sequence length.
Most of the rest can be called using the output projections
from 'full_l' or 'full_pairs' as input.

In particular, we can project onto PCA vectors:
  pca inputs: 
              n_comps: the number of components to take from 
                       PCA projection
              affinities: the affinity matrix to use for the
                          projection
  
Or random matrices:
  rnd inputs:
              vecs: vectors in l dimensional space to project
                    onto random matrices.

'''

    if ptype == 'pca':
        assert affinities != None
        assert n_comp != None
        pca_vecs = mlab.PCA(affinities).project(affinities)  
        pca_vecs = pca_vecs[:,0:n_comp]
        return pca_vecs
    elif ptype == 'l':
        assert l != None
        return array([struct_project_l(p, l) for p in structs])
    elif ptype == 'full_pairs':
        assert l != None
        return [struct_project_l2(p, l) for p in structs]
    elif ptype == 'rnd':
        assert n_comp != None
        assert vecs != None
        mat = array(np.round(random.rand(n_comp, l)),float)
        mat *= 2
        mat -= 1
        mat/= sqrt(l)
        plt.imshow(mat)
        cvecs = dot(mat,vecs.T).T
        return cvecs


    else:
        raise Exception('Projection type: {0} not yet implemented'.format(ptype)) 


def cluster_struct_affinities(affinities, 
                              ctype = 'mlpy',
                              k = None):
    ''' Return a list of cluster memberships in the form of an N-array having 
k unique elements.
'''
    if ctype == 'hcluster':
        return hcluster.fcluster(vecs,1.1,criterion='inconsistent',method = 'complete' )

    elif ctype == 'mlpy':
        import mlpy
        HC = mlpy.HCluster(method='euclidean', link='complete')
        clusts = HC.compute(vecs )
        cut = HC.cut(HC.heights[-k])
        return cut
    else: 
        raise Exception()

def cluster_structs(vecs, 
                    ctype = 'mlpy',
                    k = None):
    ''' Return a list of cluster memberships in the form of an N-array having 
k unique elements.
'''
    if ctype == 'hcluster':
        return hcluster.fclusterdata(vecs,1.1,criterion='inconsistent',method = 'complete' )

    elif ctype == 'mlpy':
        import mlpy
        HC = mlpy.HCluster(method='euclidean', link='complete')
        clusts = HC.compute(vecs )
        cut = HC.cut(HC.heights[-k])
        return cut
    else: 
        raise Exception()

def old_clusters():
    
    plf = myplots.fignum(6, (8,8))
    plf.clear()
    ax = plf.add_subplot(211)
    
    if do_rnd:
        hstart = .15
    else:
        hstart = 3.

    all_vars = []
    all_vars_n = []
    all_clusters = []
    ks = []

    all_Bvars, all_Wvars = [], []
    theights = log(HC.heights[greater(HC.heights,hstart)][-len(cvecs)/2:][:-1])

    for xval, h in enumerate(theights):
        clustering = HC.cut(exp(h))
        casrt = argsort(clustering)
        csrtd = clustering[casrt]
        d = dict([(k,array(list(g))) for k, g in it.groupby(zip(casrt,csrtd), 
                                            key = lambda x: x[1])])
        lens = array([len(v) for v in d.values()],float)
        nlens = lens / max(lens)
        
        cmeans = array([mean(cvecs[idxs[:,0],:],0)
                        for idxs in d.values()])
        Wvars = np.sum(array([np.mean(  (cvecs[idxs[:,0],:] - cmeans[i])  **2 )
                              for i, idxs in enumerate(d.values())]))

        
    
        Bvars = np.sum(lens[:,newaxis]* ( (cmeans - mean(cvecs,0)) **2)  )
    
        ks.append(len(d))

        all_Wvars.append(Wvars)
        all_Bvars.append(Bvars)
        

        cluster_vars =array([ sum(var(cvecs[idxs[:,0],:],0)) 
                              for idxs in d.values() ])

        
        cluster_vars_n =array([ sum(var(cvecs[idxs[:,0],:],0)) 
                              for idxs in d.values() ])/(lens)

        all_clusters.append([cvecs[idxs[:,0]]
                             for idxs in d.values() ])
        all_vars.append(cluster_vars)
        all_vars_n.append(cluster_vars_n)
        colors = array(argsort(argsort(lens)),float)/len(lens)
        
        ax.scatter(0*(cluster_vars) +  h , cluster_vars_n, 20,
                   color = array(array([0.,1.,0.]) * colors[:,newaxis]))

    ax2 = plf.add_subplot(212)
    all_Bvars = array(all_Bvars)
    all_Wvars = array(all_Wvars)
    
    density_based = False
    HC_based = True
    if density_based:
      #ax3 = plf.add_subplot(212,frameon = False)
      
      #COMPUTE A COALESCENCE RATE OF CLUSTERS
      #(divide the pde for heights by clustering size)
      from scipy.stats import gaussian_kde
      data = (theights)
      
      density = gaussian_kde(data)
      density_rate = 1. / array([len(v) for v in all_vars])
      xs = (theights)
      
      density.covariance_factor = lambda : .25
      density._compute_covariance()
      yvals , colors = array([density(xs), 
                              density(xs) * density_rate,
                              [sum(v/[len(c) for c in cs] )
                               for v,cs in zip(all_vars,all_clusters) ]]),\
                              [[1,0,0],[0,1,0],[0,0,1]]
      yvals /= np.max(yvals,1)[:,newaxis]
      xvals = (xs) + zeros(len(yvals))[:,newaxis]
      
      #for i in range(len(yvals)):
      #    ax2.plot(xvals[i],yvals[i], color = colors[i])
      
      dens_n = yvals[1]
      #dens_n[greater(dens_n,percentile(dens_n,60))] = percentile(dens_n,60)
      #dens_n[less(dens_n, percentile(dens_n,.25))] = percentile(dens_n,25)
      diff = dens_n - yvals[2]
      ax2.plot(xvals[0], yvals[1] - yvals[2], linewidth = 10)
      
      import scipy.signal as ss
      #diff =ss.medfilt(diff,3)
      ax2.plot(xvals[0], diff, linewidth = 5, color = 'orange')
      
      mpt = argmin(diff)
      #m#pt = len(yvals[0]) - 3
      
      ax2.scatter([xvals[0][mpt]],diff[mpt], 200, color = 'red')
      
      ax2.plot(xvals[0],yvals[1])
      
      #ax2.plot(histogram(log(HC.heights[greater(HC.heights,hstart)]))[1][:-1],
      #         histogram(log(HC.heights[greater(HC.heights,hstart)]))[0])
      #raise Exception()
    else:
        ax2.plot(ks,all_Bvars)
        ax2.plot(ks,all_Wvars)
        
        hcfun =  array(all_Bvars) / array(all_Wvars) /\
            ((array(ks,float) -1) / array(float(len(cvecs)) - array(ks,float)))

        
        hcfun = nan_to_num(hcfun)
        
        mpt = argmax(hcfun)

        
        a3 = plf.add_subplot(212, frameon = False)
        a3.plot(ks, hcfun)

        ax2.scatter([ks[mpt]]*2, [all_Bvars[mpt],all_Wvars[mpt]], 200, color = 'orange')
        
 
def suboptimals(sequence, sp_method = 'enumerate', n = 10000, name = 'NONAME'):
        fa_str = sequence.format('fasta')
        if sp_method == 'enumerate':            
            rprc = spc.Popen('RNAsubopt --deltaEnergy=2.5 ', shell = True,
                             stdin = spc.PIPE, stdout = spc.PIPE)
        elif sp_method == 'sample':
            rprc = spc.Popen('RNAsubopt --stochBT={0}'.format(n), shell = True,
                             stdin = spc.PIPE, stdout = spc.PIPE)

        out = rprc.communicate(input = fa_str)[0].splitlines()
        print 'print computing rna structures for method {0}'.format(sp_method)
        struct_pairs = [set(stk_pairs(p))
                  for p in out[3:]][::]
        return struct_pairs

def struct_energy(sequence, struct):
    strs = '''{0}
{1}'''.format(sequence.seq,pairs_stk(struct, len(sequence)))
    
    prc = spc.Popen('RNAeval ', shell = True,
                    stdin = spc.PIPE, stdout = spc.PIPE)     
    out = prc.communicate(input = strs)[0]
    energy =float(re.compile('[\d]+\.*[\d]*').search(out).group())
    return energy

def subtree_refseq(subtree, method = 'root'):
    '''
Get a reference sequence for a subtree of aligned RNA sequences.

inputs
  subtree: a biopython Tree having seqs in node.m['seq']

outputs
  node:    a biopython Clade instance in subtree
  seq:     a biopython SequenceRecord instance

'''
    
    #FIND A REFERENCE SEQUENCE (CLOSEST TO THE TREE ROOT)
    if method == 'root':
        tset = set(subtree.get_terminals())
        t_depths = [ (k,v) for k,v in subtree.depths().iteritems() if k in tset ]
        node = t_depths[argmin([t[1] for t in t_depths])][0]
        seq = node.m['seq']
    elif method == 'median':
        raise Exception('median method not yet implemented')

    return node, seq


def family_clustered_suboptimals(rfid, plots = True, num = 1000, min_count = 2,
                                 n_countsorted = 10, n_esorted = 10, 
                                 draw = False, cluster_type = 'just_list',
                                 savename = None):
    if savename == None:
        savename = rfid
    ali, tree, infos = rfam.get_fam(rfid)
    ali_ids = [a.name for a in ali]

    for i, n in enumerate(tree.get_terminals()):
        match = re.compile('_([^_]*)_').search(n.name) 
        if not match or not '/' in match.group(1):
            this_seq = []
        else:
            term_id = match.group(1)
            this_seq = ali[ali_ids.index(term_id)]
        n.m = {'seq':this_seq,
               'probs':[1 for j in range(len(this_seq))]}

    big_refnode, big_refseq = \
        subtree_refseq(tree)
    ungapped_ref = ungapped_seq(big_refseq, rfid)
    seq = ungapped_ref
    structs = suboptimals(ungapped_ref, sp_method = 'sample',name = rfid, n = num)

    stks = [pairs_stk(s,len(seq)) for s in structs]
    stk_srt = sorted([ (i,s) for i,s in enumerate(stks)], key = lambda x: x[1])
    stk_groups = [ list(g) for k, g in it.groupby(stk_srt,key =lambda x: x[1])]
    stk_unq, struct_counts = zip(*[( g[0][0] , len(g))  for g in stk_groups])
    structs  = [structs[elt] for elt in stk_unq ]
 
   
    if cluster_type == 'full_clustering':
        final_structs, final_energies = select_exemplars_from_clustering(structs,struct_counts,seq, draw = draw)
    elif cluster_type == 'just_list':
        final_structs, final_energies = select_exemplars_from_list(structs,struct_counts,seq, draw = draw)

    if draw:
        print 'DRAWING final subopts' 
        verts = struct_verts(final_structs, seq)
        show_subopts(final_structs, verts, final_energies)
        f = plt.gcf()
        f.savefig(cfg.dataPath('figs/RNAfoldz/exemplars_{0}.ps'.format(savename)))
    
    fopen = open(cfg.dataPath('RNAfoldz/subopts_{0}.pickle'.format(savename)),'w')

    return final_structs,final_energies, seq
    pickle.dump({'structs':final_structs, 'energies':final_energies, 'seq':seq},
                fopen)
    

def select_exemplars_from_clustering(structs,struct_counts,seq, draw = False):
      min_count = 2

      freq_structs = [s for i, s in enumerate(structs) if struct_counts[i] >= min_count]
      if len(freq_structs) < 10:
          min_count = 1
          freq_structs = [s for i, s in enumerate(structs) if struct_counts[i] >= min_count]
      struct_counts= [s for i, s in enumerate(struct_counts) if s >= min_count]
      structs = freq_structs
      
      struct_energies = [struct_energy(seq, s) for s in structs]
      if len(structs) > 400:
          high_e = argsort(struct_energies)[::-1][:400]
          structs =[ structs[i] for  i in high_e]
          struct_counts =[ struct_counts[i] for  i in high_e]
          struct_energies =[ struct_energies[i] for  i in high_e]
                     
          
      
      clusters = cluster_2(structs,  struct_counts, seq, ptype = 'full_pairs')
      if draw:
          print 'DRAWING Clusters'
          verts = struct_verts(structs, seq)
          cluster_2_show(clusters, verts)
          f = plt.gcf()
          f.savefig(cfg.dataPath('figs/RNAfoldz/clusters_{0}.ps'.format(savename)))
      exemplars = set(clusters)
      cluster_exemplars = []
      for e in exemplars:
          reps =array([ (i, eng) for i, eng in enumerate(struct_energies) if clusters[i] == e])
          min_rep = reps[:,0][argmax(reps[:,1])]
          cluster_exemplars.append(min_rep)
          
      cluster_exemplars = set([int(e) for e in cluster_exemplars])
      sorted_exemplars = set(argsort(struct_counts)[::-1][:n_countsorted])
      energy_exemplars = set(argsort(struct_energies)[::-1][:n_esorted])
      final_exemplars = cluster_exemplars.union(sorted_exemplars).union(energy_exemplars)

      print '''Structural exemplars found:
Clustering:     {0}  {4}
Count sorting:  {1}  {5}
Energy sorting: {2}  {6}

Total unique:   {3}'''.format(len(cluster_exemplars),len(sorted_exemplars), 
                              len(energy_exemplars),len(final_exemplars),
                              mean([struct_energies[i] for i in cluster_exemplars]),
                              mean([struct_energies[i] for i in sorted_exemplars]),
                              mean([struct_energies[i] for i in energy_exemplars]))


      final_structs, final_energies = zip(*[(structs[i],struct_energies[i]) for i in  final_exemplars])
      return final_structs, final_energies

def select_exemplars_from_list(structs, struct_counts,seq, draw = False):
      min_count = 1
      max_structs = 2000
      jacc_sim_thr = .80

      freq_structs = [s for i, s in enumerate(structs) if struct_counts[i] >= min_count]
      if len(freq_structs) < 10:
          min_count = 1
          freq_structs = [s for i, s in enumerate(structs) if struct_counts[i] >= min_count]
      struct_counts= [s for i, s in enumerate(struct_counts) if s >= min_count]
      structs = freq_structs
      
      struct_energies = [struct_energy(seq, s) for s in structs]
      if len(structs) > max_structs:
          high_e = argsort(struct_energies)[::-1][:max_structs]
          structs =[ structs[i] for  i in high_e]
          struct_counts =[ struct_counts[i] for  i in high_e]
          struct_energies =[ struct_energies[i] for  i in high_e]

      esrt = argsort(struct_energies)
      affinities = struct_affinity_matrix(structs, len(seq), aff_type = 'jaccard')
      affinities = affinities[:,esrt][esrt]
      inds_taken = array((0,), int)
      for i in range(len(affinities)):
          this_aff = np.max(affinities[i,inds_taken]) if len(inds_taken) > 0 else 0
          if this_aff < jacc_sim_thr:
              inds_taken = append(inds_taken, i)
      final_exemplars = [ esrt[i] for i in inds_taken]

      final_structs, final_energies = zip(*[(structs[i],struct_energies[i]) for i in  final_exemplars])
      return final_structs, final_energies


def struct_verts(structs, seq):
    verts = array([rna_draw(seq.seq , pairs_stk(sp,len(seq)), 'name' )
             for sp in structs])
    return verts



def struct_affinity_matrix(pairs,seq_len, aff_type = 'pairs', ss_multiplier = None):
        '''compute an affinity matrix for N RNA structures specified as 
        a list of base pairs'''


	#AFFINITY TYPE
	if aff_type == 'easy':
            #SHAPE COMPARISON
            ssm = .01 if ss_multiplier == None else ss_multiplier
            pair_vecs = zeros((len(pairs),seq_len))
            for i, spair in enumerate(pairs):
                for p_elt in spair:
                    pair_vecs[i][p_elt[0]] = 1 
                    pair_vecs[i][p_elt[1]] = -1

            nrm =  sqrt( sum(pair_vecs**2,1)[:,newaxis] )
            if sum(equal(nrm,0)) > 0: raise Exception()
            nrm[equal(nrm,0)] == 1
            
            pair_vecs /= nrm
            affinities = sum(pair_vecs[:,newaxis,:] * pair_vecs[newaxis,:,:],2)
            da = np.max(affinities) - np.min(affinities)
            ss = np.min(affinities)  + da *ssm
        
        elif aff_type == 'pairs':
            #PAIR INTERSECTION COMPARISON
	    ssm = -.1 if ss_multiplier == None else ss_multiplier
            affinities = array([[float(len(pairs[i].intersection(pairs[j])))\
					 /  sqrt((len(pairs[i])*len(pairs[j])))
                                 for i in range(len(pairs))]
                                for j in range(len(pairs))])
            da = np.max(affinities) - np.min(affinities)
            ss = np.min(affinities) + da *ssm
        elif aff_type == 'jaccard':
            affinities = array([[float(len(pairs[i].intersection(pairs[j])))\
                                     / len(pairs[i].union(pairs[j]))
                                 for i in range(len(pairs))]
                                for j in range(len(pairs))])
        else:
            raise Exception()
	return affinities
		       
