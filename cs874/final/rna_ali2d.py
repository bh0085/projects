#!/usr/bin/env python
'''
Given an rfam family, its alignment and the tree given, compute ancestral nodes
and mutations along the edges.

'''
import Bio.AlignIO as aio
import Bio.Align as ba
import Bio.Phylo.NewickIO as nio
import Bio.Phylo as phy
import Bio
import subprocess as spc, sys, re, os,inspect

import compbio.utils.memo as mem
import compbio.utils.bs_macros as bsm
import compbio.config as cfg

if __name__ != '__main__' and False:
	import compbio.utils.colors as mycolors
	import compbio.utils.seismic as seismic
	import compbio.utils.plots as myplots
	import matplotlib.pyplot as plt

from numpy import *
import numpy as np
import itertools as it

from compbio.projects.seqtree import muscle, phyml, paml

import sys
def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()


def plot_clusters(inds,
                  embeddings,
                  plot3d = False,
                  title = '',
		  ax_in =None,
		  save = False,
		  colors = None):
        exemplars = list(set(inds))
        if colors == None:
		cluster_colors = dict([(exemplars[i], col) 
                              for i, col in enumerate(mycolors.getct(len(exemplars)))]
                              )

        cols = [cluster_colors[e] for e in inds]
        try: 
		if ax == None: plt.clf()
        except Exception, e: pass
        if ax_in == None: f = plt.gcf()

        for i, k in enumerate(embeddings.keys()):
            embedding = embeddings[k]

	    #if i == 1: raise Exception()
            emb_sig = embedding[:,0:3]
            cluster_vars = [ var(emb_sig[nonzero(equal(inds, j))[0]])  for j in exemplars]
            indexed_vars = [ cluster_vars[exemplars.index(j)] for j in inds ]
	    indexed_vars[equal(index_vars,0)] = 1

            sizes = 10 *( exp( -1 * ( np.sum((emb_sig - emb_sig[inds,:])**2,1)/indexed_vars)))
            if plot3d:
                if ax_in == None: 
			ax = f.add_subplot('{1}1{0}'.format(i+1, len(embeddings)),projection = '3d')
		else: ax = ax_in
                ax.scatter(array(embedding[:,0],float)
                           ,array(embedding[:,1],float)
                           ,array(embedding[:,2],float), 
                           s = sizes,
                           color = cols)
                ax.set_xticks([])
                ax.set_yticks([])
                for tl in list(it.chain( ax.w_xaxis.get_ticklabels(),
                                    ax.w_yaxis.get_ticklabels(),
                                    ax.w_zaxis.get_ticklabels())): # re-create what autofmt_xdate but with w_xaxis
                    tl.set_visible(False)
                    tl.set_rotation(30)    
            else:
                if ax_in == None: ax = f.add_subplot('{1}1{0}'.format(i+1, len(embeddings)))
		else: ax = ax_in
                ax.scatter(array(embedding[:,0],float)
                           ,array(embedding[:,1],float),
                           s = sizes,
                           color = cols)
            ax.set_title('{0} for subopts in {1}'.format(k, title))
        
        if save: 
		f.savefig(cfg.dataPath('cs874/figs/subopt_embeddings/{0}.ps').format(title))


def parse_stk_struct(filename):
    fopen = open(filename)
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
            
    
    
def get_fam(ofs = 0, rfid = None):

    fopen = open(cfg.dataPath('rfam/Rfam.seed'))
    #fopen.readline()
    #fopen.readline()
    alis = aio.parse(fopen,'stockholm')

    if rfid != None:
	    ofs = int(re.compile('\d+').search(outputs['title']).group()) -1
	    
    for i in range(ofs):
        null = alis.next()
    
    
    infos = {}
    start = fopen.tell()
    while 1:
        l = fopen.readline()            
        if l[0] == '#':
            ukey = str(l[5:7])
            infos.update( [(ukey, infos.get(ukey,'') + l[8:])])
            
        else:
            if l.strip() != '': break
            
    fopen.seek(start)

    ali = alis.next()
    rfname = infos['AC'].strip()
    fname = cfg.dataPath('rfam/Rfam.seed_tree/{0}.seed_tree'.format(rfname))

    tree = nio.parse(open(cfg.dataPath('rfam/Rfam.seed_tree/{0}.seed_tree'.format(rfname)))).next()
    return ali, tree, info, rfname

def infernal_alignment(alis, profile,rfid):
    seqs = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([let for let in str(ali.seq)  if let in 'AUTGC' ]),
                                                Bio.Seq.Alphabet.RNAAlphabet),
                                    'S{0:03}'.format(idx))
            for idx, ali in enumerate(alis)]
    infile = cfg.dataPath('infernal/temp/{0}_{1:03}_unaligned.fa'.format(rfid,idx))
    outfile= cfg.dataPath('infernal/temp/{0}_{1:03}_aligned.stk'.format(rfid,idx))
    Bio.SeqIO.write(seqs, infile, 'fasta')
    
    cstr = 'cmalign -o {0} {1} {2}'.format(outfile, profile, infile)
    ispc = spc.Popen(cstr, shell = True,
                            stdout = spc.PIPE)
    out = ispc.communicate()[0]
    alis, stk, struct = parse_stk_struct(outfile)
    ali = ba.MultipleSeqAlignment(alis)

    for a in ali:
	    a.seq = a.seq.upper()
    return ali, stk, struct
    
def infernal_profiles(ungapped_ref, exemplar_structs):
    exemplar_stks = []
    for i, s in enumerate(exemplar_structs):
        stk = ['.'] * len(ungapped_ref)
        for p in s: stk[p[0]] ,stk[p[1]] = '(',')'
        stk = ''.join(stk)
        exemplar_stks.append('''# STOCKHOLM 1.0

{0:20} {2}
{1:20} {3}
//
'''.format(ungapped_ref.id,'#=GC SS_cons', ungapped_ref.seq, stk))

    profiles = []
    for idx, stktext in enumerate(exemplar_stks):
        stkfile = cfg.dataPath('infernal/temp/{0}_{1:03}.stk'.format(ungapped_ref.id, idx))
        cmfile = cfg.dataPath('infernal/temp/{0}_{1:03}.cm'.format(ungapped_ref.id, idx))
        fopen = open(stkfile,'w')
        fopen.write(stktext)
        fopen.close()
        cstr = 'cmbuild -F --rsearch {0} {1} {2}'.format(cfg.dataPath('infernal/matrices/RIBOSUM85-60.mat'),cmfile,stkfile)
        ispc = spc.Popen(cstr, shell = True, stdout = spc.PIPE)
        out = ispc.communicate()[0]
        profiles.append(cmfile)
    
    return profiles
        


def compute_embedding(spairs,
		      aff_type = 'pairs',
		      do_mve = False,
		      ss_multiplier = None):

	#Artificial sequence length from the maximum paired element
	seq_len = np.max(list(it.chain(*spairs))) + 1

	#AFFINITY TYPE
	if aff_type == 'easy':
            #SHAPE COMPARISON
            ssm = .01 if ss_multiplier == None else ss_multiplier
            pair_vecs = zeros((len(spairs),seq_len))
            for i, spair in enumerate(spairs):
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
	    nrm = sqrt((len(spairs[i])*len(spairs[j])))
	    nrm[equal(nrm,0)] = 1
            affinities = array([[float(len(spairs[i].intersection(spairs[j])))\
                                     / nrm
                                 for i in range(len(spairs))]
                                for j in range(len(spairs))])
            da = np.max(affinities) - np.min(affinities)
            ss = np.min(affinities) + da *ssm
	mve_dict = dict(k = 4,
			similarities = affinities,
                        do_mve = False if aff_type == 'easy' else False)
        embedding=  bsm.runmat('mve_frompy', mve_dict, 'subopt_mve')
        mve_vecs = embedding['Y'].T
        pca_vecs = embedding['YPCA'].T

	return pca_vecs, mve_vecs

		       
def setAffinities(tree = None, refseq_method = None, refseq = None,
                      sp_method = 'sample',aff_type = 'pairs',
                      **kwargs):


        fa_str = refseq.format('fasta')
        if sp_method == 'enumerate':            
            rprc = spc.Popen('RNAsubopt --deltaEnergy=2.5 ', shell = True,
                             stdin = spc.PIPE, stdout = spc.PIPE)
        elif sp_method == 'sample':
            rprc = spc.Popen('RNAsubopt --stochBT=250', shell = True,
                             stdin = spc.PIPE, stdout = spc.PIPE)

        out = rprc.communicate(input = fa_str)[0].splitlines()
        print 'print computing rna structures for method {0}'.format(sp_method)
        spairs = [set(stk_pairs(p))
                  for p in out[3:]][::]

	
        if aff_type == 'easy':
            pair_vecs = zeros((len(spairs),len(out[-1])))
            for i, spair in enumerate(spairs):
                for p_elt in spair:
                    pair_vecs[i][p_elt[0]] = 1 
                    pair_vecs[i][p_elt[1]] = -1

            nrm =  sqrt( sum(pair_vecs**2,1)[:,newaxis] )
            if sum(equal(nrm,0)) > 0: raise Exception()
            nrm[equal(nrm,0)] == 1
            pair_vecs /= nrm
            affinities = sum(pair_vecs[:,newaxis,:] * pair_vecs[newaxis,:,:],2)
            da = np.max(affinities) - np.min(affinities)
            ss = np.min(affinities)  + da/100
        
        elif aff_type == 'pairs':
            affinities = array([[float(len(spairs[i].intersection(spairs[j])))\
                                     /sqrt((len(spairs[i])*len(spairs[j])))
                                 for i in range(len(spairs))]
                                for j in range(len(spairs))])
            da = np.max(affinities) - np.min(affinities)
            ss = np.min(affinities) - da/10
        


        cluster_dict = dict(similarities=affinities,
                           self_similarity = ss)
        clusters =  bsm.runmat('ap_frompy', cluster_dict, 'subopt_clusters')
        ex_inds = list(set(list(squeeze(clusters['inds'].T))))
        exemplars = affinities[ex_inds][:,ex_inds]

        mve_dict = dict(k = 4,
                        similarities = affinities,
                        do_mve = False if aff_type == 'easy' else False)
        embedding=  bsm.runmat('mve_frompy', mve_dict, 'subopt_mve')
        mve_vecs = embedding['Y'].T
        pca_vecs = embedding['YPCA'].T


        return spairs, squeeze(clusters['inds']), pca_vecs, mve_vecs


def get_consensus(ofs = 0,     
                  mweight = .5, 
                  refseq_method = 'root',
                  sp_method = 'sample',
                  aff_type = 'pairs',
                  reset = False,
                  do_plot = False,
		  run_id = 'CONS_TEST'):

    #FETCH FAMILY DATA
    ali, tree, infos, rfid = get_fam(ofs)

    
    #SPLIT THE FAMILY INTO SUBTREES
    minterms = 10
    csize = 20                        #appx clade size for alifold
    tsize = len(tree.get_terminals()) #Tree size
    tweight=tree.total_branch_length()
    cweight=tweight*csize/tsize

    cqueue = [tree.root]
    clades = []
    #SOMETHING IS BAD IN CQUEUE
    while cqueue:
        clade = cqueue.pop()
        if clade.total_branch_length() <= cweight\
                or len(clade.get_terminals()) <= minterms:
            clades.append(clade)
        else: 
            [cqueue.append(c) for c in clade.clades]            
    ali_ids = [a.name for a in ali]

    for i, n in enumerate(tree.get_terminals()):
        term_id = re.compile('_([^\.]*\.\d*)').search(n.name).group(1) 
        this_seq = ali[ali_ids.index(term_id)]
        n.m = {'seq':this_seq,
               'probs':[1 for j in range(len(this_seq))]}


    big_refnode, big_refseq = \
        subtree_refseq(tree, method = refseq_method)
    ungapped_ref = Bio.SeqRecord.SeqRecord( Bio.Seq.Seq(''.join([let for let in str(big_refseq.seq) if let.upper() in 'AUCGT']),
                                                        big_refseq.seq.alphabet)
                                            ,rfid + '_ref', name = big_refnode.name)

    title = 'RID_{3}_{0}_subopts_for_{1}_emb_{2}'.format(sp_method,rfid,aff_type,run_id)
    spairs, inds, pca_vecs, mve_vecs = mem.getOrSet(setAffinities,**mem.rc({},reset = reset,
									   tree = tree,
									   aff_type = aff_type,
									   sp_method = sp_method,
									   refseq = ungapped_ref,
									   on_fail = 'compute',
									   hardcopy = True,
									   register = title,
									   ))

    if do_plot : plot_clusters(inds,{'pca embedding':pca_vecs},title = title,plot3d = True)
    exemplar_structs = [spairs[e] for e in set(inds)]
    struct_profiles = infernal_profiles(ungapped_ref,exemplar_structs)



    all_vecs = {'all_time':[ [ [] for i in range(len(struct_profiles))] 
			     for j in range(len(clades)) ],
		'all_mut':[ [ [] for i in range(len(struct_profiles))] 
			     for j in range(len(clades)) ],
		'fiftyfifty':[ [ [] for i in range(len(struct_profiles))] 
			     for j in range(len(clades)) ]}
    aamuts, aatimes = [], []
    for idx_clade, c in enumerate(clades):
        if len(c.get_terminals()) < 3:
		print 'SKIPPPING CUZ SUBTREE TOO SMALL'
		continue
	c_ids = [ n.m['seq'].name for n in c.get_terminals() ]
	if len(nonzero(greater([len(list(g)) for k, g in it.groupby(sorted(c_ids))],1))[0])>0:
		print 'SKIPPING CUZ THERE ARE TWO COPIES OF SOME FUCKING SEQUENCE IN TREE'
		continue
           
        all_muts = []
	all_times = []  
	print
	print 'Clade: {0}'.format(idx_clade)
        for idx_struct, struct_info in enumerate( zip( struct_profiles, exemplar_structs)):
          struct_profile, ex_struct = struct_info
	  ngaps = 0
	  sys.stdout.write('')
	  sys.stdout.flush()


          #ALTHOUGH I AM CONSTRUCTING NEW ALIGNMENTS,
          #I AM NOT RECONSTRUCTING THE TREE AFRESH.
          #THIS IS A MISTAKE?
          
          #OLD ALIGNMENTS
          calis = ba.MultipleSeqAlignment(\
              [ali[ali_ids.index(c_id)] for c_id in c_ids ])
          #NEW ALIGNMENTS AND REF STRUCTURE
          c_new_ali , stk, struct = infernal_alignment(calis, struct_profile, rfid)
          #REF STRUCTURE PAIRS
          pairs = stk_pairs(struct)
	  if len(pairs) != len(ex_struct):
		  raise Exception()
           
          cterms = c.get_terminals()
          for i2, ct in enumerate(cterms):
              lilid =  'N{0}'.format(i2)
              ct.name = lilid
              ct.m['str_seq'] = c_new_ali[i2]
              ct.m['str_seq'].id = lilid
          
          #Create a tree from the current clade and run aligner.
	  #DOES THIS ALL WORK FINE WITH 'STR_SEQ'
          tr = phy.BaseTree.Tree(c)

          paml_run_id = 'ali_anc_c{0:04}_s{0:03}'.format(idx_clade,idx_struct)
          rstfile= paml.run_paml(tr, c_new_ali, run_id = paml_run_id)
          #Get ML ancestor sequences.
          anc_tree = paml.rst_parser(rstfile) 

          #Label extent and internal nodes with sequences.
          for term in anc_tree.get_terminals():
              #Terminals have old (rfam) alis and new (infernal) alis
              term.m = filter( lambda x: x.name == term.name, cterms)[0].m
          for node in anc_tree.get_nonterminals():
              #Internals only have new alis. m['seq'] = m['str_seq']
              node.m['str_seq'] = node.m['seq']
              node.m['str_seq'].seq = node.m['str_seq'].seq.replace('T', 'U')
          subtree = anc_tree
              
 
          #Evaluate all of the structs on the first pass
          #to have access to mean frequencies of different
          #mutational types in the final score computation
          
          refnode, refseq = subtree_refseq(subtree, method = refseq_method)
          muts, times = subtree_count_struct(subtree, pairs)
          all_muts.append(muts)
          all_times.append(times)
        
        #Compute a muliplier for comp/wobble muts from their relative
        #frequencies over all suboptimal structures.
        nwob, ncomp, nucom, nreco, nbbad  =[  sum([sum(m[k]) for m in all_muts]) 
                                              for k in ['wob','comp','ucom','reco','bbad']]
	nrm = ncomp if ncomp != 0 else 1
        comp_bonus = min([10.,max([2.,nwob/nrm])])
        
        for i, struct in enumerate(exemplar_structs):          
          #Compute the vector leaving mweight at its default value
		
          vecs = [subtree_rate_struct_V0(struct,ungapped_ref,
                                         all_muts[i], all_times[i],
                                         comp_bonus = comp_bonus,
                                         mweight =mw) for mw in [0.,1.,.5]]
          for k,v in zip(['all_time','all_mut','fiftyfifty'],vecs): 
		  all_vecs[k][idx_clade][i] = v
	  

	aamuts.append(all_muts)
	aatimes.append(all_times)
    outputs = {
	    'all_vecs':all_vecs,
	    'all_muts':aamuts,
	    'all_times':aatimes,
	    'exemplar_structs':exemplar_structs,
	    'reference_seq':ungapped_ref,
	    'thermo_pairs':spairs,
	    'thermo_ex_inds':inds,
	    'thermo_embedding':pca_vecs,
	    'title':title,
	    'thermo_aff_type':aff_type,
	    'tree':tree,
	    'run_id':run_id
	    }
	 
    import pickle
    pickle.dump(outputs, open(cfg.dataPath('cs874/runs/{0}.pickle'.format(run_id)),'w'))
    return(outputs)

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


def subtree_rate_struct_V0( pairs, refseq,
                           muts, times,
                           mweight = .5,
                           comp_bonus = None):
    '''
Given the kinds of mutation that occur in a tree of sequences 
with respect to an RNA structure, compute a vector of length
l representing the likelihood that evolution has acted to
preserve base pairing according to struct.

Combines mutational and time-pairing data according to kwargs.

inputs:
  subtree: (object) biopython tree with nodes having attribute 'm',
                    a dictionary with key 'seq'.
  pairs:   []       pairs specified by struct in seq

  muts:    {}       mut types and counts for p in pairs.
  times:   {}       branch lengths paired, unpaired for p in pairs.

keywords:
  mweight: (float)  how much weight to give the mutation metric.

                    tweight = 1 - mweight

  comp_bonus:(float) how much more weight to give to double (comp)
                    mutants, compared to single (wobble) mutants.

                    comp_bonus = None => use (n_wob)/(n_comp) 
                                         clipped to [2.0, 10.0]
'''

    assert mweight <= 1 and mweight >= 0
    tvec = zeros(len(refseq))
    mvec = zeros(len(refseq))

    #In V0 we don't really penalize mispairs substantially. 
    #They only appear in the normalization terms by way of 
    #summing over all muts.values() in the denom of mval 
    #and computing tval as paired / total.

    #RATE PAIRED POSITIONS
    for i, p in enumerate(pairs):
        nrm = (sum([mtype[i]
                            for mtype in muts.values()])\
                           + muts['comp'][i] * (comp_bonus -1))
        nrm += equal(nrm,0)
        mval = float((muts['comp'][i] * comp_bonus \
                    + muts['wob'][i] \
                    + muts['reco'][i])) \
                    / nrm
	
	nrm =  times['total'][i] if times['total'][i] != 0 else 1.
        tval = (times['paired'][i] / nrm) 
        for elt in p:
            mvec[elt] = mval
            tvec[elt] = tval
    
    return mvec*mweight + tvec*(1 - mweight)


def subtree_count_struct(subtree, pairs):
   '''
Given a tree of sequences and a fixed RNA structure specified
as a list of paired bases, compute:
  
  1) muts:  the frequencies of various mutation types at
            for each pair at each tree edge.
  2) times: the total amount of time spent paired/unpaired
            for each pair over all tree edges.


inputs:
  subtree: (object) biopython tree with nodes having attribute 'm',
                    a dictionary with key 'seq'.
  pairs:   []       pairs specified by struct in seq

keywords:
  [NONE]

outputs:
  muts:    Mutation counts for each pair.  
  times:   Times paired/unpaired for each pair.

'''
    
   #CREATE A PARENTS ATTRIBUTE TO USE FOR 
   #SCANNING PAIRS.
   for n in subtree.get_nonterminals():
       if n == subtree.root: n.parent = None
       else: n.parent = ([subtree.root] + subtree.root.get_path(n))[-2] 
   for n in subtree.get_terminals():
       pass
       n.parent =([subtree.root] +  subtree.root.get_path(n))[-2]      
   
   #DATA MATRICES COUNTING SUB TYPES
   m_wob, m_ucom, m_comp, m_reco, m_bbad = \
       [zeros(len(pairs)) for i in range(5)] 
   #DATA MATRICES COMPUTING PAIRING OVER TIME
   t_up, t_p, t_tots = \
       [zeros(len(pairs)) for  i in range(3)]
   muts = dict(wob = m_wob, ucom = m_ucom, comp = m_comp, 
               reco = m_reco, bbad = m_bbad)
   times= dict(unpaired = t_up, paired = t_p,
               total = t_tots)

   gap_count = 0
   for idx,p in enumerate(pairs):
       #DEAL WITH GAPS (SKIP)
       if '-' in list(it.chain(*[(t.m['str_seq'].seq[p[0]], \
                                      t.m['str_seq'].seq[p[1]])
                                 for t in subtree.get_terminals()])):
	   gap_count += 1

           continue

       irresolvables = []
       curterms= set(subtree.get_terminals())
       curints = set(subtree.get_nonterminals())
       loops = 0
       while 1:
           loops += 1
           ct_list = sorted(list(curterms), 
                            key = lambda x: x.parent)
           preterms =[(k,list(g))
                      for k,g in it.groupby(ct_list,lambda x: x.parent)]
           preterms =filter(lambda x: len(x[1]) == 2,preterms)
           if not preterms:
               break

                     
           for pp, terms in preterms:
               tprobs = [(t.m['probs'][p[0]], t.m['probs'][p[1]]) 
                          for t in terms]
               tbases = [(t.m['str_seq'][p[0]], t.m['str_seq'][p[1]])
                          for t in terms]

               pprobs = pp.m['probs'][p[0]],pp.m['probs'][p[1]]
               pbases= pp.m['str_seq'][p[0]],pp.m['str_seq'][p[1]]
               
               def do_contract(curterms, curints,terms,pp):
                   curterms.difference_update(terms)
                   curints.difference_update([pp])
                   curterms.update([pp])
                   return array([t.branch_length for t in terms])
               def do_giveup(curterms,curints, termps, pp):
                   curterms.difference_update(terms)
                   curints.difference_update([pp])
                   
               #ALL ELTS PRESERVED
               if tbases[0] == tbases[1] == pbases:
                   lc = sum(do_contract(curterms,curints,terms,pp))
                   t_tots[idx] += lc
                   if ''.join(tbases[0]) in globals()['bps']:
                       t_p[idx] += lc
                   else:
                       t_up[idx] += lc
                       
               #ELEMENTS CHANGED, GOOD MARGINAL RECO
               elif sum([1.-prob for prob in pprobs]) < .30:
                   
                   #NEITHER CHILD MATCHES THE PARENT
                   if not sum([tb == pbases for tb in tbases]):
                       #PARENT HAS ONE BASE FROM EACH CHILD
                       if pbases in [ (tbases[0][0],tbases[1][1]),
                                      (tbases[1][0],tbases[0][1])]:
                           m_wob[idx] += 2
                           s=do_contract(curterms,curints,terms,pp)
                           for i, tb in enumerate( tbases ):
                              if tb in globals()['bps']:
                                  t_p[idx] += s[i]
                              else:
                                  t_up[idx]+= s[i]
                           t_tots[idx] += sum(s)

                       #PARENT HAS <1 BASE PER CHILD
                       else:
                           #THIS SHOULD ACTUALLY TAKE OVER....
                           #WHEN I HAVE THE TIME, ELIMINATE THE 
                           #REST OF THIS COMPUTATION AND USE THIS
                           pseq = ''.join(pbases)
                           cseqs = [''.join(tb) for tb in tbases]
                           s = do_contract(curterms, curints,terms, pp)
                           for j, cseq in enumerate(cseqs):
                               ppar = pseq in globals()['bps']
                               cpar = cseq in globals()['bps']
                               diffs= [2-len(set(z)) for z in zip(pseq,cseq)]
                               if ppar:
                                   if cpar:
                                       if sum(diffs)==0:
                                           pass #no mutation
                                       elif sum(diffs)==1:
                                           m_wob[idx] += 1
                                       else:
                                           m_comp[idx] +=1
                                       t_p[idx] += s[j];
                                   else:
                                       if sum(diffs)==0:
                                           pass #never reached
                                       elif sum(diffs)==1:
                                           m_ucom[idx] +=1
                                       else:
                                           m_ucom[idx] +=1 #SHOULD BE 2?
                                       t_up[idx] += s[j]
                               else:
                                   if cpar:
                                       if sum(diff)==0:
                                           pass #never reached.
                                       elif sum(diffs) ==1:
                                           m_reco[idx] +=1
                                       else:
                                           m_reco[idx] +=1 #SHOULD BE 2?
                                       t_p[idx] += s[j]
                                   else:
                                       if sum(diffs) == 0:
                                           pass
                                       elif sum(diffs) == 1:
                                           m_bbad[idx] += 1
                                       else:
                                           m_bbad[idx] += 1 #SHOULD BE 2?

                   
                   #PARENT MATCHES ONE CHILD EXACTLY
                   else:
                     matched = tbases.index(pbases)
                     um =[''.join(x) for x in [pbases, tbases[1-matched]]]
                     
                     #COMPENSATED MUTATION TRANSPOSING THE PARENTAL PAIR 
                     if \
                             'GC' in um and 'CG' in um \
                             or 'UA' in um and 'AU' in um \
                             or 'UG' in um and 'GU' in um:
                         lc = sum(do_contract(curterms,curints,terms,pp))
                         t_p += lc
                         t_tots[idx] += lc
                         m_comp[idx] += 1
                         
                     #MUTATION CHANGING AA IDENTITY OF PAIRED PARENT
                     elif um[0] in globals()['bps']:
                         
                         #PAIRING RETAINED
                         if um[1] in globals()['bps']:
                             #WOBBLE SILENCED (1 CHANGE)
                             if um[0][0]==um[1][0] or um[0][1] == um[1][1]:
                                 m_wob[idx] += 1
                             #COMP SILENCED   (2 CHANGES)
                             else:
                                 m_comp[idx] += 1
                             lc = sum(do_contract(curterms,curints,terms,pp))
                             t_tots[idx] += lc
                             t_p[idx] += lc

                         #PAIRING DESTROYED
                         else:
                             lens =do_contract(curterms,curints,terms,pp)
                             t_up[idx] += lens[1-matched]
                             t_p[idx] += lens[matched]
                             t_tots[idx] += sum(lens)
                             m_ucom[idx] += 1

                     #PARENT DOES NOT PAIR
                     else:
                         if um[1] in globals()['bps']:
                             #A RECOMPENSATING MUTANT
                             lens = do_contract(curterms,curints,terms,pp)
                             t_up[idx] += lens[matched]
                             t_p[idx] += lens[1-matched]
                             t_tots[idx] += sum(lens)
                             m_reco[idx] +=1
                         else:
                             #NO PAIRING ANYWHERE
                             lens = do_contract(curterms,curints,terms,pp)
                             t_up[idx] += sum(lens)
                             t_tots[idx] += sum(lens)
                             m_bbad[idx] +=1
               else:   
                   irresolvables.append(terms)
                   do_giveup(curterms,curints,terms,pp)
       
       #CLOSES SEARCH LOOP FOR P
   if len(muts.values()[0]) != len(pairs) :
	   raise Exception()

   restart_line()
   sys.stdout.write('gap count: {0:03}'.format(gap_count))
   sys.stdout.flush()
   #CLOSES SEARCH LOOP FOR ALL PAIRS                   
   return muts, times


def show_output(outputs, 
		show = 'conservation',
		save = True):
	mvecs = outputs['all_vecs']['all_time']
	tvecs = outputs['all_vecs']['all_mut']
	fvecs = outputs['all_vecs']['fiftyfifty']

	run_id = outputs['run_id']
	structs = outputs['exemplar_structs']
	ref = outputs['reference_seq']
	
	thermo_pairs = outputs['thermo_pairs']
	thermo_inds  = outputs['thermo_ex_inds']

	run_title = outputs['title']

	fig = plt.gcf()
	try: fig.clear()
	except Exception, e: print 'wonky 3d bug'

	exemplar_inds = sorted(list(set(thermo_inds)))
	struct_colors = dict([(exemplar_inds[i], col) 
			       for i, col in enumerate(mycolors.getct(len(exemplar_inds)))]
                              )


	if show == 'embeddings':

	   
	   exemplars = list(set(thermo_inds))
	   pair_embedding =  compute_embedding(thermo_pairs,
	   			       aff_type = 'pairs',
	   			       do_mve = False,
					       ss_multiplier = None)
	   
	   shape_embedding = compute_embedding(thermo_pairs,
	   			       aff_type = 'easy',
	   			       do_mve = False,
	   			       ss_multiplier = None)
	   show_3d = True
	   #shape_embedding[0] is pca
	   plot_clusters( thermo_inds, {'shape':shape_embedding[0],
	   			     'pairs':pair_embedding[0]}, 
			  plot3d = show_3d,
			  title = 'projection ({0}) '.format(run_id),
			  save = save,
			  colors = struct_colors)

	elif show == 'conservation':
		ax0 = fig.add_subplot('311')
		lstructs =  [project_lstruct(p, len(ref)) for p in structs]
		seismic.seismic([ abs(l) for l in lstructs] , 
				colors = struct_colors.values(),
				ax = ax0)

		myplots.maketitle(ax0, 'Predicted conservation patterns for {0}'.format(run_title))

		shapes = array([shape(m) for m in mvecs])
		igood = nonzero(greater(shapes[:,1],0))[0]
		clade_colors = mycolors.getct(len(igood))
		mvg, tvg, fvg = [ [vecs[i] for i in igood] for vecs in [mvecs,tvecs,fvecs]]
		cons_types = array([ mvg, tvg, tvg])
		
		for c in cons_types:
			nrm = sum(c.flatten())
			if nrm == 0: nrm = 1
			c /= sum(c.flatten())
		
		mtype_sums = np.sum(np.sum(cons_types,3),0)	
		stype_sums = np.sum(np.sum(cons_types,3),0).T


		ax1 = fig.add_subplot('312')		
		seismic.seismic(stype_sums , 
				colors = struct_colors.values(),
				ax = ax1)

		myplots.maketitle(ax1,'Observed conservation (struct v. clade) patterns for {0}'\
					  .format(run_title),
				  )

		
		ax2 = fig.add_subplot('313')
		
		seismic.seismic(mtype_sums , 
				ax = ax2, colors = clade_colors, stacked = True,
				label_y = False)

		#myplots.maketitle(ax2, 'Observed conservation (clade v. struct) patterns for {0}'\
		#			   .format(run_title)
		#		   )
		ax2.annotate('Observed conservation (clade v. struct) patterns for {0}'\
				     .format(run_title),
			     [.5,0],xycoords = 'axes fraction', ha = 'center', va = 'top',
			     size = 'x-large')

		if save: fig.savefig(cfg.dataPath('cs874/figs/cons_profiles/{0}.ps'.format(run_title)))
	       				 	



	else: raise Exception('show type not implemented: {0}'.format(show))
	
def project_lstruct(pstruct, l):
	out = zeros(l)
	for p in pstruct:
		out[p[0]] = -1.
		out[p[1]] = 1.
	return out
	
def showseq(tree, elts):
    phy.draw_graphviz(tree, lambda x:''.join([ x.m['seq'].seq[e] for e in elts]))

bps = ['GC','AU','GU','CG','UA','UG']


import compbio.utils.bsub as bsub
import compbio.utils.bsub_utils as bsu

def runmany(run_id):
	print 'TESTING WITH A LIMITED RANGE OF FAMILIES'
	inp_dicts = [dict([('ofs',r)]) for r in range(0,1493)][0:10]
	eyeball = bsub.eyeball(run_id, 
			       os.path.abspath(inspect.stack()[0][1]),
			       inp_dicts,
			       func = 'run',
			       name = 'ra2_runs_',
			       mem = 3)
	raise Exception()
	eyeball.launch()
	return 'sxs'
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
    
