'''
Given an rfam family, its alignment and the tree given, compute ancestral nodes
and mutations along the edges.

'''
import rfam
import utils as rutils
import infernal
import Bio.AlignIO as aio
import Bio.Align as ba
import Bio.Phylo.NewickIO as nio
import Bio.Phylo as phy
import Bio
import subprocess as spc, sys, re, os,inspect,pickle
import pickle
import compbio.utils.bsub as bsub
import compbio.utils.bsub_utils as bsu
import compbio.utils.memo as mem
import compbio.utils.bs_macros as bsm
import compbio.config as cfg
import matplotlib.mlab as mlab

if __name__ != '__main__' and True:
        import plots as rplots

from numpy import *
import numpy as np
import itertools as it
from compbio.projects.seqtree import muscle, phyml, paml
import sys

bps = ['GC','AU','GU','CG','UA','UG']

 
   
def compute_clusters(affinities, ss):
        cluster_dict = dict(similarities=affinities,
                           self_similarity = ss)
        clusters =  bsm.runmat('ap_frompy', cluster_dict, 'subopt_clusters')
        return squeeze(clusters['inds'])

def compute_embedding(affinities,
		      aff_type = 'pairs',
		      do_mve = False,
		      ss_multiplier = None):


	mve_dict = dict(k = 4,
			similarities = affinities,
                        do_mve = False if aff_type == 'easy' else False)
        embedding=  bsm.runmat('mve_frompy', mve_dict, 'subopt_mve')
        mve_vecs = embedding['Y'].T
        pca_vecs = embedding['YPCA'].T

	return pca_vecs, mve_vecs


def split_tree(tree):
   '''
Split a bipython tree in to similarly sized subtrees.
'''
   minterms = 10
   csize = 15                        #appx clade size for alifold
   tsize = len(tree.get_terminals()) #Tree size
   tweight=tree.total_branch_length()
   cweight=tweight*csize/tsize
   cqueue = [tree.root]
   clades = []
   while cqueue:
        clade = cqueue.pop()
        if clade.total_branch_length() <= cweight\
                or len(clade.get_terminals()) <= minterms:
            clades.append(clade)
        else: 
            [cqueue.append(c) for c in clade.clades]            
   return clades

       
def family_exemplar_structs(rfid,
                            refseq_method = None,
                            sp_method = None,
                            aff_type = None):

    suboptimals = rutils.family_suboptimals(rfid)
    c2 = rutils.cluster_2(spairs, ungapped_ref)

    arr = rutils.rna_draw(ungapped_ref.seq, 
                          rutils.pairs_stk(sp,len(ungapped_ref)),
                          'name' )

    raise Exception()
    affinities, ss = rutils.affinity_matrix(spairs, aff_type = aff_type)
    aff_shape, ss_shape = rutils.affinity_matrix(spairs, aff_type = 'easy', ss_multiplier = .5)
    

    pca_vecs = mlab.PCA(affinities).project(affinities)  
    pca_vecs_shape = mlab.PCA(aff_shape).project(aff_shape)  
    inds = compute_clusters(aff_shape, ss_shape)
    exemplars = list(set(inds))
    
    import compbio.utils.colors as mycolors
    ct = mycolors.getct(len(exemplars))
    import matplotlib.pyplot as plt
    f = plt.gcf()
    plt.clf()
    
    for idx0, embeddings in enumerate([pca_vecs, pca_vecs_shape]):
            ax = f.add_subplot('21{0}'.format(idx0 +1))

            lims =[ [min(embeddings[:,0]),max(embeddings[:,0])],
                         [min(embeddings[:,1]),max(embeddings[:,1])] ]
            lims += [-.5,.5] *squeeze(diff(lims,1))[:,newaxis]
            

            ax.set_xlim(lims[0])
            ax.set_ylim(lims[1])
    
            print sum(embeddings)
            for idx, embedding in enumerate(embeddings):
              if mod(idx,1) != 0: continue
              sp = spairs[idx]
              arr = rutils.rna_draw(ungapped_ref.seq, 
                              rutils.pairs_stk(sp,len(ungapped_ref)),
                              'name' )
              struct_emb = arr + embedding[0:2]
              #plt.plot(*struct_emb.T)
              
              pkw = {'color':ct[exemplars.index(inds[idx])],
                     'lw':8 if idx in inds else 1,
                     'alpha': 1 if idx in inds else .2}
              
              lc = rplots.show_rna(embedding, arr, pkw = pkw)
    #exemplar_structs = [spairs[e] for e in set(inds)]  
    raise Exception()

    return pca_vecs, exemplar_structs

def get_consensus(rfid = 'RF00', mweight = .5, 
                  refseq_method = 'root', sp_method = 'sample',
                  aff_type = 'pairs',  reset = True,
                  do_plot = False,  run_id = 'CONS_TEST'):

    ali, tree, infos = rfam.get_fam(rfid)
    ali_ids = [a.name for a in ali]

    for i, n in enumerate(tree.get_terminals()):
        term_id = re.compile('_([^_]*)_').search(n.name).group(1) 
        this_seq = ali[ali_ids.index(term_id)]
        n.m = {'seq':this_seq,
               'probs':[1 for j in range(len(this_seq))]}

    #if do_plot : rplots.plot_clusters(inds,{'pca embedding':pca_vecs},title = title,plot3d = True)
    

    big_refnode, big_refseq = \
        subtree_refseq(tree, method = refseq_method)
    ungapped_ref = rutils.ungapped_seq(big_refseq, rfid)
    #pca_vecs,exemplar_structs =
    return family_exemplar_structs(rfid,
                                   sp_method = sp_method,
                                   refseq_method = refseq_method,
                                   aff_type = aff_type,
                                   )
    struct_profiles = infernal.profiles(ungapped_ref,exemplar_structs, run_id)

    clades = split_tree(tree)
    all_vecs = {'all_time':[ [ [] for i in range(len(struct_profiles))] 
			     for j in range(len(clades)) ],
		'all_mut':[ [ [] for i in range(len(struct_profiles))] 
			     for j in range(len(clades)) ],
		'fiftyfifty':[ [ [] for i in range(len(struct_profiles))] 
			     for j in range(len(clades)) ]}

    aamuts, aatimes, aairr, aagaps = [], [], [], []
    for idx_clade, c in enumerate(clades):
        if len(c.get_terminals()) < 3:
		print 'SKIPPPING CUZ SUBTREE TOO SMALL'
		continue
	c_ids = [ n.m['seq'].name for n in c.get_terminals() ]
	if len(nonzero(greater([len(list(g)) for k, g in it.groupby(sorted(c_ids))],1))[0])>0:
		print 'SKIPPING CUZ THERE ARE TWO COPIES OF SOME FUCKING SEQUENCE IN TREE'
		continue          
        all_muts, all_times , all_gaps, all_irr = [], [], [], []
	print
	print 'Clade: {0}'.format(idx_clade)
        for idx_struct, struct_info in enumerate( zip( struct_profiles, exemplar_structs)):
          struct_profile, ex_struct = struct_info
	  ngaps = 0

          #OLD ALIGNMENTS
          calis = ba.MultipleSeqAlignment(\
              [n.m['seq'] for n in c.get_terminals() ])
          #NEW ALIGNMENTS AND REF STRUCTURE
          c_new_ali , stk, struct = infernal.alignment(calis, struct_profile, rfid)
          #REF STRUCTURE PAIRS
          pairs = rutils.stk_pairs(struct)
	  if len(pairs) != len(ex_struct):
		  raise Exception()
           
          cterms = c.get_terminals()
          for i2, ct in enumerate(cterms):
              lilid =  'N{0}'.format(i2)
              ct.name = lilid
              ct.m['str_seq'] = c_new_ali[i2]
              ct.m['str_seq'].id = lilid
	      ct.m['probs'] = ones(len(c_new_ali[i2]))
          
          #BUILD A TREE
          tr = phy.BaseTree.Tree(c)

          #RUN PAML
          paml_run_id = 'ali_anc_c{0:04}_s{0:03}'.format(idx_clade,idx_struct)
          rstfile= paml.run_paml(tr, c_new_ali, run_id = paml_run_id)
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
          muts, times, gaps, irresolvables = subtree_count_struct(subtree, pairs)
          all_muts.append(muts)
          all_times.append(times)
	  all_gaps.append(gaps)
	  all_irr.append(irresolvables)
        
	compute_signatures(all_vecs,idx_clade,
			   all_muts,all_times,
			   exemplar_structs,ungapped_ref )
				      
	aamuts.append(all_muts)
	aatimes.append(all_times)
	aairr.append(all_irr)
	aagaps.append(all_gaps)
    outputs = {
	    'all_vecs':all_vecs,
	    'all_muts':aamuts,
	    'all_times':aatimes,
	    'exemplar_structs':exemplar_structs,
	    'reference_seq':ungapped_ref,
	    'thermo_ex_inds':inds,
	    'thermo_embedding':pca_vecs,
	    'title':title,
	    'thermo_aff_type':aff_type,
	    'tree':tree,
	    'run_id':run_id
	    }
	 
    pickle.dump(outputs, open(cfg.dataPath('cs874/runs/{0}.pickle'.format(run_id)),'w'))
    return(outputs)
def compute_signatures(all_vecs,idx_clade,all_muts, all_times, structs, reference):
        #Compute a muliplier for comp/wobble muts from their relative
        #frequencies over all suboptimal structures.
        nwob, ncomp, nucom, nreco, nbbad  =[  sum([sum(m[k]) for m in all_muts]) 
                                              for k in ['wob','comp','ucom','reco','bbad']]
	nrm = ncomp if ncomp != 0 else 1
        comp_bonus = min([10.,max([2.,nwob/nrm])])
        
        for i, struct in enumerate(structs):          
          #Compute the vector leaving mweight at its default value
		
          vecs = [subtree_rate_struct_V0(struct,reference,
                                         all_muts[i], all_times[i],
                                         comp_bonus = comp_bonus,
                                         mweight =mw) for mw in [0.,1.,.5]]
          for k,v in zip(['all_time','all_mut','fiftyfifty'],vecs): 
		  all_vecs[k][idx_clade][i] = v




def subtree_rate_struct_V0( pairs, refseq, muts, times,
                           mweight = .5, comp_bonus = None):
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
   
   gaps = zeros(len(pairs))
   irresolvables = [[] for i in range(len(pairs))]
   muts = dict(wob = m_wob, ucom = m_ucom, comp = m_comp, 
               reco = m_reco, bbad = m_bbad)
   times= dict(unpaired = t_up, paired = t_p,
               total = t_tots, total_incl_unresolved = subtree.total_branch_length() -\
		       subtree.root.branch_length)

   gap_count = 0
   for idx,p in enumerate(pairs):
       #DEAL WITH GAPS (SKIP)
       if '-' in list(it.chain(*[(t.m['str_seq'].seq[p[0]], \
                                      t.m['str_seq'].seq[p[1]])
                                 for t in subtree.get_terminals()])):
	   gap_count += 1
	   gaps[idx] = 1
           continue
       
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
                   irresolvables[i].append(terms)
                   do_giveup(curterms,curints,terms,pp)   
   return muts, times, gaps, irresolvables


def get_remote_runs(run_range):
	outdir = cfg.dataPath('batch/outputs')
	files = [os.path.join(outdir, f) for f in os.listdir(outdir) if 'ra2_' in f]
	outs = [ pickle.load(open(f)) for f in files[0:20] ]
	return outs
def draw_remote_runs(show = 'conservation'):
	outdir = cfg.dataPath('batch/outputs')
	files = [os.path.join(outdir, f) for f in os.listdir(outdir) if 'ra2_' in f][1:]
	for idx, f in enumerate(files):
		print '{0} of {1} files'.format(idx,len(files))
		print f[-100::]
		fopen = open(f)
		out = pickle.load(fopen) 
		if transform:
			'''Fix stuff'''
			out_t = out
		else: out_t = out
		rplots.show_output(out_t)
		fopen.close()
	return outs	
	
