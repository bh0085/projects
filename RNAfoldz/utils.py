import subprocess as spc, sys, re, os,inspect, itertools as it
import matplotlib.mlab as mlab
import Bio.AlignIO as aio
import Bio
import compbio.utils.colors as mycolors
import compbio.utils.seismic as seismic
import compbio.utils.plots as myplots
import compbio.config as cfg


import rfam

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
        
        print lines
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
    return project_lstruct(pairs, l)

def struct_project_l3(pairs,l):
    return stk_mat(pairs_stk(pairs,l))

def cluster_2(structs, polys, seq):
    mats =array([struct_project_l3(p, len(seq)) for p in structs])
    vecs =array([struct_project_l(p, len(seq)) for p in structs])

    polys = array(polys)
    sortorder = argsort([sum(mats[0] * mats[i]) for i in range(len(polys))])[::-1]

    affinities, ss = affinity_matrix(structs, aff_type = 'easy')
    #aff_shape, ss_shape = affinity_matrix(structs, aff_type = 'easy', ss_multiplier = .5)
    
    pca_vecs = mlab.PCA(affinities).project(affinities)  
    #pca_vecs_shape = mlab.PCA(aff_shape).project(aff_shape)  

    sortorder2 = argsort(pca_vecs[:,0])
    sortorder = sortorder2

    import mlpy

    HC = mlpy.HCluster(method='euclidean', link='complete')

    #cvecs = vecs
    cvecs = pca_vecs[:,0:5]
    
    #plt.gcf().clear()
    #plt.plot(cvecs[sortorder])
    #return
    
    clusts = HC.compute(cvecs )#/ sum(cvecs **2, 1)[:,newaxis])
    cut = HC.cut(HC.heights[::-1][3])

    ct = mycolors.getct(len(mats))
    rplots.grid_rnas(polys[sortorder], 
                     colors = [ct[i] for i in cut[sortorder]],
                     size = (5,5))

    return HC

    
def suboptimals(sequence, sp_method = 'enumerate', n = 400):
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


def family_suboptimals(rfid, plots = True):
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

    spairs = suboptimals(ungapped_ref, sp_method = 'sample', n = 400)
    if plots:
        verts = [rna_draw(ungapped_ref.seq, 
                              pairs_stk(sp,len(ungapped_ref)),
                              'name' )
                 for sp in spairs]
    else:
        verts = None
    return spairs, verts, ungapped_ref


def affinity_matrix(pairs, aff_type = 'pairs', ss_multiplier = None):

	#Artificial sequence length from the maximum paired element
	seq_len = np.max(list(it.chain(*pairs))) + 1

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
	return affinities, ss
		       
