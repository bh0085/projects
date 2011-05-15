import compbio.utils.bsub_utils as butils
import compbio.utils.colors as mycolors
import matplotlib.pyplot as plt
import compbio.utils.plots as myplots

import utils, tree_conservation
import compbio.utils.memo as mem
import rfam, infernal, phase
from numpy import *
import numpy as np
import matplotlib.mlab as mlab

import compbio.projects.seqtree.phyml as phyml
import compbio.projects.seqtree.paml as paml
import compbio.projects.seqtree.muscle as muscle
import compbio.utils.seismic as seismic

import Bio.Align as align
import Bio.Phylo as phylo
from scipy.stats import linregress
import compbio.config as cfg
import hcluster 
import itertools as it

draw_all_easy= False
draw_all_hard =False
easy_way = True
if easy_way:
    clade_tree_method = 'bionj'
    clade_alignment_method = 'cm'
    clade_ancestor_method = 'independent'
else:
    raise Exception('Sorry but the hard way is not done yet')
def seq_dists(ali,run_id, tree = True):
    import Levenshtein
    n = len(ali)
    dists = zeros((n,n))

    if tree:
        ali_named = align.MultipleSeqAlignment(ali)
        maps = {}
        for idx, a in enumerate(ali_named):
            a.id = 'S{0:05}'.format(idx) 
            maps[a.id] = idx
        tree = phyml.tree(ali_named, run_id = run_id, bionj = True)
        for n1 in tree.get_terminals():
            for n2 in tree.get_terminals():
                dists[maps[n1.name],maps[n2.name]] = \
                    tree.distance(n1,n2)
    else:
        for i in range(n):
            for j in range(i):
                dists[i,j] = Levenshtein.distance(str(ali[i].seq), str(ali[j].seq))
                dists[j,i] = dists[i,j]
    return dists


def setDistances(ali = None,run_id = None, tree = True, **kwargs):
    assert ali != None ; assert run_id != None
    dists = seq_dists(ali, run_id, tree = True)
    return dists

def maxclust_dists(dists, k, method = 'complete'):
    d2 = hcluster.squareform(dists)
    Z = hcluster.linkage(d2, method = method)
    fcl = hcluster.fcluster(Z, t = k, criterion = 'maxclust')
    return fcl
                   

def get_seq_groups(rfid = 'RF00167', reset = True, tree = True,
        draw_distances = draw_all_easy,
        draw_clusters = draw_all_easy,
        draw_single_cluster = draw_all_hard):
    '''
Run the tree computation for each clsuter in the rfam family.
(Or just one)

1) Compute clusters using a distance measure derived either 
   phyml or a simple levenshtein dist.

   kwds:
     tree          [True]  Use a tree or just a levenshtein 
                           distance to get distances for
                           init clustering.

2) Choose a cluster of well related sequences and for this 
   this cluster, compute an alignment (For each structure 
   using phase or for sequences using MUSCLE)
  
   kwds:
     struct_align  [True]   Whether to compute structural 
                            alignments or use MUSCLE

'''
    rutils = utils

    ali, tree, infos = rfam.get_fam(rfid)
    n = len(ali)

    if draw_distances:
        dists_t = seq_dists(ali,rfid, tree = True)
        dists_l = seq_dists(ali,rfid, tree = False)
        dtf = dists_t.flatten()
        dlf = dists_l.flatten()
        lin = linregress(dtf, dlf)
        rsquared = lin[2]**2

        f = myplots.fignum(5, (8,8))
        ax = f.add_subplot(111)
        ax.annotate('Levenshtein distance vs. BioNJ branch lengths',
                    [0,1], xycoords = 'axes fraction', va = 'top',
                    xytext = [10,-10],textcoords = 'offset pixels')
        ax.annotate('R-Squared: {0}'.format(rsquared),
                    [1,0], xycoords = 'axes fraction', ha = 'right',
                    xytext = [-10,10],textcoords = 'offset pixels')
        ax.set_xlabel('BIONJ Tree ML Distance')
        ax.set_ylabel('Levenshtein Distance')

        ax.scatter(dtf, dlf, 100)
        
        datafile = cfg.dataPath('figs/gpm2/pt2_lev_tree_dists.svg')
        f.savefig(datafile)
        
    dists = mem.getOrSet(setDistances, ali = ali, tree = tree, run_id = rfid,
                         register = rfid, 
                         on_fail = 'compute',
                         reset = reset)
    
    clusters = maxclust_dists(dists, k = 5, method = 'complete')
    clusters -= 1

    if draw_clusters:

        ct = mycolors.getct(len(set(clusters)))
        colors = [ct[elt] for elt in clusters]
        pca_vecs = mlab.PCA(dists).project(dists) 
        
        f = myplots.fignum(5, (8,8))
        ax = f.add_subplot(111)
        ax.annotate('Rfam sequence clusters in first 2 PC of sequence space.',
                    [0,1], xycoords = 'axes fraction', va = 'top',
                    xytext = [10,-10],textcoords = 'offset pixels')
        ax.annotate('Number of Clusters: {0}'.format(len(ct)),
                    [1,0], xycoords = 'axes fraction', ha = 'right',
                    xytext = [-10,10],textcoords = 'offset pixels')
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')

        ax.scatter(pca_vecs[:,0],pca_vecs[:,1], 20, color = colors)
        
        datafile = cfg.dataPath('figs/gpm2/pt2_all_seqs_clustered.svg')
        f.savefig(datafile)        

    #now take the largest cluster and do the analysis.
    
    cgrps = dict([ (k, list(g)) 
              for k , g  in it.groupby(\
                sorted( list(enumerate(clusters)),key = lambda x: x[1]),
                key = lambda x: x[1])])
    cbig = argmax([len(x) for x in cgrps.values()])
    cluster_seqs = [ elt[0] for elt in cgrps.values()[cbig] ] 
    csize = len(cluster_seqs)
    seqs =[ali[c] for c in cluster_seqs]

    
    
    if draw_single_cluster:

        ct = mycolors.getct(2)
        pca_vecs = mlab.PCA(dists).project(dists) 
        colors =[ct[1] if elt in cluster_seqs else ct[0] for elt in range(len(pca_vecs))] 
        
        f = myplots.fignum(5, (8,8))
        ax = f.add_subplot(111)
        ax.annotate('Inter and intra cluster distances vs. PC0 component for chosen cluster.',
                    [0,1], xycoords = 'axes fraction', va = 'top',
                    xytext = [10,-10],textcoords = 'offset pixels')
        ax.annotate('Number of cluster sequences: {0}, Number of total sequences'.format(csize, n  - csize),
                    [1,0], xycoords = 'axes fraction', ha = 'right',
                    xytext = [-10,10],textcoords = 'offset pixels')
        ax.set_xlabel('PC 0')
        ax.set_ylabel('Distance')


        for s in cluster_seqs:
            ax.scatter(pca_vecs[:,0],dists[s,:] ,200 *exp(-(dists[s,:] / .5) **2),  color = colors, alpha = .2)
        
        datafile = cfg.dataPath('figs/gpm2/pt2_focused_cluster_dists.svg')
        f.savefig(datafile)        
        
    clusters_final  = [ [ elt[0] for elt in cgrps.values()[i] ] for i in range(len(cgrps.values()))]
    seqs_final = [ [ ali[idx] for idx in clust ] for clust in clusters_final]
    return seqs_final



    #Get alignments for each structure with muscle/structure
  
def tree_similarity(dist1, dist2, run_id,criterion = 'knn', k = 6):
    if criterion == 'knn':
        nq = len(dist1)
        nb1 = argsort(dist1, 1)[:,1:k+1]
        nb2 = argsort(dist2, 1)[:,1:k+1]
        all_nbs = [set(n1).union(set(n2)) for n1, n2 in zip(nb1, nb2)]
        nb_intersection = [set(n1).intersection(set(n2)) for n1, n2 in zip(nb1, nb2)]
        nb_dists = [ array([[dist1[i, n], dist2[i,n]]for n in nbs ]) for i,nbs in enumerate(all_nbs)]
        #take the first k distances.
        n_disagreements = [len(nbd) - k for nbd in nb_dists]
        nb_dists = array([ sorted(nbd, key = lambda x: min(x))[:k] for nbd in nb_dists])

        frac_diffs = [abs(diff(elt, 1).flatten()) / mean(elt,1) for  elt in nb_dists]
        abs_diffs = [abs(diff(elt, 1).flatten()) for  elt in nb_dists]
        
        ct = mycolors.getct(nq)
        f = myplots.fignum(4, (10,8))
        ax = f.add_axes([.05,.08,.25,.87])
        seismic.seismic(abs_diffs, ax = ax, colors = ct)
        
        jaccard = mean([float(len(nb_intersection[i])) / float(len(all_nbs[i])) for i in range(nq)])

        ax2 = f.add_axes([.34,.08,.6,.87])
        for i,d in enumerate(nb_dists):
            ax2.scatter(d[:,0], d[:,1], 20, alpha = .5,color =ct[i])

        
        lin = linregress(nb_dists[:,:,0].flatten(),nb_dists[:,:,1].flatten())
        rsquared = lin[2]**2

        ax2.annotate('NN dists for multi/struct-aligned trees.\nK = {0}'.format(k),
                    [0,1], xycoords = 'axes fraction', va = 'top',
                    xytext = [10,-10],textcoords = 'offset pixels')
        ax2.annotate('R-Squared: {0:3.3}\nJaccard Index: {1:3.3}'.format(rsquared, mean(jaccard)),
                    [1,0], xycoords = 'axes fraction', ha = 'right',
                    xytext = [-10,10],textcoords = 'offset pixels')
        ax2.set_xlabel('Muscle aligned tree distances')
        ax2.set_ylabel('Struct algined tree distances')
        
        datafile = cfg.dataPath('figs/gpm2/pt2_mus_cm_tree_dists_{0}_k={1}.svg'.format(run_id, k))
        f.savefig(datafile, dpi = 200, format = 'svg', bbox_inches = 'tight')
            
  
def draw_cm_muscle_congruencies(seqs, profiles, run_id, reset = True):
    print 'computing alignments...'
    print '  ...using muscle'
    malis, mrefs, mpairs =\
            mem.getOrSet(setAlignments, 
                         **mem.rc({},
                                  seqs = seqs, profiles = profiles, 
                                  run_id = run_id, ali_type = 'muscle',
                                  reset = reset,
                                  on_fail = 'compute', 
                                  register = 'tuali_musc_{0}'.format(run_id))) 
    print '  ...using cmalign.'
    salis, srefs, spairs  =\
        mem.getOrSet(setAlignments, 
                     **mem.rc({},
                              seqs = seqs, profiles = profiles, 
                              run_id = run_id, ali_type = 'struct',
                              reset = reset,
                              on_fail = 'compute', 
                              register = 'tuali__struct_{0}'.format(run_id)))
 
    print '  ...making trees.'
    
    for idx, alis in enumerate(zip(malis, salis)):
        m, s = alis
        mtree  = phyml.tree(m,run_id, bionj = True)
        stree  = phyml.tree(s,run_id, bionj = True)
        
        maps = dict([(elt.id,i) for i, elt in enumerate(m)])
        mdists = zeros((len(maps),len(maps)))
        sdists = zeros((len(maps),len(maps)))
        for n1 in mtree.get_terminals():
            for n2 in mtree.get_terminals():
                mdists[maps[n1.name],maps[n2.name]] = \
                    mtree.distance(n1,n2)
        
        for n1 in stree.get_terminals():
            for n2 in stree.get_terminals():
                sdists[maps[n1.name],maps[n2.name]] = \
                    stree.distance(n1,n2)
        tree_similarity(sdists, mdists, '{0}_struct_{1}'.format(run_id,idx), k = len(sdists - 1))
        tree_similarity(sdists, mdists, '{0}_struct_{1}'.format(run_id,idx), k = 6)

        f = myplots.fignum(4, (8,10))
        ct = mycolors.getct(len(mtree.get_terminals()))

        import networkx

        for t, sp, ttype in zip([mtree, stree], [211,212], ['sequence', 'structural']):
            a = f.add_subplot(sp)
            layout = 'neato'
            G = phylo.to_networkx(t)
            Gi = networkx.convert_node_labels_to_integers(G, discard_old_labels=False)
            posi = networkx.pygraphviz_layout(Gi, layout, args = '')
            posn = dict((n, posi[Gi.node_labels[n]]) for n in G)


            networkx.draw(G, posn, labels = dict([(n, '') for n in G.nodes()]),
                      node_size = [100 if  n.name in maps.keys() else 0 for n in G.nodes()],
                      width = 1, edge_color = 'black',
                      ax = a,
                      node_color = [ct[maps.get(n.name, -1)] for n in G.nodes()] )
        

            a.annotate('Embedded tree for {0} alignment.'.format(ttype),
                    [0,1], xycoords = 'axes fraction', va = 'top',
                    xytext = [10,0],textcoords = 'offset pixels')
            a.annotate('Total branch length is {0}'.format(t.total_branch_length()),
                    [1,0], xycoords = 'axes fraction', ha = 'right',
                    xytext = [-10,10],textcoords = 'offset pixels')            

        #phylo.draw_graphviz(  mtree,  label_func = lambda x: '', 
        #                      node_color = [ct[maps.get(n.name, -1)] for n in G.nodes()] +\
        #                          [ct[0] for n in mtree.get_nonterminals()], axes = ax)

        datafile = cfg.dataPath('figs/gpm2/pt2_mus_cm_tree_embeddings_{0}_struct_{1}.svg'.format(run_id, idx))
        f.savefig(datafile, dpi = 200, format = 'ps')
                


def run(rfid,run_id, inp_run_id, reset = True,
        draw_alis = draw_all_hard):

    sgs = get_seq_groups(rfid = rfid, **mem.sr({},reset = reset))
    all_seq_group_datas = []
    for s in sgs:
        all_seq_group_datas.append(eval_seq_group(s,rfid, run_id, inp_run_id, reset = reset,
                                                  draw_alis = draw_all_hard))
    return all_seq_group_datas
    
def eval_seq_group(gap_seqs, rfid, run_id, inp_run_id, reset = True,
                   draw_alis = draw_all_hard,
                   clade_alignment_method = clade_alignment_method,
                   max_structs = 10):

    rutils = utils
    data = butils.load_data(inp_run_id, 'output')
    structs = data['structs']
    energies = data['energies']
    esrt = argsort(energies)[::-1]
    s_inds = esrt[:max_structs]
    structs, energies = [structs[i] for i in s_inds], [energies[i] for i in s_inds]

    refseq = data['seq']
    
    nq = len(gap_seqs)
    ns = len(structs)

    names = ['N{1:04}'.format(rfid, idx) for idx in range(nq)]
    seqs = [rutils.ungapped_seq(gap_seqs[i], names[i]) for i in range(nq)]
    


    profiles = mem.getOrSet(setProfiles, 
                            **mem.rc({},
                                     seq = refseq, structs = structs, run_id = rfid,
                                     reset = reset,
                                     on_fail = 'compute', 
                                     register = 'tuprof_{0}'.format(rfid)))
    
    if draw_alis: 
        draw_cm_muscle_congruencies(seqs, profiles, 
                                    run_id, reset = reset)
    

    if clade_alignment_method == 'cm':
        alis, refs, all_pairs  =\
            mem.getOrSet(setAlignments, 
                         **mem.rc({},
                                  seqs = seqs, profiles = profiles, 
                                  run_id = rfid, ali_type = 'struct',
                                  reset = reset,
                                  on_fail = 'compute', 
                                  register = 'tuali_struct_{0}'.format(rfid)))
    else:
        raise Exception('No methods besides cm are yet implemented')
    

    seq_group_data = {}
    seq_group_data['seqs'] = gap_seqs
    seq_group_data['structs'] = []
    for i, struct in enumerate(structs):
        struct_data = {}
        ali = alis[i]
        ref = refs[i]
        pairs = all_pairs[i]
        
        #NOTE THAT DUE TO AN AWKWARD SYNTAX DECISION,
        #I AM ALLOWING FOR THE POSSIBILITY THAT EACH
        #ALI ELT HAS DIFFERENT PAIRS.
        #
        #ALL OF MY ROUTINES SO FAR ONLY USE A SINGLE 
        #PAIR SET AND SO I USE PAIRS[0] EXCLUSIVELY
        struct_data.update(ref = ref[0], 
                           pairs = pairs[0],
                           ali = ali)
                        
        rid = '{0}_{1}'.format(run_id, i)

        if clade_tree_method ==  'bionj': 
            tree = phyml.tree(ali, run_id = rid, bionj = True)
        else: tree = get_phase_tree(ali, pairs[0], run_id)

        for i, ct in enumerate(tree.get_terminals()):
            seq = filter(lambda x: x.id == ct.name, ali)[0]
            ct.m = {'seq':seq,
                    'probs':array([1 for j in range(len(seq))])}

        if clade_ancestor_method == 'independent':
            ml_tree = get_ml_ancestor_tree(tree, ali, 
                                           '{0}_paml{1}'.format(run_id, i))
        else:
            ml_tree = get_structure_ancestor_tree(\
                tree, ali,'{0}_stree{1}'.format(run_id, i))
        
        muts, times, gaps, irresolvables = tree_conservation.count_struct(ml_tree, pairs[0])

        struct_data.update(muts = muts, times = times, 
                        gaps = gaps, irresolvables = irresolvables)
        seq_group_data['structs'].append(struct_data)

    return seq_group_data

def get_structure_ancestor_tree(tree, ali, run_id):
    raise Exception('Get phase working!!! This does not do anything')

def get_ml_ancestor_tree(tree, ali, run_id):
    
    print 'Running ancestor inference in PAML'
    #RUN PAML
    rstfile= paml.run_paml(tree, ali, run_id = run_id)
    anc_tree = paml.rst_parser(rstfile) 
          
    #Label extent and internal nodes with sequences.
    for term in anc_tree.get_terminals():
        #COPY OLD TERMINAL NODES TO THE NEW TREE
        term.m = filter( lambda x: x.name == term.name, 
                         tree.get_terminals())[0].m

    for node in anc_tree.get_nonterminals():
        #REPLACE Ts WITH US IN THE TERMINALS OF THE NEW TREE
        node.m['seq'].seq = node.m['seq'].seq.replace('T', 'U')

    return anc_tree
              

def get_phase_tree(ali, struct, rid):
    '''Once completed, this routine will
return trees computed using a dependent
mutation model from phase. 
'''
    phase.make_mcmc(rid,
                    struct, 
                    [str(a.seq)  for a in ali], 
                    [a.name for a in ali])

        
    phase.make_ml(rid,
                  struct, 
                  [str(a.seq)  for a in ali], 
                  [a.name for a in ali])
    raise Exception('No methods besides bionj are yet implemented')


def setProfiles( seq = None, structs = None, run_id = None,
                 **kwargs ):
    '''
Make profiles for structs in a single sequence using 
infernal.
'''
    assert seq; assert structs; assert run_id
    profiles = infernal.profiles(seq, structs, run_id)
    return profiles

def setAlignments( seqs = None, profiles = None, run_id = None, ali_type = 'struct',
                   **kwargs ):
    '''
Make N structural alignemnts for each of N profiles
having the same M sequences apiece. 
'''
    assert seqs; assert run_id; assert profiles
    
    ns = len(profiles)
    nq = len(seqs)
    rutils = utils

    if ali_type == 'struct':
        ali_out = [infernal.alignment(seqs, p, run_id) for p in profiles]
        alignments = [a[0] for a in ali_out]
        refs = [[a[1]]*nq for a in ali_out]
        stks =  [a[2] for a in ali_out]
        pairs = [[rutils.stk_pairs(stk)]*nq for stk in stks] 
    elif ali_type == 'muscle':
        alignment = muscle.align(seqs)
        alignments = [alignment] * ns     
        
        raise Exception('MUSCLE ALIGNMENT NOT YET IMPLEMENTED... SORRY :(')
        #ONCE IMPLEMENTED, WILL FIND HOMOLOGOUS STRUCTURES FOR EACH
        #PROFILE IN THE GIVEN MULTIPLE SEQUENCE ALIGNMENT.

        #IT REMAINS UNCLEAR HOW I WILL PUT THESE INTO A TREE AND TRACK THE
        #STRUCTURAL ELEMENTS THROUGH IT!

        #all_refs = []
        #all_pairs = []
        #for i, profile  in enumerate(profiles):
        #    these_refs = []
        #    these_pairs = []
        #    for j, q in enumerate(seqs):
        #        ali_out = [infernal.alignment([q],profile, run_id)]
        #        these_refs.append(ali_out[0][1])
        #        stk = ali_out[0][2]
        #        
        #        this_ali = ali_out[0][0]
        #        this_seq = this_ali[0]
        #        all_seqs.append(this_seq)
        #        these_pairs.append( rutils.stk_pairs(stk))
        #        raise Exception()
        #    all_refs.append(these_refs)
        #    all_pairs.append(these_pairs)
        refs = None
        pairs = None
                
    else:  
        raise Exception('sorry not implemented')
    return alignments, refs, pairs
