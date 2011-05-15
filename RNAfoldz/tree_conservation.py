from numpy import *
import numpy as np
import itertools as it

bps = ['GC','AU','GU','CG','UA','UG']


def create_parents(tree):
    '''Create a parent attribute for each node in
tree.

returns: None
'''
    #CREATE A PARENTS ATTRIBUTE TO USE FOR 
    # SCANNING PAIRS.
    for n in tree.get_nonterminals():
       if n == tree.root: n.parent = None
       else: n.parent = ([tree.root] + tree.root.get_path(n))[-2] 
    for n in tree.get_terminals():
       n.parent =([tree.root] +  tree.root.get_path(n))[-2]      

def get_preterminals(curterms):
    '''Get a list of all nodes with exactly two terminals as children.

returns: Preterminal nodes.
'''
    ct_list = sorted(list(curterms), 
                     key = lambda x: x.parent)
    preterms =[(k,list(g))
               for k,g in it.groupby(ct_list,lambda x: x.parent)]
    preterms =filter(lambda x: len(x[1]) == 2,preterms)
    return preterms

def do_contract(curterms, curints,terms,pp):
    '''Contract a parental node, returning the branch lengths 
of edges connecting it to its children, 'terms' 

returns: Branch lengths'''
    curterms.difference_update(terms)
    curints.difference_update([pp])
    curterms.update([pp])
    return array([t.branch_length for t in terms])

def do_giveup(curterms,curints, terms, pp):
    '''Give up contraction of the current line of descent.
Simply remove the parent the the set of eligible internals
and remove the terminals from the set of eligible terms.

returns: None'''
    curterms.difference_update(terms)
    curints.difference_update([pp])

def compute_contraction(idx, pbases, tbases, 
                        curterms, curints,terms,  pp,
                        times, muts):

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
                   muts['wob'][idx] += 1
               else:
                   muts['comp'][idx] +=1
               times['paired'][idx] += s[j];
           else:
               if sum(diffs)==0:
                   pass #never reached
               elif sum(diffs)==1:
                   muts['ucom'][idx] +=1
               else:
                   muts['ucom'][idx] +=1 #SHOULD BE 2?
               times['unpaired'][idx] += s[j]
       else:
           if cpar:
               if sum(diffs)==0:
                   pass #never reached.
               elif sum(diffs) ==1:
                   muts['reco'][idx] +=1
               else:
                   muts['reco'][idx] +=1 #SHOULD BE 2?
               times['paired'][idx] += s[j]
           else:
               if sum(diffs) == 0:
                   pass
               elif sum(diffs) == 1:
                   muts['bbad'][idx] += 1
               else:
                   muts['bbad'][idx] += 1 #SHOULD BE 2?
               times['unpaired'][idx] += s[j]
   times['total'][idx] += sum(s)
         

def count_struct(tree, pairs):
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
    
   create_parents(tree)
   #DATA MATRICES COUNTING SUB TYPES (INTS)
   m_wob, m_ucom, m_comp, m_reco, m_bbad = \
       [zeros(len(pairs)) for i in range(5)] 
   #DATA MATRICES COMPUTING PAIRING OVER TIME (FLOATS)
   t_up, t_p, t_tots = \
       [zeros(len(pairs)) for  i in range(3)]

   #ALIGNMENT GAP COUNTS PER PAIR
   gaps = zeros(len(pairs))
   #IRRESOLVABLE TERMINAL NODES
   irresolvables = [[] for i in range(len(pairs))]

   #DICTS HOLDING SUBS/TIMES
   muts = dict(wob = m_wob, ucom = m_ucom, comp = m_comp, 
               reco = m_reco, bbad = m_bbad)
   times= dict(unpaired = t_up, paired = t_p,
               total = t_tots,
               total_incl_unresolved = tree.total_branch_length() -\
                   tree.root.branch_length)

   for idx,p in enumerate(pairs):
       #DEAL WITH GAPS (SKIP)
       if '-' in list(it.chain(*[(t.m['seq'].seq[p[0]], \
                                      t.m['seq'].seq[p[1]])
                                 for t in tree.get_terminals()])):
	   gaps[idx] = 1
           continue
       
       curterms= set(tree.get_terminals())
       curints = set(tree.get_nonterminals())
       while 1:
           #GET PRETERMINALS, IF NONE THEN BREAK
           preterms = get_preterminals(curterms)
           if not preterms: break

           for pp, terms in preterms:
               tprobs = [(t.m['probs'][p[0]], t.m['probs'][p[1]]) 
                          for t in terms]
               tbases = [(t.m['seq'][p[0]], t.m['seq'][p[1]])
                          for t in terms]

               pprobs = pp.m['probs'][p[0]],pp.m['probs'][p[1]]
               pbases= pp.m['seq'][p[0]],pp.m['seq'][p[1]]

               if  sum([1.-prob for prob in pprobs]) > .30:
                   irresolvables[i].append(terms)
                   do_giveup(curterms,curints,terms,pp)   
               else:
                   compute_contraction(idx,pbases, tbases, 
                                       curterms, curints, terms,pp,
                                       times, muts)

   return muts, times, gaps, irresolvables
