import cb.utils.memo as mem
import mincost, blast
from numpy import *
import networkx as nx
import os
import features as feat

def get_flows(**kwargs):
    def set_flows(**kwargs):
        return mincost.flow_all()
    return mem.getOrSet(set_flows,**kwargs)

def get_features(g,ngs,fgs,ftype = 'id'):
    '''compute a set of features on the consensus graph and
    assign them appropriately to nodes in the others'''

    g0_feats = feat.compute_cons_features(g, ftype = ftype)
    nfeats = []
    ffeats = []
    for ng in ngs:
        nfeats.append(feat.assign_cons_features(g0_feats,ng))
    for fg in fgs:
        ffeats.append(feat.assign_cons_features(g0_feats,fg))
    return g0_feats, nfeats, ffeats
        


def get_feats(ftype = 'id'):
    g, ngs, fgs, communities = get_flows()
    return g0_feats, nfeats,ffeats
    
def get_paths(graph, nodelist):
    
    fg1 = nx.DiGraph()
    fg1.add_nodes_from(nodelist)
    fg1.add_weighted_edges_from([
            [e[0], e[1], graph[e[0]][e[1]]['weight']]
            for e in graph.edges()])

    paths = walk0(fg1)
    return paths

def feature_paths(paths, feats):
    ftype = type(feats.values()[0])
    strings = empty(shape(paths), ftype)
    for i in range(shape(paths)[0]):
        for j in range(shape(paths)[1]):
            strings[i,j] = feats[paths[i][j]]
    return strings
    
def write_walks(g, ngs, fgs, ftype = 'id'):
    g0_feats, nfeats, ffeats = get_features(g, ngs, fgs, ftype = ftype)
    
    outs = []
    nruns = 4
    for i, fg in enumerate(fgs):
        feats = ffeats[i]
        fpaths = []
        for j in range(nruns):
            print 'run ({0},{1})'.format(i, j)
            
            paths = get_paths(fg, g.nodes())
            return paths
            strings = feature_paths(paths, feats)
            gid = '{0}_{1}'.format(i,j)
            fpath = blast.write_fa(strings, gid)

            #fpath = os.path.join(blast.temp_dir, '{0}_seqlines.fa'.format(gid))
            fpaths.append(fpath)

        prefix = 'randwalk{0}'.format(i)
        out =blast.make_db(fpaths,prefix)
        outs.append((fpaths,prefix))

    return outs

def meta_0():
    g, ngs, fgs, communities = get_flows()
    return write_walks(g, ngs, fgs)

def meta_0q(outs):
    outfiles = []
    for i in outs:
        prefix = i[1]
        for j in i[0]:
            fpath = j
            print prefix, fpath
            outfiles.append(blast.query_db(fpath, prefix))
    return outfiles

def meta_0p(outfiles):
    blast.query_parse(outfiles)

def walk0(fg1, nlist):

    for e in fg1.edges():
        fg1[e[0]][e[1]]['flow_remaining'] = \
            fg1[e[0]][e[1]]['weight']        
    cont = True
    
    nwalks = 10000
    walk_len =500
    
    paths = zeros((nwalks, walk_len))
    np = nwalks
    ns = walk_len

    #nlist = fg1.nodes()
    print 'building matrices'
    adj_matrix = zeros((len(nlist),len(nlist)))
    for i,e1 in enumerate(nlist):
        for k,e2 in fg1[e1].iteritems():
            adj_matrix[i,k] = e2['flow_remaining']
    nrm = sum(adj_matrix,1)[:,newaxis]
    nrm += equal(nrm,0)
    transition_matrix = adj_matrix/nrm
    transition_csum = cumsum(transition_matrix,1)

    degs = sum(adj_matrix,1)
    cum_degs = cumsum(degs)
    stoch_degs = cum_degs / cum_degs[-1]
    print 'computed transition matrices'

    for i in range(np):
        current = searchsorted(stoch_degs, random.rand())
        for j in range(ns):
            next_node = searchsorted(transition_csum[current],random.rand())
            paths[i,j] = current
            current = next_node
    return paths

