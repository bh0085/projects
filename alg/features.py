from numpy import *
def compute_cons_features(g0,
                          ftype = 'id'):
    feats = {}
    nlist = g0.nodes()
    if ftype == 'id':
        for n in nlist:
            feats[n] = n

    if ftype in ['nt', 'aa8', 'aa20']:
        if ftype == 'nt':
            n_basemap =dict([(n, ['A','T','G','C'][int(floor(random.rand()*4))]) 
                             for i,n in enumerate(nlist)])
        elif ftype == 'aa8':
            n_basemap =dict([(n, list('RHKDESTN')[int(floor(random.rand()*8))]) 
                             for i,n in enumerate(nlist)])
        elif ftype =='aa20':
            n_basemap =dict([(n, list('RHKDESTNQCGPAVILMFYW')\
                                  [int(floor(random.rand()*20))]) 
                             for i,n in enumerate(nlist)])

        for n in nlist:
            feats[n] = n_basemap[n]

    return feats
    

def assign_cons_features(g0_feats,g1):
    '''
assign features to a graph g1 derived from a consensus g0.

ftype == 'id': the features are computed simply as the id
               of a given node in the original graph.
'''

    feats ={}
    for n in g1.nodes():
        feats[n] = g0_feats[n]
    return feats
