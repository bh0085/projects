import networkx as nx
from numpy import *
import numpy as np
import subprocess as spc
import compbio.config as cfg
import os
import modularity, cs2, aba_info as aba
import cb.utils.graphs.draw as gd
import cb.utils.colors as mycolors
import nx_addons.community as community

def edge_weight_fun():
    min_edge = 5
    edge_randomness = 5
    return int(floor(random.rand() * edge_randomness) + min_edge)

def make_sf(n = 200,
            alpha = .41,
            beta = .54,
            gamma = .05):
    g0 = nx.scale_free_graph(n = n,
                             alpha = alpha,
                             beta = beta,
                             gamma = gamma)
    g = nx.DiGraph()
    g.add_weighted_edges_from([[e[0],e[1],edge_weight_fun()]
                                for e in g0.edges()])
    return g

def make_sf_mesh(n =5,
                 sf_params = None,
                 mesh_params = None):
    '''
    For default values of sf_params, sf_mesh see flow_all.
'''
    
    cxn_rule =mesh_params['cxn_rule']
    if cxn_rule == 'random_matrix':
        block_odds = reshape(random.rand(square(n)), (n,n))*.1
        block_counts = [sf_params['n'] for i in range(n)]
    elif cxn_rule == 'mouse_brain':
        centers,volumes = aba.brain_regions(n)
        dists = sum(square(centers),1)[:,newaxis] \
            + sum(square(centers),1)[newaxis,:] \
            - 2 * dot(centers,centers.T)
        ranks = argsort(argsort(dists,1),1)
        block_odds = 1.0 / square((ranks +1))
        for i in range(len(block_odds)): block_odds[i,i] = 0
        block_counts = array([v  * sqrt(sf_params['n']/ mean(volumes) )
                              for v in volumes],int)    
    else:
        raise Exception('connection rule {0} not understood'.format(cxn_rule))
    
    sf_graphs = []
    for i in range(n):
        sfp ={}
        sfp.update(sf_params)
        sfp['n'] = block_counts[i]
        g = make_sf(**sfp)
        sf_graphs.append(g)
        
    ofs= 0
    gbig = nx.DiGraph()
    communities = []
    for i,g in enumerate(sf_graphs):
        nodes = [n + ofs for n in g.nodes()]
        edges = [(e[0] + ofs, e[1] + ofs , g.get_edge_data(e[0],e[1])['weight'])
                 for e in g.edges()]
        communities.append(nodes)
        gbig.add_nodes_from(nodes)
        gbig.add_weighted_edges_from(edges)
        ofs += len(g)
        

    for i,ci in enumerate(communities):
       for j,cj in enumerate(communities):
           in_degs = gbig.in_degree(ci)
           for node_name in ci:
               indeg = in_degs[node_name]
               out_idxs = np.unique(array(
                       np.floor(random.rand(
                               indeg*block_odds[i,j])*len(cj))
                       ,int))
               out_nodes = [cj[oidx] for oidx in out_idxs]
               out_edges = [(node_name, onode, edge_weight_fun()) 
                            for onode in out_nodes ]
               gbig.add_weighted_edges_from(out_edges)
               
    return gbig, communities

def noisy(g0,n= 2, efrac = .8):
    graphs = []
    for i in range(n):
        g = nx.DiGraph()
        g.add_nodes_from(g0)
        ge = list(g0.edges_iter())
        rnd_edges = [(ge[i][0], ge[i][1], g0[ge[i][0]][ge[i][1]]['weight']) 
                     for i in random.permutation(len(ge))\
                         [:len(ge) * efrac]]
        g.add_weighted_edges_from(rnd_edges)
        graphs.append(g)
        
    return graphs


default_sf_params = dict(alpha =.2, 
                         beta = .75,
                         gamma = .05,
                         n = 300)

default_mesh_params = dict(cxn_rule = 'mouse_brain',
                           modules = 10)

def flow_all(name = 'sf_mesh',
             modules = 10,
             sf_params = default_sf_params,
             mesh_params = default_mesh_params):
    communities = None
    if name == 'sf_small':
        g = make_sf(**sf_params)
    elif name == 'sf_mesh':
        g,communities = make_sf_mesh(n = modules, 
                         sf_params = sf_params,
                         mesh_params = mesh_params)
        
    #cs2.run_flow(g,1)
    ngs = noisy(g, 4)
    fgs = []
    for i,ng in enumerate(ngs):
        fg = cs2.run_flow(ng, i)
        fgs.append(fg)
        
    return g, ngs, fgs, communities


def show_brain_flows(g,ngs,fgs, communities):

    centers,volumes, voxels = aba.brain_regions(len(communities),
                                                return_voxels = True)
    p0 = {}
    n_communities = {}
    for i, c in enumerate(communities):
        nv = len(voxels[i])
        for j, n in enumerate(c):
            p0[n] = voxels[i][int(floor(random.rand()*nv))][:2] + random.rand()*.4
            n_communities[n] = i
        
    comm_ct = mycolors.getct(len(n_communities))
    nodelist = g.nodes()
    
    ckw = dict([(k,dict(facecolor = 'gray',
                        alpha = 1,
                        linewidth = .5,
                        arrowstyle = '-|>',
                        edgecolor = 'black',
                        color = comm_ct[n_communities[k[0]]],
                        shrinkA = 0,
                        shrinkB = 0))
                for k in g.edges() ])

    skw =dict(facecolor = 'none',
              edgecolor = [comm_ct[n_communities[n]] 
                           for n in nodelist],
              s = 1)

                      
    gd.draw(g,p0,g.edges()[::10], 
            scatter_nodes = nodelist,
            ckw = ckw,
            skw = skw,
            ckalpha = .8,
            cktype = 'simple')
    
    
    return
    

    colors = mycolors.getct(len(fgs))
    for i,fg in enumerate(fgs):
        edges = fg.edges()
        weights = [fg[e[0]][e[1]]['weight'] for e in edges]
        ckw = dict([(k, dict(color = colors[i],
                             linewidth = weights[j]/3))
                    for j, k in enumerate(edges)
                    ])
                            
        gd.draw(fg, 
                p0, 
                edges,
                ckw = ckw,
                scatter_nodes = [],
                cktype = 'simple',
                ckalpha = .25,
                )
        
                
                
            


        


def show_flows(g, ngs, fgs):

    p0 = gd.getpos(g)

    ckw = dict([(k,dict(facecolor = 'gray',
                        alpha = 1,
                        linewidth = .5,
                        arrowstyle = '-|>',
                        edgecolor = 'black',
                        color = 'gray',
                        shrinkA = 0,
                        shrinkB = 0))
                for k in g.edges() ])
    skw =dict(facecolor = 'none',
              edgecolor = 'black',
              s = 20)
                      
    gd.draw(g,p0,g.edges(), 
            ckw = ckw,
            skw = skw,
            cktype = 'simple')
    
    
    colors = mycolors.getct(len(fgs))
    for i,fg in enumerate(fgs):
        edges = fg.edges()
        weights = [fg[e[0]][e[1]]['weight'] for e in edges]
        ckw = dict([(k, dict(color = colors[i],
                             linewidth = weights[j]/3))
                    for j, k in enumerate(edges)
                    ])
                            
        gd.draw(fg, 
                p0, 
                edges,
                ckw = ckw,
                scatter_nodes = [],
                cktype = 'simple')
        
                
                
            

