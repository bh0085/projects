import cb.utils.plots as myplots
from numpy import *
def plot_mers(mer_cts):
    f = myplots.fignum(3, (8,8))
    
    
    ax = f.add_subplot(111)
    hist,bin_edges = histogram(mer_cts.values(), 20)
    ax.fill_between(bin_edges[:-1], log(hist),
                    edgecolor = 'black',
                    linewidth = 5)
    ax.set_xlabel('mer rediscovery rate')
    ax.set_ylabel('$log(n)$')
    ax.set_title('Frequencies of 5-mer reoccurence across 10,000 walks.')

    f.savefig(myplots.figpath('mer_ct_hist'))
    
    return 
import numpy as np
from numpy import *
import itertools as it


def align_heatmap(parsed):
    p0 = parsed.values()[0]
    bitlens = array( [sorted([e['expect'] for e in val.values()])[:-1] 
                      for val in p0.values()])
    f = myplots.fignum(3, (8,8))

    ax = f.add_subplot(111)
    nodes = set(p0.keys())
    for v in p0.values():
        nodes = nodes.union(set(v.keys()))
        
    nmap = dict([(i, k) for i,k in enumerate(nodes)])
    r_nmap = dict([(k,i) for i,k in nmap.iteritems()])

    z = zeros((len(nodes),len(nodes)))
    for k,v in p0.iteritems():
        i = r_nmap[k]
        for k2,v2 in v.iteritems():
            j = r_nmap[k2]
            z[i,j] = 1 / (.0001 + v2['expect'])
    ax.imshow(z[argsort(sum(z,1)),:][::-1,::-1][:100,:100])

def align_len_histogram(parsed):

    p0 = parsed.values()[0]
    bitlens = array( [sorted([e['bits'] for e in val.values()])[:-1] 
                      for val in p0.values()]).flatten()
    bitlens = array(list(it.chain(*bitlens)))
    bitlens = 2 * (bitlens - np.min(bitlens.flatten())) + 8
    
    mind=  8 #min(deg_c.values())+.00001
    maxd = max(bitlens) #max(deg_c.values())/3
    bins = linspace(mind,maxd,8)
    
    h_paths,bin_edges = histogram(bitlens,bins)

    h_paths = array(h_paths,float)
    h_paths/= sum(h_paths)

    f = myplots.fignum(3, (8,8))

    ax = f.add_subplot(111)

    ax.plot(bins[:-1], h_paths, color = 'red')
    ax.set_xlabel('alignment hit length')
    ax.set_ylabel('frequency')
    ax.set_title('best matched substring lengths')

    raise Exception()

    f.savefig(myplots.figpath('walk_centrality'))
     

    paths_cat = paths.flat
    n = len(paths_cat)
    
    degs = [deg_c[p]for p in paths_cat[::10]]
    
    mind=  min(deg_c.values())+.00001
    maxd = max(deg_c.values())/3
    bins = linspace(mind,maxd,8)
    
    h_paths,bin_edges = histogram(degs,bins)
    h_rand,bin_edges  = histogram(deg_c.values(), bins)

    h_paths = array(h_paths,float)
    h_rand = array(h_rand,float)
    h_paths/= sum(h_paths)
    h_rand/= sum(h_rand)

    f = myplots.fignum(3, (8,8))

    ax = f.add_subplot(111)

    ax.plot(bins[:-1], h_paths, color = 'red')
    ax.plot(bins[:-1], h_rand, color = 'black')
    ax.set_xlabel('node centrality')
    ax.set_ylabel('frequency')
    ax.set_title('distribution of centrality in walks vs. random')
