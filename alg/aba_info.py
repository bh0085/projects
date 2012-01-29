import networkx as nx
import cb.config as cfg
import os, itertools as it
from numpy import *
import numpy as np
import cb.utils.colors as mycolors

aba_path = cfg.dataPath('brain_atlas/mouse_aba')
structures_path = os.path.join(aba_path, 'brainstructures.csv')
coords_path = os.path.join(aba_path, 'AtlasAnnotation200.sva')

def structure_voxels():
    structs_open = open(structures_path)
    slines = structs_open.readlines()
    structs_open.close()
    
    coords_open = open(coords_path)
    clines = coords_open.readlines()
    coords_open.close()
    
    crows = [[int(e) for e in  l.strip().split(',')] for l in clines[2:]]

    scols = slines[0].strip().split(',')
    srows = [dict( zip(scols,l.strip().split(','))) for l in slines[1:]]
    
    cgroups =[(k,list(g)) for k, g in  it.groupby(
        sorted(crows, key =lambda x: x[3])
        , key = lambda x: x[3])]
    
    coords = dict([(k, array([[e[0],e[1],e[2]] 
                    for e in v]))
                   for k,v in cgroups])
    
    return coords

def show_brain():
    n = 10
    import matplotlib.pyplot as plt
    f = plt.figure(3)
    plt.clf()
    ax = f.add_subplot(111)

    coords = structure_voxels()
    struct_keys =[e[0] for e in sorted( coords.iteritems(),
                                        key = lambda x:len(x[1]))[::-1][:n]
                  ]
    coords = dict([(k,coords[k]) for k in struct_keys])
    
    f = plt.figure(3)
    
    
    ct = dict([(k,v) 
               for k,v in zip(coords.keys(), mycolors.getct(len(coords)))])
    for k, v in coords.iteritems():
        ax.scatter(*v[::10,[0,1]].T,s = 50,c=ct[k],edgecolor = 'none')
        
    return

def brain_regions(n, return_voxels = False):
    coords = structure_voxels()
    volumes = []
    centers = []
    struct_keys =[e[0] for e in sorted( coords.iteritems(),
                                        key = lambda x:len(x[1]))[::-1][:n]
                  ]
    voxels = []
    for k in struct_keys:
        c = coords[k]
        voxels.append(array(c))        
        volumes.append(len(c))
        centers.append(np.mean(c,0))

    volumes = array(volumes)
    centers = array(centers)
    if not return_voxels:
        return centers, volumes
    else:
        return centers, volumes, voxels
