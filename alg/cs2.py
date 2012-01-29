import networkx as nx
from numpy import *
import subprocess as spc
import compbio.config as cfg
import os, re

default_flow_dir = cfg.dataPath('graph_flows')
def run_flow(g,gid):
    if not os.path.isdir(default_flow_dir):
        os.mkdir(default_flow_dir)
    write_flow(default_flow_dir, g, gid)
    compute_flow(default_flow_dir, gid)
    return parse_flow(default_flow_dir, gid)

def flow_inp_file(flow_dir, gid):
    return os.path.join(flow_dir, 'flow_{0}.inp'.format(gid))
def flow_out_file(flow_dir, gid):
    return os.path.join(flow_dir, 'flow_{0}.out'.format(gid))
def write_flow(flow_dir,g,gid):
    nn = len(g.nodes())
    ne = len(g.edges())
    lines = []
    lines.append('p min {0} {1}'.format(nn, ne))
    lines.append('n 1 0')
    for e in g.edges():
        w = g[e[0]][e[1]]['weight']
        lines.append('a {0} {1} {2} {3} {4}'
                     .format(e[0] + 1, e[1] + 1,
                             0,int(w), -1 ))
    
    fin = flow_inp_file(flow_dir, gid)
    fopen = open(fin,'w')
    fopen.write('\n'.join(lines))
    fopen.close()
def compute_flow(flow_dir,gid):
    fin = flow_inp_file(flow_dir,gid)
    fout = flow_out_file(flow_dir,gid)
    cmd = 'cat {0} | cs2 > {1}'.format(fin, fout)
    prc = spc.Popen(cmd, shell = True)
    out = prc.communicate()
def parse_flow(flow_dir, gid):
    fout = flow_out_file(flow_dir, gid)
    fopen = open(fout)
    g = nx.DiGraph()
    flow_lines = [ re.compile('\s+').split(l.strip()) for l in fopen.readlines() if l[0] == 'f']
    fopen.close()

    flow_edges = [ (int(l[1])-1, int(l[2])-1, int(l[3])) for l in flow_lines]
    g.add_weighted_edges_from(flow_edges)
    return g

    
