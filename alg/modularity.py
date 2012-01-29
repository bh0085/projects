import nx_addons.community as community
import networkx as nx


def partition(g):
    return community.best_partition(g)

def modularity(g, partition):
    return community.modularity(partition, g) 
