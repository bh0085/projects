from numpy import *

def compute_occurences(strings, 
                       nl = inf, np = inf,
                       mer = 5):
    np = min(np, len(strings))
    nl = min(nl, len(strings[0]) - mer)
    
    mer_cts = {}
    for i in range(np):
        for j in range(nl): 
            subs = tuple(strings[i,j:j+mer])
            mer_cts[subs] = mer_cts.get(subs, 0) + 1
    
    return mer_cts


