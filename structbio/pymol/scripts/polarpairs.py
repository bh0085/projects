'''
(c) 2011 Thomas Holder, MPI for Developmental Biology
'''
 
from pymol import cmd
 
def polarpairs(sel1, sel2, cutoff=4.0, angle=63.0, name='', state=1, quiet=1):
    '''
ARGUMENTS
 
    sel1, sel2 = string: atom selections
 
    cutoff = float: distance cutoff
 
    angle = float: h-bond angle cutoff in degrees. If angle="default", take
    "h_bond_max_angle" setting. If angle=0, do not detect h-bonding.
 
    name = string: If given, also create a distance object for visual representation
 
SEE ALSO
 
    cmd.find_pairs, cmd.distance
    '''
    cutoff = float(cutoff)
    quiet = int(quiet)
    state = int(state)
    if angle == 'default':
        angle = cmd.get('h_bond_max_angle', cmd.get_object_list(sel1)[0])
    angle = float(angle)
    mode = 1 if angle > 0 else 0
    x = cmd.find_pairs('(%s) and donors' % sel1, '(%s) and acceptors' % sel2,
            state, state,
            cutoff=cutoff, mode=mode, angle=angle) + \
        cmd.find_pairs('(%s) and acceptors' % sel1, '(%s) and donors' % sel2,
            state, state,
            cutoff=cutoff, mode=mode, angle=angle)
    x = sorted(set(x))
    if not quiet:
        print 'Settings: cutoff=%.1fangstrom angle=%.1fdegree' % (cutoff, angle)
        print 'Found %d polar contacts' % (len(x))
    if len(name) > 0:
        for p in x:
            cmd.distance(name, '(%s`%s)' % p[0], '(%s`%s)' % p[1])
    return x
 
cmd.extend('polarpairs', polarpairs)
