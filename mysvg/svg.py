#!/usr/bin/env python


def usage():
    print '''
Usage str [func_name] [input]

func_name:

  rs2pdb:   converts stdin rasmol file to a PDB file. 
  fix:   converts a broken pdb to a less broken pdb
  cent:     center a pdb file on a residue 
  rup:      rotate a pdb file so that a residue is on top.

'''
    return 0
  

if __name__ == '__main__':
    args = sys.argv[1:] if len(sys.argv) > 1 else exit(usage())
    prog = args[0]
                  
    if prog == 'rs2pdb':
        rs2pdb()
    elif prog == 'cent':
        struct_center()
    elif prog == 'rup':
        num = 1 if len(args) < 2 else int(args[1])
        residue_up(num)
    elif prog == 'fix':
        fix_pdb()
    exit(0)
