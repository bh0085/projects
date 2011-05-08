import subprocess as spc, sys, re, os,inspect
import compbio.config as cfg
import Bio.Align as ba
import Bio
import utils as rutils


def alignment(seqs, profile,run_id):
    '''Compute an alignment of multiple sequences to a given 
covariance model profile such as constructed by cmbuild
via infernal.profiles.

input:
  seqs:    a list of biopython SeqRecord objects
  profile: the filename of a covariance model profile
  run_id:  a run id to use for naming temporary files to avoid collisions

output:
  ali:     an rfam multiple sequence alignment
  ref:     the profile reference sequence aligned to ali
  struct:  the profile reference structure aligned to ali

'''
    seqs = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''.join([let for let in str(ali.seq)  if let in 'AUTGC' ]),
                                                Bio.Seq.Alphabet.RNAAlphabet),
                                    'S{0:03}'.format(idx))
            for idx, ali in enumerate(seqs)]
    infile = cfg.dataPath('infernal/temp/{0}_{1:03}_unaligned.fa'.format(run_id,idx))
    outfile= cfg.dataPath('infernal/temp/{0}_{1:03}_aligned.stk'.format(run_id,idx))
    Bio.SeqIO.write(seqs, infile, 'fasta')
    
    cstr = 'cmalign -o {0} {1} {2}'.format(outfile, profile, infile)
    ispc = spc.Popen(cstr, shell = True,
                            stdout = spc.PIPE)
    out = ispc.communicate()[0]
    fopen = open(outfile)
    seqs, ref, struct = rutils.stk_parse(fopen)
    fopen.close()
    ali = ba.MultipleSeqAlignment(seqs)

    for a in ali:
	    a.seq = a.seq.upper()
    return ali, ref, struct
    


def profiles(seq, structs, run_id):
    '''Compute a sequence profile using cmbuild with --rsearch
from a single sequence and fixed secondary structure. 

The reason to call profiles for several structures at
once is to avoid filename collisions by automatically
generating filenames for each of n structs.

input:
  seq:      a biopython SeqRecord object.
  structs:  an array of biopython.
  run_id:   a run id to avoid collisions of temporary files.

output:
  profiles: paths to files containing cm profiles for each struct

'''
    exemplar_stks = []
    for i, s in enumerate(structs):
        stk = ['.'] * len(seq)
        for p in s: stk[p[0]] ,stk[p[1]] = '(',')'
        stk = ''.join(stk)
        exemplar_stks.append(rutils.stk_format(seq,stk))
    profiles = []
    for idx, stktext in enumerate(exemplar_stks):
        stkfile = cfg.dataPath('infernal/temp/{0}_{1:03}_{2}.stk'.format(seq.id, idx,run_id))
        cmfile = cfg.dataPath('infernal/temp/{0}_{1:03}_{2}.cm'.format(seq.id, idx,run_id))
        fopen = open(stkfile,'w')
        fopen.write(stktext)
        fopen.close()
        cstr = 'cmbuild -F --rsearch {0} {1} {2}'.format(cfg.dataPath('infernal/matrices/RIBOSUM85-60.mat'),cmfile,stkfile)
        ispc = spc.Popen(cstr, shell = True, stdout = spc.PIPE)
        out = ispc.communicate()[0]
        profiles.append(cmfile)
    
    return profiles
        

