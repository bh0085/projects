import Bio.SeqIO as SeqIO
import Bio.Seq as seq
import Bio.Alphabet as Alphabet
import Bio.Align as al
import Bio.AlignIO as aio
import Bio.SeqRecord as sr

from numpy import *
import os, re
import cb.config as cfg
temp_dir = cfg.dataPath('alg')
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)
import subprocess as spc

def writelines(lines,gid):
    ali = al.MultipleSeqAlignment([sr.SeqRecord(seq.Seq(''.join(line)),
                                                '{1}{0}i'.format(i,gid),
                                                description = '{1}{0}d'.format(i,gid)
                                                ) 
                                   for i, line in enumerate(lines[::50])])
    return ali
    
def write_fa(lines, gid):
    fpath = os.path.join(temp_dir, '{0}_seqlines.fa'.format(gid))
    fopen = open(fpath,'w')
    alignment = writelines(lines,gid)
    SeqIO.write(alignment,fopen,'fasta')
    fopen.close()
    return fpath

def make_db(fpaths, prefix):
    '''makeblastdb -in good_proteins.fasta -dbtype prot -out my_prot_blast_db'''

    path = os.path.split(fpaths[0])[0]
    fnames = [ os.path.split(fp)[1] for fp in fpaths]
    cmd =  'formatdb -i "{0}" -p T -o T -n {1} '.format(' '.join(fnames),    prefix)
    prc = spc.Popen(cmd, 
                    shell = True,
                    cwd = path,
                    stdout = spc.PIPE)
    out =  prc.stdout.read()
    return out

def query_db(fpath, prefix):
    '''blastp -db my_prot_blast_db -query good_proteins.fasta 
    -outfmt 6 -out all-vs-all.tsv'''
    
    path,fname = os.path.split(fpath)
    #prefix = '.'.join(os.path.splitext(fname)[:-1])
    outfile =  'q{0}_d{1}.out'.format(fname, prefix)
    cmd = '''blast -p blastp -d {1} -i {0} -m 0 -o q{0}_d{1}.out'''\
        .format(fname,  prefix)

    prc = spc.Popen(cmd, 
                    shell = True,
                    cwd = path,
                    stdout = spc.PIPE)
    out = prc.stdout.read()
    print out
    return os.path.join(path,outfile)

def query_parse(outpaths):
    path_hits = {}
    for p in outpaths:
        
        fopen = open(p)
        content = fopen.readlines()
        fopen.close()
        nlines = len(content)
        qidxs = [i for i ,l in enumerate(content) if l[:7] == 'Query= ']
        qups = roll(qidxs,-1)
        qups[-1] = nlines
        hitlens = []
        
        hits = {}
        for i,q in enumerate(qidxs):
            rng = content[q:qups[i]]
            ptr = 0
            sig_matches = []

            qname = re.search(re.compile('Query= (.*)'), rng[0]).group(1)
            hits[qname] = {}
            while ptr <= len(rng)-1:
                ptr += 1;
                if ptr >= len(rng): break
                if rng[ptr][0:5] == 'Seque': break
            
            while ptr <= len(rng)-1:
                ptr += 1
                if rng[ptr][0] != '>': sig_matches.append(rng[ptr])
                else: break
            
            while ptr <= len(rng)-1:
                l = rng[ptr]
                title = l.strip()
                hit = {}
                hits[qname][title] = hit
                ptr += 3
                l = rng[ptr]
                bits = float(re.compile('Score =\s*(\d+)').search(l).group(1))

                expect = float(re.compile('Expect = ([\S]+),').search(l).group(1))
                hit['bits'] = bits
                hit['expect'] = expect
                while ptr < len(rng) -1:
                    l = rng[ptr]
                    ptr += 1
                    if rng[ptr][0] == '>': break
                if not ptr < len(rng) -1:
                    break
        path_hits[p] = hits
    return path_hits
