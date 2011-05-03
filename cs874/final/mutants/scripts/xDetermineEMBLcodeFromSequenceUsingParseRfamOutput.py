#! /usr/bin/env python

# xDetermineEMBLcodeFromSequenceUsingParseRfamOutput.py
# determine Rfam ID for target sequence in faa file

import sys,os
from misc import _RT_
import math

def main(target,filename):
  file = open(filename)
  line = file.readline()
  while line:
    FASTA = line.strip()[1:]
    rna   = file.readline().strip()
    if rna == target:
      print FASTA
    else:
      line = file.readline()
      line = file.readline()
  file.close()

if __name__ == '__main__':
  if len(sys.argv) < 3:
    print "Usage: %s target filename" % sys.argv[0]
    sys.exit(1)
  target   = sys.argv[1]
  filename = sys.argv[2]
  main(target,filename)


