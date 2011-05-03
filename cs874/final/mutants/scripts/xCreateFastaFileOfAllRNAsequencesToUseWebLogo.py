#! /usr/bin/env python

# xCreateFastaFileOfAllRNAsequencesToUseWebLogo.py
# P.Clote

import sys,os
from misc import _RT_
import math

def main(filename):
  file = open(filename)
  line = file.readline()
  k    = 0
  while line:
    if line[0]=='>':
      line = file.readline()
      continue
    k     += 1
    rna    = line.strip().upper()
    secStr = file.readline().strip()
    print "> %d" % k
    print rna
    line = file.readline()
  file.close()

if __name__ == '__main__':
  if len(sys.argv) < 2:
    print "Usage: %s filename" % sys.argv[0]
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)


