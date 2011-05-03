#! /usr/bin/env python

# xCreateInputForGraphDotPlotForAllPtwiseMutantsAllK.py
# P.Clote

import sys,os
from misc import _RT_
import math

def main(filename):
  file = open(filename)
  file.readline() #discard line
  rna = file.readline().strip().upper()
  secStr = file.readline().strip()
  print rna
  print secStr
  line = file.readline()
  while line:
    if line[0]=='>':
      line = file.readline()
      continue
    rna    = line.strip()
    secStr = file.readline().strip()
    print secStr
    line = file.readline()
  file.close()

if __name__ == '__main__':
  if len(sys.argv) < 2:
    print "Usage: %s filename" % sys.argv[0]
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)

