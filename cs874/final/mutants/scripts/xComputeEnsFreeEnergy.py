#! /usr/bin/env python

# xComputeEnsFreeEnergy.py
# P.Clote

import sys,os
from misc import _RT_
import math

def main(filename):
  file = open(filename)
  line = file.readline()
  while line:
    words = line.split()
    k     = int(words[0])
    x     = float(words[1])
    y     = - _RT_ * math.log(x)
    print "%d\t%f" % (k,y)
    line = file.readline()
  file.close()

if __name__ == '__main__':
  if len(sys.argv) < 2:
    print "Usage: %s filename" % sys.argv[0]
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)


