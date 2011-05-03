#! /usr/bin/env python

# xSuperoptimalFreeEnergy.py
# P.Clote

import sys,os
from misc import _RT_
import math

def main(filename):
  file = open(filename)
  line = file.readline()
  k    = 0
  while line:
    words = line.split()
    x     = float(words[-1])
    print "%d\t%f" % (k,x)
    line = file.readline()
    line = file.readline()
    line = file.readline()
    k     += 1
  file.close()

if __name__ == '__main__':
  if len(sys.argv) < 2:
    print "Usage: %s filename" % sys.argv[0]
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)


