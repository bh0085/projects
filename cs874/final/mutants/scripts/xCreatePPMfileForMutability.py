#! /usr/bin/env python

# xCreatePPMfileForMutability.py
# P.Clote

import sys,os
from misc import _RT_
import math

#output ppm either in red/green or shades of gray

def main(filename):
  file = open(filename)
  rna0 = file.readline().strip()
  D    = {}; n = len(rna0); num = 0
  line = file.readline()
  numMut = 0
  while line:
    if line[0]=='>':
      words  = line.split()
      numMut = int(words[-2])
      D[numMut] = {}
      for i in range(n): D[numMut][i] = 0.0
         #D[numMut][i] is number of mutations in position i
      num  = 0 #start over counter of number of mutations
      line = file.readline()
      continue
    num    += 1
    rna    = line.strip().upper()
    secStr = file.readline().strip()
    for i in range(n):
      if rna0[i]!=rna[i]: D[numMut][i] += 1.0
    line = file.readline()
  file.close()
  X = n; Y = numMut
  if COLOR:
    print "P3"
    print X,Y
    numGrayShades = 255
    red = 0; green = 0; blue = 0
    print numGrayShades
    for x in range(numMut,0,-1):
      for i in range(n):
        D[x][i] = D[x][i]/num #mutability on scale of [0,1]
        red    = int(round((1-D[x][i])*numGrayShades,0)) #red means not mutable
        green  = int(round(D[x][i]*numGrayShades,0))     #green means mutable
        sys.stdout.write("%d %d %d " % (red,green,blue))
      print
  else: #no color, but rather black and white    
    print "P2"
    numGrayShades = 100
    for x in range(numMut,0,-1):
      for i in range(n):
        D[x][i] = D[x][i]/num #mutability on scale of [0,1]
        sys.stdout.write("%d " %int(round(D[x][i]*100,0)))
      print
  return D   

 
if __name__ == '__main__':
  if len(sys.argv) < 2:
    text = """Usage: %s filename [color (0/1)]
    1) file contains first line of RNA sequence (unmutated), then
       all subsequent lines contain RNAmutants samples, where each sample
       is given by RNA seq and RNA sec str, each on separate lines.
    2) color=1 (red/green), color=0 (gray scale), default is 1 """
    text = text % sys.argv[0] 
    print text
    sys.exit(1)
  filename = sys.argv[1]
  if len(sys.argv)>2:
    COLOR = int(sys.argv[2])
  else:
    COLOR = 1
  main(filename)


