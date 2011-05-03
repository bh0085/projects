#! /usr/bin/env python


# consensusSequence.py
# P.Clote

# Computes consensus sequences for each k-point mutant sample set
# as well as for all sample structures over all mutants

import sys,os,math
from vienna import computeViennaSecStr

NUCL = ['A','C','G','U']

def main(filename):
  file    = open(filename)
  rna0    = file.readline().strip()
  secStr0 = computeViennaSecStr(rna0)
  D    = {}
  n = len(rna0); num = 0
  line = file.readline()
  numMut = 0; numSamples = 0
  while line:
    if line[0]=='>':
      if numMut == 0: #first time, get number of samples
        numSamples = int(line.split()[2])
      #Now add information for new number of mutations   
      words      = line.split()
      numMut     = int(words[-2]) #new number of mutations
      D[numMut]  = {}
      for i in range(n):
        D[numMut][i] = {}
        for ch in NUCL: D[numMut][i][ch] = 0.0
      num  = 0 #start over counter of number of mutations
      line = file.readline()
      continue
    num    += 1
    rna    = line.strip().upper()
    secStr = file.readline().strip()
    for i in range(n):
      D[numMut][i][rna[i]] += 1.0
    line = file.readline()
  file.close()
  #Normalize D values
  print "Consensus sequences for k-mutants, each k"
  for num in range(1,numMut+1):
    sys.stdout.write(">%d\n" % num)
    for i in range(n): 
      consensusNucl = ''; maxFreq = 0
      for ch in NUCL:
        D[num][i][ch] /= numSamples
        nuclFreq       = D[num][i][ch] #nucleotide frequency
        if nuclFreq>maxFreq:
          maxFreq = nuclFreq; consensusNucl = ch
      sys.stdout.write("%s" % consensusNucl)
    sys.stdout.write("\n")
  E = {} #E[i] is consensus nucleotide in position i over ALL mutations
  print "Consensus sequences over all mutants"
  print "> %d" % (numMut+1)
  for i in range(n):
    E[i] = {}
    for ch in NUCL: E[i][ch] = 0.0
    for num in range(1,numMut+1):
      for ch in NUCL:
	E[i][ch] += D[num][i][ch]
    maxFreq = 0; consensusNucl = ''
    for ch in NUCL:
      E[i][ch] = E[i][ch]/numMut  #normalize
      if E[i][ch]>maxFreq:
	maxFreq = E[i][ch]; consensusNucl = ch
    sys.stdout.write("%s" % consensusNucl)
  sys.stdout.write("\n")
	  
	 
if __name__ == '__main__':
  if len(sys.argv) < 2:
    text = """Usage: %s filename
    1) file contains first line of RNA sequence (unmutated), then
       all subsequent lines contain RNAmutants samples, where each sample
       is given by RNA seq and RNA sec str, each on separate lines."""
    text = text % sys.argv[0]
    print text
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)

