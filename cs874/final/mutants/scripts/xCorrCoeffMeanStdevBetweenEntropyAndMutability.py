#! /usr/bin/env python


# xCorrCoeffMeanStdevBetweenEntropyAndMutability.py
# P.Clote

import sys,os,math
from misc import _RT_,basePairList
from vienna import computeViennaSecStr
from corrCoeff import corrCoeff
from stats import getSampleStats

DEBUG = 0

def xLogX(x):
  if x == 0:
    return 0
  else:
    return x*math.log(x)
  

def main(filename):
  file    = open(filename)
  rna0    = file.readline().strip()
  secStr0 = computeViennaSecStr(rna0)
  D    = {}; SS = {}; H = {} #H is entropy
  n = len(rna0); num = 0
  line = file.readline()
  numMut = 0; numSamples = 0
  while line:
    if line[0]=='>':
      if numMut == 0: #first time, get number of samples
        numSamples = int(line.split()[2])
      else: #update entropy H
        H[numMut] = {}
        for i in range(1,n+1): #indices 1<=i<=n
          sumEntropyForI = 0.0  #compute H[i]
          probII         = 1.0  #compute p(i,i) by 1.0 - p(i,j) all j
          for j in range(1,n+1):
            if i!=j: 
              probIJ = SS[numMut][i][j]/float(numSamples)
              sumEntropyForI += -xLogX(probIJ)
              probII         -= probIJ
          H[numMut][i] = sumEntropyForI + -xLogX(probII)
      #Now add information for new number of mutations   
      words      = line.split()
      numMut     = int(words[-2]) #new number of mutations
      D[numMut]  = {}
      SS[numMut] = {}
      for i in range(1,n+1):
        SS[numMut][i] = {}
        for j in range(1,n+1): SS[numMut][i][j] = 0.0
        D[numMut][i-1] = 0.0 
         #D[numMut][i] is number of mutations in position i
         #WARNING: indices in D are 0<=i<n, while those in SS are 1<=i<=n
      num  = 0 #start over counter of number of mutations
      line = file.readline()
      continue
    num    += 1
    if DEBUG: print numMut,num 
    rna    = line.strip().upper()
    secStr = file.readline().strip()
    bps    = basePairList(secStr)
    for i in range(n):
      if rna0[i]!=rna[i]: D[numMut][i] += 1.0
      baseI = i+1 #warning 1<=baseI<=n, but 0<=i<n for accessing RNA string
      for j in range(1,n+1):
        if baseI<j and (baseI,j) in bps:
          SS[numMut][baseI][j] += 1
        elif j<baseI and (j,baseI) in bps:
          SS[numMut][baseI][j] += 1
    line = file.readline()
  file.close()
  #Must complete computation of H for last value of numMut
  #As well, we normalize D values to be between [0,1]
  H[numMut] = {}
  for i in range(1,n+1): #indices 1<=i<=n
    sumEntropyForI = 0.0  #compute H[i]
    probII         = 1.0  #compute p(i,i) by 1.0 - p(i,j) all j
    for j in range(1,n+1):
      if i!=j:
        probIJ = SS[numMut][i][j]/float(num)
        sumEntropyForI += -xLogX(probIJ)
        probII         -= probIJ
    H[numMut][i] = sumEntropyForI + -xLogX(probII)
  #Normalize D values
  for num in range(1,numMut+1):
    for i in range(n): D[num][i] /= numSamples
  #Now compute correlation coefficient
  corrCoeffList = []
  print "k\tcorrCoeff\tmean1\t\tstdev1\t\tmean2\t\tstdev2"
  for num in range(1,numMut+1):
    L1 = []; L2 = []
    for i in range(n):
      L1.append(D[num][i])
      L2.append(H[num][i+1]) #WARNING: indices in H are from 1 to n
    mean1,stdev1,max1,min1 = getSampleStats(L1)
    mean2,stdev2,max2,min2 = getSampleStats(L2)
    print "%d\t%f\t%f\t%f\t%f\t%f" % (num,corrCoeff(L1,L2),mean1,stdev1,mean2,stdev2)
   

 
if __name__ == '__main__':
  if len(sys.argv) < 2:
    text = """Usage: %s filename 
    1) file contains first line of RNA sequence (unmutated), then
       all subsequent lines contain RNAmutants samples, where each sample
       is given by RNA seq and RNA sec str, each on separate lines.  """
    text = text % sys.argv[0]
    print text
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)

