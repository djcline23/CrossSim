#!/usr/bin/env python
# encoding: utf-8
"""
WormIndividual.py

Zifan Xiang
Copyright (c) 2014 Northwestern University. All rights reserved.
"""

from Chromosomes import *
from Individual import*
from numpy import *
from WormIndividual import *
from WormUtils import *
import itertools
import operator
import os.path
import sys
import thread
import threading

lock = threading.Lock();

def backCrossTillLimit(wormASet, wormB, physLoc, chromNumber, parent, limit):
  """Crosses each worm in set A with worm B until the limit number of offspring that keep the parent segment at the desired location has been met"""
  generation = []
  loc = Chromosome.getLoc(physLoc, chromNumber)
    
  while (len(generation) < limit):
    for wormA in wormASet:
      if (wormA.sex == "hermaphrodite"):
        curWorm = wormA.mate(wormB)
        
        if curWorm.chromosome_set[0][chromNumber].getParentAtLocation(loc) == parent or curWorm.chromosome_set[1][chromNumber].getParentAtLocation(loc) == parent:
          generation.append(curWorm)

    return generation

def backCrossTillLimitDiploid(diploidASet, diploidB, physLoc, chromNumber, parent, limit):
  """Crosses each worm in set A with worm B until the limit number of offspring that keep the parent segment at the desired location has been met"""
  generation = []
  loc = Chromosome.getLoc(physLoc, chromNumber)
    
  while (len(generation) < limit):
    for diploidA in diploidASet:
      curDiploid = diploidA.mate(diploidB)[0]
                
    if curDiploid.chromosome_set[0][chromNumber].getParentAtLocation(loc) == parent or curDiploid.chromosome_set[1][chromNumber].getParentAtLocation(loc) == parent:
      generation.append(curDiploid)
    
  return generation

def roundRobinCrossTillLimitDiploid(diploidSet, physLoc, chromNumber, parent, limit):
  generation = []
  loc = Chromosome.getLoc(physLoc, chromNumber)
    
  while (len(generation) < limit):
    for i in range(len(diploidSet) - 1):
      curDiploid = diploidSet[i + 1].mate(diploidSet[i])[0]
                
    if curDiploid.chromosome_set[0][chromNumber].getParentAtLocation(loc) == parent or curDiploid.chromosome_set[1][chromNumber].getParentAtLocation(loc) == parent:
      generation.append(curDiploid)
           
  return generation

def averagePercentages(diploidSet, targetChrom, targetName):
  length = len(diploidSet)
  totalSelected = 0
  totalGenome = 0
    
  for diploid in diploidSet:
    for chrSet in diploid.chromosome_set:
      totalSelected += chrSet[targetChrom].getPercentageOfParent(targetName)
        
    totalGenome += diploid.getPercentageOfGenome(targetName)
    
  return [totalSelected / (length * 2), totalGenome / length]

def writeGeneralStatistics(crossNumber, physLoc, diploidSet, targetChrom, targetName, bucketSize, g):
  indNumber = 1;
  genLoc = Chromosome.getLoc(physLoc, chromNumber)
  
  for diploid in diploidSet:
    totalSelected = 0
    curIntervals = []
  
    for chrSet in diploid.chromosome_set:
      percent = chrSet[targetChrom].getPercentageOfParent(targetName)
      print(percent)
      totalSelected += chrSet[targetChrom].getPercentageOfParent(targetName);
        
      if chrSet[chromNumber].getParentAtLocation(genLoc) == targetName:
        curIntervals.append(chrSet[chromNumber].physicalLocsOfInterval(genLoc, chromNumber))
            
    
    totalLower = 0;
    totalUpper = 0;
    for interval in curIntervals:
      totalLower += interval[0];
      totalUpper += interval[1]
      
    perGenome = diploid.getPercentageOfGenome(targetName)
    avgLower = totalLower / len(curIntervals)
    avgUpper = totalUpper / len(curIntervals)
    avgSelected = totalSelected / 2
    g.write('%d,%d,%d,%d,%f,%f, %d, %d\n' % (crossNumber, indNumber, targetChrom + 1, physLoc, avgSelected, perGenome, avgLower, avgUpper))
    #calculateAveragePhysicalIntervals(curIntervals, physLoc, targetChrom, g)
    #putIntervalsIntoBuckets(curIntervals, bucketSize, g)
    indNumber += 1
    

def writeGroupSegments(fileName, group):
  f = open(fileName, 'wb')
  t = 1 #Counter for the different individuals that have to printed within each cross 
  for diploid in group:
    l = 1 #Counter for the two sets of chromosomes
    f.write("Individual %d\n" % t)
    for chrSet in diploid.chromosome_set:
      f.write("Set %d Chr 1: %s\n" % (l, chrSet[0].segments))
      f.write("Set %d Chr 2: %s\n" % (l, chrSet[1].segments))
      f.write("Set %d Chr 3: %s\n" % (l, chrSet[2].segments))
      f.write("Set %d Chr 4: %s\n" % (l, chrSet[3].segments))
      f.write("Set %d Chr 5: %s\n" % (l, chrSet[4].segments))
      f.write("Set %d Chr 6: %s\n" % (l, chrSet[5].segments))
      l += 1
        
    t += 1
    
  f.close()

#num represents if the desired intervals are the lower or higher intervals, num = 0/1
def separatePhysicalInterval(selectedPhysInterval, num):
  toReturn = [0 for x in range(len(selectedPhysInterval))]
  
  for i in range(len(selectedPhysInterval)):
    toReturn[i] = selectedPhysInterval[i][num]
    
  return toReturn

def bucketPhysicalIntervals(selectedPhysInterval, low, high, bucketSize):
  size = int(high - low) / bucketSize
  buckets = [0 for i in range(size)]
  
  if (len(buckets) == 0):
    buckets = [0]
  
  for i in range(len(selectedPhysInterval)):
    bucketPos = int((selectedPhysInterval[i] - low) / bucketSize)
    
    if bucketPos == len(buckets):
      bucketPos -= 1

    buckets[bucketPos] += 1

  return buckets

def putIntervalsIntoBuckets(physInterval, bucketSize, filePath):
    
  filePath.write(',%d' % bucketSize)
  
  for i in range(2):
    sepPhysicalInterval = separatePhysicalInterval(physInterval, i)
    low = min(sepPhysicalInterval)
    high = max(sepPhysicalInterval)
    toWrite = bucketPhysicalIntervals(sepPhysicalInterval, low, high, bucketSize)
    numUnique = 0;
    
    for bucket in toWrite:
      if bucket != 0:
        numUnique += 1;
        
    filePath.write(',%d,%d,%d' % (low, high, numUnique))
        
  filePath.write('\n')
    
def calculateAveragePhysicalIntervals(physIntervals, physLoc, chromNumber, filePath):
  lowerIntervalSums = 0
  upperIntervalSums = 0
  lowerIntervalAverages = 0
  upperIntervalAverages = 0
    
  for x in range(len(physIntervals)):
    interval = physIntervals[x]
    lowerIntervalSums += interval[0]
    upperIntervalSums += interval[1]

    lowerIntervalAverages = lowerIntervalSums / len(physIntervals)
    upperIntervalAverages = upperIntervalSums / len(physIntervals)

  filePath.write(',%d,%d' % (lowerIntervalAverages, upperIntervalAverages))

def backCrossSimulation(physLoc, chromNumber, numCrosses, numIndividuals, bucketSize):
  if os.path.isfile('general_statistics_%d_%d.csv' % (physLoc, chromNumber + 1)):
    g = open('general_statistics_%d_%d.csv' % (physLoc, chromNumber + 1), 'a')
  else:
    g = open('general_statistics_%d_%d.csv' % (physLoc, chromNumber + 1), 'wb')
    g.write('Number of Back Crosses,Individual Number,Selected Chromosome,Selected Base Pair,Percent Selected Chromosome,Percent Genome,Left Physical Loc,Right Physical Loc,Bucket Size, Left Distal, Left Proximal, Left Unique, Right Proximal, Right Distal, Right Unique \n')

  #Runs through the number of crosses specified and makes the individuals
  for crossNumber in numCrosses:
    physIntervals = []
    for indNumber in numIndividuals:
      Aparent = [ Diploid(name = "A", newChr = 6) ]
      Bparent = Diploid(name = "B", newChr = 6)
      targetNameDip = Aparent[0].name
      genLoc = Chromosome.getLoc(physLoc, chromNumber)

      for k in range(crossNumber):
        Aparent = backCrossTillLimitDiploid(Aparent, Bparent, physLoc, chromNumber, targetNameDip, indNumber)
    
      # Format of the output files is as follows: Number of Crosses_ Number Of Individuals per Cross _ Target Chromosome _ Physical Location on the Target Chromosome
      fileName = "%d_%d_%d_%d_crossConfig.csv" % (crossNumber, indNumber, chromNumber + 1, physLoc)
      writeGroupSegments(fileName, Aparent)

      for diploid in Aparent:
        for chrSet in diploid.chromosome_set:
          if chrSet[chromNumber].getParentAtLocation(genLoc) == targetNameDip:
            physIntervals.append(chrSet[chromNumber].physicalLocsOfInterval(genLoc, chromNumber))
                    
      
      writeGeneralStatistics(crossNumber, physLoc, Aparent, chromNumber, targetNameDip, bucketSize, g)
      #hold = averagePercentages(Aparent, chromNumber, targetNameDip)
      #averageTarget = hold[0]
      #averageGenome = hold[1]
      #calculateAveragePhysicalIntervals(physIntervals, physLoc, chromNumber, g)
      #putIntervalsIntoBuckets(physIntervals, bucketSize, g)

  g.close()

class CrossThread (threading.Thread):
    def __init__(self, threadID, physLoc, chromNumber, numCrosses, numIndividuals, bucketSize):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.physLoc = physLoc
        self.chromNumber = chromNumber
        self.numCrosses = numCrosses;
        self.numIndividuals = numIndividuals
        self.bucketSize = bucketSize
    
    def run(self):
        print(self.threadID)
        backCrossSimulation(self.physLoc, self.chromNumber, self.numCrosses, self.numIndividuals, self.bucketSize)

#Parameters: physLoc chromNumber numCrossStart numCrossEnd numCrossStep numIndStart numIndEnd numIndStep bucketSize
if __name__ == '__main__':
  physLoc = int(sys.argv[1])
  chromNumber = int(sys.argv[2]) - 1;
  numCrosses = range(int(sys.argv[3]), int(sys.argv[4]) + int(sys.argv[5]), int(sys.argv[5]))
  #TODO(zifanxiang): Will keep the numIndividuals as a range but the inputs will just have one value in the range
  #TODO(zifanxiang): Potentially change to just one number
  numIndividuals = range(int(sys.argv[6]), int(sys.argv[7]) + int(sys.argv[8]), int(sys.argv[8]))
  bucketSize = int(sys.argv[9])
  numIter = int(sys.argv[10])
  #numThreads = int(sys.argv[11])
  
  for i in range(numIter):
    backCrossSimulation(physLoc, chromNumber, numCrosses, numIndividuals, bucketSize)
  
  #threads = []
  
  #for i in range(numThreads):
  #  threads.append(CrossThread('%d' % i, physLoc, chromNumber, numCrosses, numIndividuals, bucketSize))
  
  #for i in range(numIter / numThreads):
  #  for i in range(numThreads):
  #    threads[i].start()
    
  #for t in threads:
  #  t.join()
    
