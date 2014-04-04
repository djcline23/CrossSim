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
import getopt
import itertools
import operator
import os.path
import sys

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

#num represents if the desired intervals are the lower of higher intervals, num = 0/1
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

def backCrossSimulation():
  physLoc = int(sys.argv[1])
  chromNumber = int(sys.argv[2])
  numCrosses = range(int(sys.argv[3]), int(sys.argv[4]) + int(sys.argv[5]), int(sys.argv[5]))
  numIndividuals = range(int(sys.argv[6]), int(sys.argv[7]) + int(sys.argv[8]), int(sys.argv[8]))
  bucketSize = int(sys.argv[9])
  
  if os.path.isfile('general_statistics.csv'):
    g = open('general_statistics.csv', 'a')
  else:
    g = open('general_statistics.csv', 'wb')
    g.write('Number of Back Crosses,Number of Individuals,Selected Chromosome,Selected Base Pair,Percent Selected Chromosome,Percent Genome,Left Physical Loc,Right Physical Loc,Bucket Size, Left Distal, Left Proximal, Left Unique, Right Proximal, Right Distal, Right Unique \n')
   
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
                    
      hold = averagePercentages(Aparent, chromNumber, targetNameDip)
      averageTarget = hold[0]
      averageGenome = hold[1]
      g.write('%d,%d,%d,%d,%f,%f' % (crossNumber, indNumber, chromNumber + 1, physLoc, averageTarget, averageGenome))
      calculateAveragePhysicalIntervals(physIntervals, physLoc, chromNumber, g)
      putIntervalsIntoBuckets(physIntervals, bucketSize, g)

  g.close()

#Parameters: physLoc chromNumber numCrossStart numCrossEnd numCrossStep numIndStart numIndEnd numIndStep
if __name__ == '__main__':
  backCrossSimulation()
