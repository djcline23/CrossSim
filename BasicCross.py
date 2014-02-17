#!/usr/bin/env python
# encoding: utf-8
"""
WormIndividual.py

Zifan Xiang
Copyright (c) 2014 Northwestern University. All rights reserved.
"""

from WormIndividual import *
from Individual import*
from Chromosomes import *
import operator
import itertools
from numpy import *

def backCrossTillLimit(wormASet, wormB, physLoc, chromNumber, parent, limit):
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

if __name__ == '__main__':
  #A = [ Worm(name="A", sex="hermaphrodite") ]
  #B = Worm(name="B", sex="male")
  Aparent = [ Diploid(name = "A", newChr = 6) ]
  Bparent = Diploid(name = "B", newChr = 6)
  #targetName = A[0].name
  targetNameDip = Aparent[0].name
  print(targetNameDip)
  loc = Chromosome.getLoc(6500, 2)
  print(loc);
  #for i in range(50):
    #A = backCrossTillLimit(A, B, 8000, 0, targetName, 20)

  for i in range(75):
    Aparent = backCrossTillLimitDiploid(Aparent, Bparent, 6500, 2, targetNameDip, 20)
    
  for i in range(10):
    Aparent = roundRobinCrossTillLimitDiploid(Aparent, 6500, 2, targetNameDip, 20)
    
  for diploid in Aparent:
    i = 1
    for chrSet in diploid.chromosome_set:
        print "Set %d Chr 2: %s" % (i, chrSet[1].segments)
        print "Set %d Chr 3: %s" % (i, chrSet[2].segments)
        i += 1
  
  #for worm in A:
  #  i = 1
  #  for chrSet in worm.chromosome_set:
  #          print "Set %d Chr 1: %s" % (i, chrSet[0].segments)
  #          i += 1
        
