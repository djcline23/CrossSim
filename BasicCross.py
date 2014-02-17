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

"""The following four arrays contain the cross configuration parameters for each back cross"""
numCrosses = range(2, 15, 1)
numIndividuals = range(4, 2, 20)
physLocs = []
chromNumbers = range(0, 5, 1)

def setUpPhysLocs():
    i = 0
    
    for cM in Chromosome.cM_max:
        for j in range(10, 1, -1):
            physLocs[i][10 - j] = cM / j
    
    
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

if __name__ == '__main__':
  
  for crossNumber in numCrosses:
    for indNumber in numIndividuals:
        for physLoc in physLocs:
            for chromNumber in chromNumbers:
                Aparent = [ Diploid(name = "A", newChr = 6) ]
                Bparent = Diploid(name = "B", newChr = 6)
                targetNameDip = Aparent[0].name

                for i in range(crossNumber):
                    Aparent = backCrossTillLimitDiploid(Aparent, Bparent, physLoc, chromNumber, targetNameDip, numIndividuals)
    
                fileName = "%d_%d_%d_%d_crossConfig", (crossNumber, indNumber, physLoc, chromNumber)
                f = open('')
                for diploid in Aparent:
                    i = 1
                    for chrSet in diploid.chromosome_set:
                        print "Set %d Chr 1: %s" % (i, chrSet[0].segments)
                        print "Set %d Chr 2: %s" % (i, chrSet[1].segments)
                        print "Set %d Chr 3: %s" % (i, chrSet[2].segments)
                        print "Set %d Chr 4: %s" % (i, chrSet[3].segments)
                        print "Set %d Chr 5: %s" % (i, chrSet[4].segments)
                        print "Set %d Chr 6: %s" % (i, chrSet[5].segments)
                        i += 1
  