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

chromosome_phys_max = [15072, 15279, 13784, 17494, 20920, 17719];
cM_max = [47.0507, 53.92552, 53.84778, 47.44498, 51.69473, 52.22193]

"""The following four arrays contain the cross configuration parameters for each back cross"""
numCrosses = range(2, 7, 1)
numIndividuals = range(4, 6, 1)
physLocs = [];
chromNumbers = range(0, 6, 1)

def setUpPhysLocs():
    i = 0
    
    for kB in chromosome_phys_max:
        for j in range(1, 11, 1):
            physLocs.append((kB / 10) * j)
        
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
    print(loc)
    
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

if __name__ == '__main__':
  setUpPhysLocs();
  g = open('General_statisics', 'w')
  g.write('{:<30}'.format('Number of Back Crosses'))
  g.write('{:<35}'.format('Number of Individuals'))
  g.write('{:<40}'.format('Selected Chromosome'))
  g.write('{:<26}'.format('Selected Base Pair'))
  g.write('{:<40}'.format('Percent Selected Chromosome'))
  g.write('{:<26}'.format('Percent Genome'))
  g.write('\n')

  for crossNumber in numCrosses:
    for indNumber in numIndividuals:
        i = 0
        for chromNumber in chromNumbers:
            for j in range(0, 10, 1):
                Aparent = [ Diploid(name = "A", newChr = 6) ]
                Bparent = Diploid(name = "B", newChr = 6)
                targetNameDip = Aparent[0].name

                for k in range(crossNumber):
                    Aparent = backCrossTillLimitDiploid(Aparent, Bparent, physLocs[(i * 10) + j], chromNumber, targetNameDip, indNumber)
    
                # Format of the output files is as follows: Number of Crosses_ Number Of Individuals per Cross _ Target Chromosome _ Physical Location on the Target Chromosome
                fileName = "%d_%d_%d_%d_crossConfig" % (crossNumber, indNumber, chromNumber + 1, physLocs[(i * 10) + j])
                f = open(fileName, 'w')
                t = 1 #Counter for the different individuals that have to printed within each cross 
                for diploid in Aparent:
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
                    
                    t = t + 1
                
                f.close()
                hold = averagePercentages(Aparent, chromNumber, targetNameDip)
                averageTarget = hold[0]
                averageGenome = hold[1]
                g.write('{:<38}'.format('%d' % crossNumber))
                g.write('{:<35}'.format('%d' % indNumber))
                g.write('{:<40}'.format('%d' % (chromNumber + 1)))
                g.write('{:<26}'.format('%d' % physLocs[(i * 10) + j]))
                g.write('{:<35}'.format('%f' % averageTarget))
                g.write('{:<26}'.format('%f' % averageGenome))
                g.write('\n')

            i = i + 1
  
  g.close()