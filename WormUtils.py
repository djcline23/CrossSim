#!/usr/bin/env python
# encoding: utf-8
"""
WormUtils.py
"""

chromosome_phys_max = [15072000, 15279000, 13784000, 17494000, 20920000, 17719000]
cM_max = [47.0507, 53.92552, 53.84778, 47.44498, 51.69473, 52.22193]

#Checks for a valid chromsome number
def check_chromosome_number(chromNumber):
  if chromNumber < 1 or chromNumber > 6:
    raise ValueError, "Chromosome number must be within 1 to 6"

#Checks for a valid physical location
def check_physLoc(physLoc, chromNumber):
  if physLoc > chromosome_phys_max[chromNumber] or physLoc < 0:
    raise ValueError, "Physical location must be within the range of the chromosome"

#Gets the average percentages of the target chromsome and genome from a set of individuals  
def averagePercentages(diploidSet, targetChrom, targetName):
  length = len(diploidSet)
  totalSelected = 0
  totalGenome = 0
    
  for diploid in diploidSet:
    for chrSet in diploid.chromosome_set:
      totalSelected += chrSet[targetChrom].getPercentageOfParent(targetName)
        
    totalGenome += diploid.getPercentageOfGenome(targetName)
    
  return [totalSelected / (length * 2), totalGenome / length]