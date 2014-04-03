#!/usr/bin/env python
# encoding: utf-8
"""
WormUtils.py
"""

chromosome_phys_max = [15072, 15279, 13784, 17494, 20920, 17719]
cM_max = [47.0507, 53.92552, 53.84778, 47.44498, 51.69473, 52.22193]

def check_chromosome_number(chromNumber):
  if chromNumber < 1 or chromNumber > 6:
    raise ValueError, "Chromosome number must be within 1 to 6"

def check_physLoc(physLoc, chromNumber):
  if physLoc > chromosome_phys_max[chromNumber] or physLoc < 0:
    raise ValueError, "Physical location must be within the range of the chromosome" 