#!/usr/bin/env python
# encoding: utf-8
"""
WormIndividual.py

Zifan Xiang
Copyright (c) 2014 Northwestern University. All rights reserved.
"""

from WormUtils import *
import os
  
if __name__ == '__main__':
  chromNumber = int(input("Enter the target chromosome: "))
  check_chromosome_number(chromNumber)
  chromNumber = chromNumber - 1;
  physLoc = int(input("Enter the target base pair (KB): "))
  check_physLoc(physLoc, chromNumber)
  crossStart = int(input("Enter the starting number of crosses: "))
  crossEnd = int(input("Enter the ending number of crosses: "))
  crossStep = int(input("Enter the number of steps between the start and end crosses: "))
  indStart = int(input("Enter the starting number of individuals (per cross): "))
  indEnd = int(input("Enter the ending number of individuals (per cross): "))
  indStep = int(input("Enter the number of steps between the start and end individuals: "))
  numIter = int(input("Enter the number of iterations of these crosses: "))
  
  for i in range(numIter):
    os.system("BasicCross.py %d %d %d %d %d %d %d %d" % (physLoc, chromNumber, crossStart, crossEnd, crossStep, indStart, indEnd, indStep))