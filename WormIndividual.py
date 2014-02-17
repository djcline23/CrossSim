#!/usr/bin/env python
# encoding: utf-8
"""
WormIndividual.py

Created by Joshua Shapiro on 2008-08-11.
Copyright (c) 2008 Princeton University. All rights reserved.
"""
from Individual import Diploid, newChromosomes
from Chromosomes import *
import operator
import itertools
from numpy import *

class Worm(Diploid):
    """ Makes diploid karyotypes"""
    def __init__(self,  name = None, chromosome_set = None, newChr = None, cM = 200, chrNames = None, sex= None):
        self.name = name
        self.sex = sex
        if (chromosome_set == None and sex == None) or (chromosome_set != None and sex != None):
            raise ValueError, "Must specify only one of either sex or a list of the chromosomes themselves."
        elif chromosome_set != None:
            self.chromosome_set = chromosome_set
            length = len(self.chromosome_set[1])
            if (length == 5):
                self.sex = "male"
            elif (length == 6):
                self.sex = "hermaphrodite"
            else:  raise ValueError, "sex must be either male or hermaphrodite"
        elif (newChr == None):
            if (self.sex == "hermaphrodite"):
                self.chromosome_set = (newChromosomes(parent = self.name, n = 6, cM = cM, chrNames = chrNames, interference="complete"),
                                       newChromosomes(parent = self.name, n = 6, cM = cM, chrNames = chrNames, interference="complete"))
            elif (self.sex == "male"):
                self.chromosome_set = (newChromosomes(parent = self.name, n = 6, cM = cM, chrNames = chrNames, interference="complete"),
                                      newChromosomes(parent = self.name, n = 5, cM = cM, chrNames = chrNames, interference="complete"))
          
           
        for chr_list in self.chromosome_set:
            chr_list.sort(key= operator.attrgetter("name"))
       
       
    def make_gamete(self):
        gamete = []
        recombines = random.binomial(1, 0.5, len(self.chromosome_set[1]))
        if (self.sex == "hermaphrodite"):
            for maternal, paternal, crossover in itertools.izip (self.chromosome_set[0], self.chromosome_set[1], recombines):
                if (crossover):
                    gamete.append(maternal.recombine(paternal, interference = "complete")[0])
                else:  gamete.append((maternal, paternal)[random.binomial(1, 0.5)])
            return gamete        
        elif (self.sex == "male"):
            x = 0
            while (x < 5):
                if (recombines[x] == 1):
                   gamete.append(self.chromosome_set[0][x].recombine(self.chromosome_set[1][x], interference = "complete")[0])
                else:  
                   gamete.append((self.chromosome_set[0][x], self.chromosome_set[1][x])[random.binomial(1, 0.5)])
                x += 1
            gamete_sex = random.binomial(1, 0.5)
            if (gamete_sex == 1):
                gamete.append(self.chromosome_set[0][5])
                return gamete
            else:
                return gamete
    
    # TODO(zifanxiang): Easy refactor of this code with mate, this is a just a special case of mate
    def mate_self(self, nOffspring = 1):
        offspring = []
        if (self.sex == "male"):
            raise ValueError, "Only hermaphrodites can self."
        else:
            if nOffspring == 1:
                egg = self.make_gamete()
                sperm = self.make_gamete()
                offspring = Worm(chromosome_set = (egg, sperm))
            else:
                for i in range(nOffspring):
                  egg = self.make_gamete()
                  sperm = self.make_gamete()
                  offspring.append(Worm(chromosome_set = (egg, sperm)))
            return offspring
    
    def mate(self, dad, nOffspring = 1):
        offspring = []
        if (self.sex == "male"):
            raise ValueError, "First parent must be a hermaphrodite."
        elif (dad.sex == "hermaphrodite"):
            raise ValueError, "Second parent must be a male."
        else:
            if nOffspring == 1:
                egg = self.make_gamete()
                sperm = dad.make_gamete()
                offspring = Worm(chromosome_set = (egg, sperm))
            else:
                for i in range(nOffspring):
                  egg = self.make_gamete()
                  sperm = dad.make_gamete()
                  offspring.append(Worm(chromosome_set = (egg, sperm)))
            return offspring
    
    # TODO(zifanxiang): Figure out inheritance in python
    def get_all_genos(self, interval= 0.1):
      return self.getAllGenos(interval = interval, cM = False)
    
if __name__ == '__main__':
  N2 = Worm(name="N2", sex="hermaphrodite")
  CB = Worm(name="CB", sex="male")
  F1_herm = N2.mate(CB)
  F1_male = N2.mate(CB)
  while (F1_herm.sex == "male"):
    F1_herm = N2.mate(CB)
  while (F1_male.sex == "hermaphrodite"):
    F1_male = N2.mate(CB)
  F2 = F1_herm.mate(F1_male, nOffspring=1)
  
  #print F2.get_all_genos(interval = .1)
  print "Offspring Chromosomes:"
  for set in F1_herm.chromosome_set: 
    for chr in set:
      print "Chr %s: %s" % (chr.name, chr.segments)
