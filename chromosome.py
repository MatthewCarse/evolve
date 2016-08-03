# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 2016

@author: Matthew Carse
"""

#@ Chromosome object which holds the unique, random anchor-chain combinations
#@ and their resultant fitness metric (precision, recall, f1 score or accuracy)
#@ Chromosome objects are created and curated by a chromosomeList

class chromosome:      
    # store the ids (anchor-chain combinations) and fitness metric score
    def __init__(self, ids, fitness):
        self.ids = ids
        self.fitness = fitness   
    
    def __repr__(self):
        return repr((self.ids, self.fitness))