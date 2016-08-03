# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 2016

@author: Matthew Carse
"""

#@ Class with functions to mutate and recombine sequences.
#@ Mutation randomly replaces one of the anchor-chain combinations
#@ in a chromosome with one from the list of 400 possibilities.
#@ Recombination involves cleaving two sequences using a common
#@ index, taking from the beginning of the chromosome to the index
#@ from the first sequence, and from the index to the end of the
#@ chromosome from the second sequence, and forming a new, third sequence
#@ from these two portions. The remaining portions of each sequence
#@ are used to create a fourth sequence.

import vpGenList as vpGenList
import random

class mutateRecombine:      
    def __init__(self):
        pass 
        
    # frequency of mutation - (slim) possibility of the same sequence being mutated twice
    # p = 0.96 (1/25 chance of mutation)
    # p = 0.98 (1/50 chance of mutation)
    # p = 0.99 (1/100 chance of mutation)
    # p = 0.995 (1/200 chance of mutation)   
    # function which takes a sequence and randomly replaces an anchor-chain combination
    def mutate(self, seq):
        # randomly pick an identifier
        # randomly sample a new identifier from the list
        # replace the old with the new
        self.p = 0.98
        if random.random() >= self.p:
            i = random.randint(0, len(seq)-1)
            n = ''.join(random.sample(vpGenList.residueCombos,1))
            seq[i] = n
            return [seq, True] # return True if mutation occurs
        else:
            return [seq, False] # return False if mutation did not occur
        
    # function which takes two sequences and creates a further two through recombination
    # a cleavage point is chosen at random between the first and length-2 indices
    def recombine(self, cSome1, cSome2):
        # pick an index to cleave the first sequence
        # using indices 1 to length-2, otherwise would be duplicating the original sequences
        seq1, seq2 = cSome1.ids, cSome2.ids
        i = random.randint(1, len(seq1)-2)
        seq3 = seq1[:i] + seq2[i:]
        seq4 = seq2[:i] + seq1[i:]
        return seq3, seq4
