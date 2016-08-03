# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 2016

@author: Matthew Carse
"""
#@ Class to read peptide sequences from a file.
#@ Each sequence in the dataset is checked for numbers, length and composition.
#@ If the sequence contains >25% abundance of a base (A, C, T or G) the user is 
#@ asked whether the sequence should be included.
#@ The sequences which pass inspection each have a dictionary with 400 keys 
#@ (each residue-residue combination) and a feature set for each key,
#@ generated using the featureGen function of the vpGen class.
#@ The dictionary for each sequence is added to a list and returned to vpGenList 

import os, sys

class seqRead:
    def __init__(self, vpGens, residues):
        self.allvpGens = vpGens
        self.residues = residues
    
    # use the sequences from the dataset along with all possible anchor/chain combinations to create 400 feature sets
    # the feature sets are stored in the form of a dictionary, allowing look-up using the id combination
    def seqRead(self, dataset):
        allTrainFSets = allValidFSets = allTestFSets = list()
        
        # read file from current directory
        f = open(("{0}/{1}").format(os.getcwd(),dataset)).read()
        for line in f.split('\n'):
            # error check sequence, including:
                # empty line
                # length (10 residues or greater)
                # containing numbers
                # non-peptide sequence
                # potential nucleotide sequence
                # upper case
            # discounting empty lines
            if len(line) > 0:
                # get result of sequence composition check for later use
                res = self.seqComp(line.upper())
                if len(line) < 10:
                    print ('Sequence: {}\nIs not of sufficient length (10). Disregarding.').format(line)
                    raw_input("Press Enter to continue...")
                elif self.hasNumbers(line):
                    print ('Sequence: {}\nContains numbers. Disregarding.').format(line)
                    raw_input("Press Enter to continue...")
                elif not self.residueCheck(line.upper()):
                    print ('Sequence: {}\nContains non-amino acid code(s). Disregarding.').format(line)
                    raw_input("Press Enter to continue...")
#                elif not res[0]:
#                    print ('Sequence: {0}\nAppears to contain an overabundance (>25%) of residue {1}.').format(line.upper(),res[1])
#                    var = 'd'
#                    while str(var.lower().strip()) not in ['y', 'n', 'q']:
#                        var = raw_input('Is this a nucleotide sequence? (y to include sequence / n to disregard sequence / q to quit)\n')
#                    else:
#                        if str(var.lower().strip()) == 'y':
#                            d = {i: v.featureGen(line) for i, v in self.allvpGens.items()}
#                            if dataset == "training_seqs.txt":
#                                allTrainFSets.append(d)
#                            elif dataset == "validation_seqs.txt":
#                                allValidFSets.append(d)
#                            elif dataset == "testing_seqs.txt":
#                                allTestFSets.append(d)
#                        elif str(var.lower().strip()) == 'n':
#                            continue
#                        else:
#                            sys.exit()
                else:
                    # variable d is a dictionary with the vpGen id as a key
                    # each of the 400 vpGen anchor-chain combinations is iterated through
                    # the line is passed to each and the feature set added to the dictionary
                    # the dictionary is then appended to the list for all sequences
                    d = {i: v.featureGen(line) for i, v in self.allvpGens.items()}
                    if dataset in ("training_seqs.txt", "all_seqs.txt"):
                        allTrainFSets.append(d)
                    elif dataset == "validation_seqs.txt":
                        allValidFSets.append(d)
                    elif dataset == "testing_seqs.txt":
                        allTestFSets.append(d)
                        
        if dataset in ("training_seqs.txt", "all_seqs.txt"):
           print 'allTrainFSets created!'
           return allTrainFSets
        if dataset == "validation_seqs.txt":
           print 'allValidFSets created!'
           return allValidFSets
        if dataset == "testing_seqs.txt":
           print 'allTestFSets created!'
           return allTestFSets
      
      
    # before using sequence, check that it does not contain any numbers
    def hasNumbers(self, seq):
        return any(char.isdigit() for char in seq)  


    # function to check sequence composition
    # if the sequence contains an overabundance (>25%) of a base (A, C, T or G)
    # the test is 'failed' (though the user can override this)
    def seqComp(self, seq):
        if (seq.count('A')/float(len(seq))*100) >= 25:
            return [False, 'A']
        elif (seq.count('C')/float(len(seq))*100) >= 25:
            return [False, 'C']
        elif (seq.count('T')/float(len(seq))*100) >= 25:
            return [False, 'T']
        elif (seq.count('G')/float(len(seq))*100) >= 25:
            return [False, 'G']
        else:
            return [True]
            
            
    # function to check whether sequence contains only amino acid residues
    # UniProt appears to include X if residue is unknown
    def residueCheck(self, seq):
        for i, c in enumerate(seq):
            if c not in self.residues:
                return False
        return True       
           