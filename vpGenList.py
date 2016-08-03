# -*- coding: utf-8 -*-
"""
Created on Wed May 18 2016

@author: Matthew Carse
"""
#@ The vpGenList class produces a list of all possible anchor-chain residue combinations
#@ and accordingly creates 400 vpGen objects which are used to produce every possible
#@ feature set for each of the sequences within the training set, read in from the
#@ file: training_seqs.txt using the seqRead class, and likewise for the validation dataset
#@ and test dataset.
#@ This approach allows the feature sets to be stored in memory and not recalculated
#@ each time. When a series of random anchor-chain combinations are produced,
#@ the feature sets for each sequence can be looked up and combined into a 
#@ 6n-feature/n-feature-set chromosome (by the chromosomeList class).
#@ This chromosome is standardised and the scaling factors retained (using the lda class).


from vpGen import vpGen
from machineLearning import machineLearning
from seqRead import seqRead
from chromosomeList import chromosomeList
from chromosome import chromosome
import time


class vpGenList:
    def __init__(self, trainDataset, validDataset, testDataset, mode):
        self.train = trainDataset   
        self.valid = validDataset
        self.test = testDataset
        self.mode = mode
        time1 = None
        
        # first step - produce list of anchor-chain combinations and resultant vpGens
        # residueCombos and allvpGens are class variables so operation should only need to be performed once
        if (len(residueCombos) != 400) & (len(allvpGens) != 400):
            self.residueGeneration()
            time1 = time.time()
            
        # second step - read sequences, error check and create feature sets
        # all*FSets needs only be created once as it is a class variable
        # if mode = cross-validation, all seq feature sets are stored in allTrainFSets
        sR = seqRead(allvpGens, residues)
        global allTrainFSets
        if len(allTrainFSets) == 0:
            allTrainFSets = sR.seqRead(self.train)   
            
        if self.mode == "split":
            # send allTrainFSets for preprocessing (standardisation), returning scaled feature sets dictionary list
            l = machineLearning(self.mode)
            global allTrainFSetsScaled
            if len(allTrainFSetsScaled) == 0:
                allTrainFSetsScaled = l.preprocess()
    
            # read in validation dataset sequences using same method
            global allValidFSets
            if len(allValidFSets) == 0:
                allValidFSets = sR.seqRead(self.valid)
            # send validation set for scaling
            global allValidFSetsScaled
            if len(allValidFSetsScaled) == 0:
                allValidFSetsScaled = l.scale('valid')
                print ('Time for creation of vpGens and scaling: {}s').format(time.time() - time1)
                
            # read in testing dataset sequences
            global allTestFSets
            if len(allTestFSets) == 0:
                allTestFSets = sR.seqRead(self.test)
            # send test set for scaling
            global allTestFSetsScaled
            if len(allTestFSetsScaled) == 0:
                allTestFSetsScaled = l.scale('test')
                print ('Time for creation of test set and scaling: {}s').format(time.time() - time1)    
           
        # once all sequences have been read from files and preprocessing has occurred:   
        # create chromosomeList object to produce chromosomes and start genetic algorithm
        c = chromosomeList(mode)
        
        ids = c.geneticAlgorithm()

        # create chromosome using highest-scoring ids and report fitness when used to classify test dataset
        if self.mode == "split":
            print "Prediction fitness for highest-scoring ids using linear discriminant analysis and random forest classification:"
            print self.evaluate(ids)
            print "Producing confusion matrices..."
            print time.time()-time1      
            
        if self.mode == 'cv':
            print "Producing Receiver Operating Characteristic curve..."
            self.evaluate(ids)
            print time.time()-time1
        
       
        
        
    # method to generate all possible residue combinations and according vpGens    
    def residueGeneration(self):
        # create all possible anchor-chain combinations
        for residue in residues[:len(residues)-1]:
            for residue2 in residues[:len(residues)-1]:
                residueCombos.append(residue+residue2)
        print 'residueCombos created!'
                
        # create all possible vpGen objects using anchor-chain combinations, AA ... YY
        for combo in residueCombos:
            allvpGens[combo[0]+combo[1]] = vpGen(combo[0],combo[1])
        print 'allvpGens created!'
        
    # function to evaluate the best model produced by the genetic algorithm
    # if split mode is active, the test set is analysed, producing a fitness value
    # for both linear discriminant analysis and random forest classification and
    # a confusion matrix for each classifier
    # if using cross-validation, the program produces a roc curve and a confusion matrix
    def evaluate(self, ids):
        if self.mode == 'split':
            trainChromosomeList = list()
            validChromosomeList = list()
            testChromosomeList = list()
            # recreate the models, using the ids, from the scaled feature sets
            for seq in allTrainFSetsScaled:
                fset = []
                for i in range(len(ids)):                
                    fset += seq[ids[i]]
                trainChromosomeList.append(fset)
                
            for seq in allValidFSetsScaled:
                fset = []
                for i in range(len(ids)):
                    fset += seq[ids[i]]
                validChromosomeList.append(fset)
                
            for seq in allTestFSetsScaled:
                fset = []
                for i in range(len(ids)):
                    fset += seq[ids[i]]
                testChromosomeList.append(fset)
                
            l = machineLearning(self.mode)
            ### CHANGE: testFitness/testFitness2 to validFitness/validFitness2 in return statement       ###
            ###         to instead evaluate validation dataset using highest-scoring ids                 ###
#            l.cm = True # change variable to true to toggle plotting of confusion matrix
#            validFitness = l.classify(trainChromosomeList, validChromosomeList) # lda model for validation set
#            validFitness2 = l.classifyRFC(trainChromosomeList, validChromosomeList) # random forest model for validation set
            testFitness = l.testEvaluateLDA(trainChromosomeList, testChromosomeList) # lda model for test set
            testFitness2 = l.testEvaluateRFC(trainChromosomeList, testChromosomeList) # random forest model for test set
            return chromosome(ids, testFitness), chromosome(ids, testFitness2) # return both lda and random forest chromosomes for highest-scoring ids
#            return chromosome(ids, validFitness), chromosome(ids, validFitness2) # lda and rfc chromosomes for validation set
            
        if self.mode == 'cv':
            trainChromosomeList = list()
            # recreate models from unscaled feature sets
            for seq in allTrainFSets:
                fset = []
                for i in range(len(ids)):                
                    fset += seq[ids[i]]
                trainChromosomeList.append(fset) 
       
            l = machineLearning(self.mode)
            l.rocCurve(trainChromosomeList)
            l.cm = True # toggle variable to true to plot confusion matrix
            l.classifyLDA(trainChromosomeList, list())
       
# class variables - accessible by every instance of the class when run in parallel
    
# amino acid residue codes used to create anchor-chain combinations
# can be used to check for non-peptide sequences, anchors or chains      
# includes X as this is an ambiguity code for unknown residues
residues = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']
residueCombos = list()
# variable to hold all 400 possible vpGens in order to produce every feature set for each sequence
# dictionary used to look up vpGen objects using their identifier
allvpGens = dict()
# empty list to contain (number of sequences) * dictionary of feature sets produced using allvpGens
allTrainFSets = list()
allTrainFSetsScaled = list()
allValidFSets = list()
allValidFSetsScaled= list()
allTestFSets = list()
allTestFSetsScaled = list()

