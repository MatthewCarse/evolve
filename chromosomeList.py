# -*- coding: utf-8 -*-
"""
Created on Sun Jul 03 2016

@author: Matthew Carse
"""
#@ ChromosomeList class handles chromosome objects for use with the genetic algorithm.
#@ The first step involves producing a number of chromosome objects with random anchor-chain 
#@ residue combinations which are replaced until unique.
#@ The unique, random anchor-chain combinations (ids) are used to create a machine-learning
#@ model (linear discriminant analysis) for prediction of the classification of the validation
#@ dataset, and the resultant fitness metric is recorded and stored within the chromosome object.
#@ The chromosomes are then ranked according to their prediction fitness and only a proportion
#@ of the top chromosomes are kept while the rest are discarded. This retained population is used
#@ to generate new chromosomes until the original population size is reached, with the likelihood
#@ of a chromosome being picked for recombination dependent upon its fitness. 
#@ The new and retained population of chromosomes are sorted and assessed to determine whether
#@ any of the stopping criteria are met: fitness threshold, lack of diversity, lack of progress.
#@ If any criterion is met, an html graph is produced using plotly and the program quits. 


from machineLearning import machineLearning
from mutateRecombine import mutateRecombine
from chromosome import chromosome
from plotly.graph_objs import Scatter
import vpGenList as vpGenList
import random, sys, time, numpy as np
import plotly as py


class chromosomeList:
    def __init__(self, mode):
        self.mode = mode
    
    
    def geneticAlgorithm(self):
        # first step - generate self.origLength number of chromosomes with self.vpGenNum random, unique anchor-chain combinations
        # each self.vpGenNum-id set is used to create an lda model and predict the classifications of the validation dataset
        # using the same anchor-chain combinations
        # a fitness metric (precision, recall, f-statistic, accuracy) is calculated from the prediction and
        # attached to the chromosome object
        # all chromosome objects are stored in a chromosomes list variable
    
        ### VARIABLES TO CHANGE ###
        # self.origLength = number of chromosomes in population
        # self.retLength = number of 'top' chromosomes (by fitness) kept for recombination
        # self.vpGenNum = number of vpGen objects per model
        self.origLength, self.retLength, self.vpGenNum, self.time1 = 1000, 125, 10, time.time()
        print ('Original population: {}, retained population: {}, vpGens: {}').format(self.origLength, self.retLength, self.vpGenNum)

        # create self.origLength number of chromosomes by calling the idsGenerator function
        chromosomes = [self.chromosomeGenerator(self.idsGenerator(), True) for i in range(self.origLength)]
        # once all chromosomes have been created - sort chromosomes according to fitness metric
        chromosomes = sorted(chromosomes, key = lambda chromosome: chromosome.fitness, reverse = True)
        # print the top 10 chromosomes by fitness for the user
        print chromosomes[0:10]
        print ''
        
        
        # store fitness metric for each chromosome for plotting
        fit = [c.fitness for c in chromosomes]

        # plot first iteration using plotly subplots, with ids and fitness as annotation
        self.fig = py.tools.make_subplots(shared_xaxes=True, shared_yaxes=True, print_grid=False)
        trace1 = Scatter(y=fit,x=range(1,len(fit)+1), text=[str(c) for c in chromosomes], name = ('GA pass: {}').format(0))
        self.fig.append_trace(trace1, 1, 1)  

       
       
        meanFitness = list() # empty list to hold the meanFitness for the top 1% of each genetic algorithm pass
        self.gaPasses = 0 # counter for number of runs of genetic algorithm
       
       
       
        # create weightings for recombination such that the higher the fitness of a chromosome, the more likely it
        # is to be picked - the first is thrice as likely (*2.998) as the last with 1000 chromosomes retained
        # calculate the last value and subtract this divided by retLength from all values in order to fit all weightings in
        # otherwise results in a probability of greater than 1
        last = ((1.0/self.retLength)/2 + (1.0/self.retLength ** 2))
        weightings = [((1.0/self.retLength)/2 + (float(self.retLength - i)/self.retLength ** 2) - float(last)/self.retLength) for i in range(self.retLength)] 
        # add the difference from 1 to the first weighting (a negligible amount, equal to (last - (last/self.retLength))/(self.retLength/2))
        # or ~ 2 * (last/self.retLength)
        weightings[0] += (1 - sum(weightings))



        print "Time to start of genetic algorithm: ", time.time()-self.time1

        # second step - use genetic algorithm by which the top self.retLength of chromosomes, according to the fitness metric,
        # are selected and used to repopulate to the original population size using recombination
        # the genetic algorithm runs until one of the stopping criteria are met:
            # fitness metric threshold (e.g. 0.8)
            # lack of chromosome diversity (all of the top 1% of chromosomes are identical)
            # lack of progression (the top 1% of chromosomes stay the same for two runs of the genetic algorithm)
              
        ### STOPPING CRITERION TO BE CHANGED ###
        # run genetic algorithm until mean 1.0 f1 score for the top 10 chromosomes (both can be changed for a different
        # threshold or a different number of chromosomes)
        while (sum(c.fitness for c in chromosomes[0:10]))/10 < 1.0:
            # increment counter for number of times genetic algorithm has been called          
            self.gaPasses += 1
            gaStart = time.time()
            m = mutateRecombine() # mutate object
            
            
            # take top self.retLength of original chromosomes by fitness
            chromosomes = chromosomes[:self.retLength]
            
            
            # mutate remaining original population            
            for i, c in enumerate(chromosomes):
                # store result (sequence and True or False)
                res = m.mutate(c.ids)
                # if mutation occurs
                if res[1] == True:
                    # recalculate fitness and replace chromosome
                    # pass False so mutation only occurs once
                    chromosomes[i] = self.chromosomeGenerator(res[0], False)

            print "Time to start of recombination: ", time.time()-gaStart
            
            # repopulate list to original length using recombination
            while len(chromosomes) < self.origLength:
                # choose two chromosomes from retained population, with higher-fitness chromosomes more likely
                # to be picked, according to the weightings (replace=False so the same chromosome is not picked)
                choice1, choice2 = np.random.choice(chromosomes[:self.retLength], size=2, replace=False, p=weightings)
                ids1, ids2 = m.recombine(choice1, choice2) # generate new ids using recombination method in mutate class
                chromosomes.append(self.chromosomeGenerator(ids1, True))
                chromosomes.append(self.chromosomeGenerator(ids2, True))

            print "Time to start of sorting: ", time.time()-gaStart
            
            # sort according to fitness and print top 10 chromosomes for user
            chromosomes = sorted(chromosomes, key = lambda chromosome: chromosome.fitness, reverse = True)
            print chromosomes[0:10]
            print ''

            print "Time to fitness extraction and plotting: ", time.time()-gaStart
            
            # extract each fitness metric value and plot on the graph
            fit = [c.fitness for c in chromosomes]
            
            # plot chromosomes by fitness with ids/fitness annotation and line name dynamically produced using self.gaPasses
            trace1 = Scatter(y=fit,x=range(1,len(fit)+1), text=[str(c) for c in chromosomes], name = ('GA pass: {}').format(self.gaPasses))
            self.fig.append_trace(trace1, 1, 1)

            print "Time to convergence check: ", time.time()-gaStart
            
            
            ### STOPPING CRITERION TO BE CHANGED ###
            # note time of convergence - if top 1% of chromosomes are identical then diversity has been lost
            # can change the percentage of chromosomes to check for diversity
            # produce a list of the top 1% of ids
            sets = [c.ids for c in chromosomes[0:(len(chromosomes)/100)]]
            # iterate through all sets of ids, comparing them to the first
            # if all are True, diversity has been lost and the time to convergence is noted and the program exits
            if all(sets[0] == sets[i+1] for i in range((len(chromosomes)/100)-1)):
                convTime = ('Time to convergence: {}s\n').format(time.time()-self.time1)
                print convTime
                self.savePlot(m.p) # call function to save graph
                                
                return chromosomes[0].ids # return highest-scoring ids

            print "Time to plateau check: ", time.time()-gaStart
            
            # if ids are unique for recombination, change to use set of fitness in top 1%
            
            ### STOPPING CRITERION TO BE CHANGED ###
            # append the mean fitness over the top 1% of chromosomes to a list in order to check performance
            # over previous genetic algorithm passes - if current performance is the same as the performance
            # from 1 pass ago, quit (the number of genetic algorithm passes over which there is no change can
            # be altered - recommended one pass for large populations due to long running time)
            meanFitness.append((sum(c.fitness for c in chromosomes[0:(len(chromosomes)/100)]))/(len(chromosomes)/100))  
            # only checks for mean fitness performance if there have been at least 2 algorithm passes 
            # (also requires changing in line with the number below otherwise will fail)
            if self.gaPasses >= 2:
                # change the second number to alter the number of algorithm iterations
                # note: self.gaPasses-1 is the index of the current algorithm pass
                if meanFitness[self.gaPasses-1] == meanFitness[self.gaPasses-2]:
                    platTime = ('Time to plateau: {}s\n').format(time.time()-self.time1)               
                    print platTime
                    self.savePlot(m.p) # call function to save graph
                    
                    return chromosomes[0].ids # return highest-scoring ids
                                        
            print "Genetic algorithm pass time: ", time.time()-gaStart
            print ''
            
        # if threshold f-score is reached
        self.savePlot(m.p) # call function to save graph
        return chromosomes[0].ids # return highest-scoring ids
        
    
    # function to give plotly plot a dynamically produced title and filename and save as html
    def savePlot(self, mutationRate):
        title = ('Prediction fitness - {}s, {} ga passes').format(time.time() - self.time1, self.gaPasses)
        xaxis = dict(title='Chromosome')
        yaxis = dict(title='Fitness')
        self.fig['layout'].update(title=title,xaxis=xaxis,yaxis=yaxis)
        saveName = ('pred_f1_{}_{}_sorted_{}vpGen_{}_{}').format(self.origLength, self.retLength, self.vpGenNum, int(round(1/(1-mutationRate))), self.mode)
        py.offline.plot(self.fig, filename=saveName)
        
    
    # function to check whether the id of each scaled feature set is unique
    # if the id is not unique, fetch another id combination and update ids
    def unique(self, ids):
        # default set
        seen = set()
        seen_add = seen.add
        # read through list of anchor-chain identifiers
        for i, id in enumerate(ids):
            # if id has already been used:
            # fetch another id
            # update list of identifiers
            if id in seen:
                ids[i] = self.idGenerator()
            # if id hasn't previously been seen, add to set
            else:
                seen_add(id)
        return ids
        
        
    # return two random combinations from the available peptide residues
    def idGenerator(self):
        # use all but one residue code as last is for ambiguity (X)
        anchor = ''.join(random.sample(vpGenList.residues[:len(vpGenList.residues)-1],1))
        chain = ''.join(random.sample(vpGenList.residues[:len(vpGenList.residues)-1],1))
        return anchor+chain
        
        
    # function to produce self.vpGenNum (defined above) random, unique anchor-chain combinations     
    def idsGenerator(self):
        # call idGenerator function for the defined number of times to produce id combinations
        ids = [self.idGenerator() for i in range(self.vpGenNum)]
        # check that all anchor-chain identifiers are unique
        while len(set(ids)) < len(ids):
            ids = self.unique(ids)
        return ids
        
        
    # function which creates a chromosome - an object which stores self.vpGenNum random, unique anchor-chain
    # combinations used to create an lda classification model, and the resultant fitness metric value
    # mutable is a True or False variable indicating whether mutation should occur
    # new chromosomes are mutated by the function whereas existing chromosomes are not
    # and are instead mutated with the same chance but in a way that does not require 
    # remaking all 'retained' chromosomes in the genetic algorithm
    def chromosomeGenerator(self, ids, mutable):
        if mutable:
            # chance to mutate chromosome before classification
            m = mutateRecombine()
            # store result (sequence + True or False depending on mutation)
            res = m.mutate(ids)
            # if mutation occurs
            if res[1] == True:
                ids = res[0]
        # check that all anchor-chain identifiers are unique
        while len(set(ids)) < len(ids):
            ids = self.unique(ids)
        # sort ids into alphabetical order
        ids.sort()    
          
        # if cross-validation mode - training set = all sequences
        if self.mode == "cv":
            trainChromosomeList = list()
            for seq in vpGenList.allTrainFSets:
                fset = []
                for i in range(len(ids)):                
                    fset += seq[ids[i]]
                trainChromosomeList.append(fset)  
          
#        print ids
            
#        # use the self.vpGenNum unique, random combinations to produce a chromosome for each sequence
#        # storing the feature set chromosome for each sequence in a list, the scaled feature sets are
#        # fetched from the class variable in vpGenList
        if self.mode == "split":
            trainChromosomeList = list()
            validChromosomeList = list()
            for seq in vpGenList.allTrainFSetsScaled:
                fset = []
                for i in range(len(ids)):                
                    fset += seq[ids[i]]
                trainChromosomeList.append(fset)
                
            for seq in vpGenList.allValidFSetsScaled:
                fset = []
                for i in range(len(ids)):
                    fset += seq[ids[i]]
                validChromosomeList.append(fset)
            
        # classification  
        # if length of trainChromosomeList is 0, no valid sequences were parsed
        # if sequences are present:
        # if split mode - pass to lda object to create model and classify validation chromosome
        # if cv mode - create and predict training chromosome using cross-validation
        if len(trainChromosomeList) == 0:
            print "Error: no valid sequences"
            sys.exit()
        else:
            l = machineLearning(self.mode)
            if self.mode == "cv":
                fitness = l.classifyLDA(trainChromosomeList, list()) # lda model for whole dataset (cross-validation)
            else:
                fitness = l.classifyLDA(trainChromosomeList, validChromosomeList) # lda model for validation set
            return chromosome(ids, fitness)
