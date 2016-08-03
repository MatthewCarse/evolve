# -*- coding: utf-8 -*-
"""
Created on Mon May 23 12:44:34 2016

@author: Matthew Carse
"""
#@ Class used by datasetGet.py to query UniProt with UniRef representative cluster id,
#@ by way of xml, extracting the peptide sequence for the cluster
#@ The object takes as parameters a dataset and the name of this dataset
#@ The program uses the dataset name to write the appropriate file containing the sequences
#@ The seed-randomised groupings are written to similarly named files, space-separated

import urllib2, os, xml
from xml.dom import minidom

class query:
    # takes dataset (identifiers), name (dataset name) and groupings (protein groups) as input
    def __init__(self, dataset, name, groupings):
        self.dataset = dataset
        self.name = name
        self.groupings = groupings
        # could require email as option
        self.email = "2213142C@student.gla.ac.uk"
        if self.name == 'training':
            fo = open((('{}/training_seqs.txt').format(os.getcwd())), 'w')
            self.writer(fo)
            fo.close()
            print (('{}/training_seqs.txt').format(os.getcwd())) + " created"
            
            fo = open((('{}/training_groupings.txt').format(os.getcwd())), 'w')
            for g in self.groupings:
                fo.write(('{} ').format(g))
            fo.close()
            print (('{}/training_groupings.txt').format(os.getcwd())) + " created"
        elif self.name == 'validation':
            fo = open((('{}/validation_seqs.txt').format(os.getcwd())), 'w')
            self.writer(fo)
            fo.close()
            print (('{}/validation_seqs.txt').format(os.getcwd())) + " created"
            
            fo = open((('{}/validation_groupings.txt').format(os.getcwd())), 'w')
            for g in self.groupings:
                fo.write(('{} ').format(g))
            fo.close()
            print (('{}/validation_groupings.txt').format(os.getcwd())) + " created"
        elif self.name == 'testing':
            fo = open((('{}/testing_seqs.txt').format(os.getcwd())), 'w')
            self.writer(fo)
            fo.close()
            print (('{}/testing_seqs.txt').format(os.getcwd())) + " created"
            
            fo = open((('{}/testing_groupings.txt').format(os.getcwd())), 'w')
            for g in self.groupings:
                fo.write(('{} ').format(g))
            fo.close()
            print (('{}/testing_groupings.txt').format(os.getcwd())) + " created"
        elif self.name == 'all':
            fo = open((('{}/all_seqs.txt').format(os.getcwd())), 'w')
            self.writer(fo)
            fo.close()
            print (('{}/all_seqs.txt').format(os.getcwd())) + " created"
            
            fo = open((('{}/all_groupings.txt').format(os.getcwd())), 'w')
            for g in self.groupings:
                fo.write(('{} ').format(g))
            fo.close()
            print (('{}/all_groupings.txt').format(os.getcwd())) + " created"
            
    # method to call query for each id and write result to file
    def writer(self, file):
        delCounter = 0 # count number of deletions
        for i, id in enumerate(self.dataset):
            res, seq = self.query(id)
            # if parsing xmldoc was a success, write sequence to file
            if res:
                file.write(('{}\n').format(seq))
            # if parsing failed, remove grouping
            else:
                # as the index is not being removed from the dataset
                # past the first removal, the indices no longer align
                # so remove the index minus the number of deletions
                del self.groupings[i - delCounter]
                delCounter += 1
    
    # function to parse xml document and extract peptide sequence  
    def parser(self, xmldoc):
        # get sequence
        seq = xmldoc.getElementsByTagName('sequence')
        # split and join to remove inherent xml sequence whitespace
        return str(''.join(str(seq[0].firstChild.nodeValue).split()))
        
    # query uniprot with uniref representative cluster id (repeat for training/validation/testing) to obtain xml document    
    def query(self, id):
        success = True # variable to track success of xml parsing
        url = 'http://www.uniprot.org/uniref/?query={}&format=xml'.format(id)
        request = urllib2.Request(url)
        contact = self.email 
        request.add_header('User-Agent', 'Python %s' % contact)
        response = urllib2.urlopen(request)
        # attempt to parse xmldoc, catching error if does not exist
        try:
            xmldoc = minidom.parse(response)
            return success, self.parser(xmldoc)
        except xml.parsers.expat.ExpatError:
            print "Unable to find xml document for id: {} - removing from dataset".format(id) # discard id
            success = False # change to false if failure
            return success, None
        
