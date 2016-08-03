# -*- coding: utf-8 -*-
#!C:/Python27/
#!/usr/bin/python
"""
Created on Wed May 18 2016

@author: Matthew Carse
"""
#@ Main class responsible for running the programs and checking whether
#@ all the required files exist (in the current directory):
#@ If run in 'split' mode:
#@      1) training_seqs.txt
#@      2) validation_seqs.txt
#@      3) testing_seqs.txt
#@      4) training_groupings.txt
#@      5) validation_groupings.txt
#@      6) testing_groupings.txt
#@ If run in cross-validation mode:
#@      1) all_seqs.txt
#@      2) all_groupings.txt
#@ The *_seqs.txt files contain the newline-delimited peptide sequences downloaded from UniProt
#@ The *_groupings.txt files contain the 'protein class' groupings for sequences, following randomisation,
#@ to validate the outcomes of machine learning classification
#@ If all files are found, the program creates an object which error checks the sequences

import sys, getopt, os
from vpGenList import vpGenList

class main:
    def main(self, argv):
        # defaults
        self.mode = None
        # parse options and arguments, giving error message if unknown options
        try:
            opts, args = getopt.getopt(argv,"hm:",["mode="])
        except getopt.GetoptError as err:
            print err
            print "Usage: python evolve.py -h <gives requirements to run>"
            sys.exit(2)
        # if valid option - help
        # help option gives requirements
        for opt, arg in opts:
            if opt == '-h':
                print "If mode is 'split':"
                print "This program requires the following files: \n\t1) training_seqs.txt\n\t2) validation_seqs.txt\n\t3) testing_seqs.txt\n" \
                "\t4) training_groupings.txt\n\t5) validation_groupings.txt\n\t6) testing_groupings.txt\n"
                print "If mode is 'cv':"
                print "This program requires the following files: \n\t1) all_seqs.txt\n\t2) all_groupings.txt\n"
                print "As produced by the randomiser.py script:"
                print "The *_seqs.txt files contain the newline-delimited peptide sequences downloaded from UniProt"
                print "The *_groupings.txt files contain the space-delimited 'protein group' groupings for sequences following " \
                "randomisation\n"
                print "Usage: python evolve.py -m <'split' or 'cv' for 60/30/10 training/validation/test set approach or cross-validation>"
                sys.exit()
            elif opt in ("-m", "--mode"):
                if arg.lower() in ('split', 'cv'):
                    self.mode = arg.lower()
                    self.fileCheck()
                else:
                    print "Mode must be one of 'split' or 'cv' (no quotes)"
                    print "Usage: python evolve.py -h <gives requirements to run>"
                    sys.exit()

        if self.mode == None:
            print "Mode must be one of 'split' or 'cv' (no quotes)"
            print "Usage: python evolve.py -h <gives requirements to run>"
            sys.exit()
        
            
    # check whether each of the required files exists
    def fileCheck(self):
        if self.mode == "split":
            if os.path.isfile("training_seqs.txt"):
                print "file: training_seqs.txt found"
            else:
                print "Error: file: training_seqs.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            if os.path.isfile("validation_seqs.txt"):
                print "file: validation_seqs.txt found"
            else:
                print "Error: file: validation_seqs.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            if os.path.isfile("testing_seqs.txt"):
                print "file: testing_seqs.txt found"
            else:
                print "Error: file: testing_seqs.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            if os.path.isfile("training_groupings.txt"):
                print "file: training_groupings.txt found"
            else:
                print "Error: file: training_groupings.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            if os.path.isfile("validation_groupings.txt"):
                print "file: validation_groupings.txt found"
            else:
                print "Error: file: validation_groupings.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            if os.path.isfile("testing_groupings.txt"):
                print "file: testing_groupings.txt found"
            else:
                print "Error: file: testing_groupings.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            # if all tests pass
            vpGenList("training_seqs.txt", "validation_seqs.txt", "testing_seqs.txt", "split")
        if self.mode == "cv":
            if os.path.isfile("all_seqs.txt"):
                print "file: all_seqs.txt found"
            else:
                print "Error: file: all_seqs.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            if os.path.isfile("all_groupings.txt"):
                print "file: all_groupings.txt found"
            else:
                print "Error: file: all_groupings.txt not found"
                print "Usage: python evolve.py -h <gives requirements to run>"
                sys.exit()
            # if all tests pass
            vpGenList("all_seqs.txt", None, None, "cv")
            
            
if __name__ == '__main__':
    m = main()
    m.main(sys.argv[1:])
