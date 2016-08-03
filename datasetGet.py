# -*- coding: utf-8 -*-
#!C:/Python27/
#!/usr/bin/python
"""
Created on Wed May 18 14:12:48 2016

@author: Matthew Carse
"""
#@ Class to read in UniRef identifiers and randomise using seed for reproducibility.
#@       !Uses UniRef cluster ids for representative sequence¡ e.g. UniRef90_Q587C9
#@ Randomised dataset identifiers sectioned into training (60%), validation (30%) and testing (10%) datasets.
#@ Sequence identifiers used to query UniProt and sequences downloaded and written to files with set names:
#@       training_seqs.txt
#@       validation_seqs.txt
#@       testing_seqs.txt
#@ Command line input requires the full name of a file in the current working directory (including file type e.g. .txt)
#@ Delimiter parameter is optional (default: \n)

import random, sys, getopt, os
from query import query

class randomiser:
    def main(self, argv):
        # defaults
        delimiter = "\n"
        mode = "split"
        pathname = None
        self.groupTotals = list()
        
        # evaluate parameters, giving error message and help message if incorrect option
        try:
            opts, args = getopt.getopt(argv,"hf:d:g:m:",["file=","delimiter=","groupings=","mode="])
        except getopt.GetoptError as err:
            print err
            self.usage()
            sys.exit(2)
            
        # iterate through options
        # if -h give help message showing correct usage
        # -f / --file mandatory - name of file containing identifiers
        # -d / --delimiter optional - non-default delimiter for parsing identifiers
        # -g / --groupings mandatory - total number of sequences according to each grouping
        # -m / --mode optional - traditional 60/30/10 split or cross-validation
        # additional check for if file is found
        for opt, arg in opts:
            if opt == '-h':
                self.usage()
                sys.exit()
            elif opt in ("-d", "--delimiter"):
                delimiter = arg
            elif opt in ("-f", "--file"):
                # check for existence of file
                if os.path.isfile(arg):
                    # get file path using current directory and option
                    pathname = ('{0}/{1}').format(os.getcwd(),arg)
                # if file does not exist, quit with error message
                else:
                    print "Error: file not found (please include in same working directory)"
                    sys.exit()
            elif opt in ("-g", "--groupings"):
                # separate list of totals
                g = arg.split(',')
                # use enumerate to get index and value
                for i, c in enumerate(g):
                    # use try/except to catch non-integer inputs
                    try:
                        self.groupTotals.extend([int(i)]*int(c))
                    except ValueError:
                        print "Error: groupings input must be numerical"
                        sys.exit()
            elif opt in ("-m", "--mode"):
                mode = arg.lower()         
        
        # quit with error message if no groupings are included in the options
        if len(self.groupTotals) == 0:
            print "Error: number of sequences per groupings required"
            self.usage()
            sys.exit()
            
        # once options and arguments have been parsed:
        # 1) if a filename has been entered then the file is opened and parsed
        # 2) the identifiers are added to a list and shuffled using a seed (6593)
        # 3) the method to split the dataset into training/validation/testing subsets is then called
        # 4) the identifiers from each set are then used to query for sequences
        # 5) the sequences are then written to files
        # the user is shown output of the first 5 identifiers prior to randomisation in order to check
        # whether the file parsing worked, and again after randomisation to see if the operation was performed
        # the user is also asked whether they want to continue as the process of querying for sequences may take a while
        if pathname is not None:
            ids = []
            print "Reading file " + pathname + " ...\n"
            f = open(pathname).read()
            for line in f.split(delimiter):
                ids.append(line)
            print "First 5 identifiers: " + ','.join(ids[0:5]) + "\n"
            if len(ids) == 1:
                print "Parsing appears not to have worked (this should only be 5 identifiers), did you use the correct delimiter? "
                var = raw_input("Enter y to continue / q to quit\n")
                if var.lower() == 'y':
                    pass
                else:
                    sys.exit()
            print "Randomising identifiers...\n"
            self.randomise(ids)
            print "First 5 identifiers again: " + ','.join(ids[0:5]) + "\n"
            print "Creating datasets..."
            # if mode = split then perform dataSplit
            if mode == "split":
                self.dataSplit(ids)
            print "Querying UniProt, parsing sequences..."
#            print "This may take a while for larger datasets"
#            var = raw_input("Enter y to continue / q to quit\n")
#            if var.lower() == 'y':
            if mode == "split":
                self.seqQuery("split", ids)
            if mode == "cv":
                self.seqQuery("cv", ids)
#            # the option says 'q to quit' but any input other than 'y' will cause the system to exit
#            # so that querying of UniProt does not accidentally begin
#            else:
#                sys.exit()                    
        else:
            print "Error: file name required"
            self.usage()
            sys.exit()
            
    def usage(self):
        print "Usage: python datasetGet.py -f <file name (including file type e.g. .txt)> -d <delimiter separating IDs (default \\n)> " \
        "-g <comma-separated number of proteins in each grouping (must be in the same order as sequence identifiers e.g. 386,410,393)> " \
        "-m <'split' or 'cv' for 60/30/10 training/validation/test set approach or cross-validation (default 'split')>"
        
       
    # randomise cluster ids using seed for reproducibility
    # also randomise groupings concurrently
    def randomise(self, ids):
        random.seed(6593)
        random.shuffle(ids)
        random.shuffle(self.groupTotals)
   
    # take 60%, 30% and 10% of ids for each set respectively
    # similarly section the groupings
    def dataSplit(self, ids):
        self.training = ids[0:int((len(ids)/10.0)*6)]
        self.validation = ids[int((len(ids)/10.0)*6):int((len(ids)/10.0)*9)]
        self.testing = ids[int((len(ids)/10.0)*9):len(ids)]
        self.training_groupings = self.groupTotals[0:int((len(self.groupTotals)/10.0)*6)]
        self.validation_groupings = self.groupTotals[int((len(self.groupTotals)/10.0)*6):int((len(self.groupTotals)/10.0)*9)]
        self.testing_groupings = self.groupTotals[int((len(self.groupTotals)/10.0)*9):len(self.groupTotals)]
           
    # query and store seqs and groupings
    # if m = split:
    # write sequences to files:
    #       training_seqs.txt
    #       validation_seqs.txt
    #       testing_seqs.txt 
    # if m = cv:
    # write sequences to files:
    #       all_seqs.txt
    def seqQuery(self, m, ids):
        if m == "split":
            query(self.training, 'training', self.training_groupings)
            query(self.validation, 'validation', self.validation_groupings)
            query(self.testing, 'testing', self.testing_groupings)
        if m == "cv":
            query(ids, 'all', self.groupTotals)
        
        
if __name__ == '__main__':
    r = randomiser()
    r.main(sys.argv[1:])
    
    
