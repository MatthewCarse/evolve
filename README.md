# evolve
Genetic algorithm for heuristic optimisation of machine-learning classification of protein sequences

The .py scripts above contain two separate programs:

- `datasetGet.py` - dataset acquisition using UniRef protein cluster accession numbers
- `evolve.py` - genetic algorithm program for heuristic classification optimisation

The `datasetGet.py` program produces files for use by `evolve.py` but is not run as part of a pipeline. The machine learning for the genetic algorithm can use the standard approach with held-out data - training/validation/testing, or cross-validation. 
As the program is still in development, currently parameters must be changed within the code itself, as detailed in comments and in the parameters section.

The dataset acquisition program is called through the command line, taking:
- -f/--file argument containing the protein accession numbers
- -d/--delimiter argument (default \n) for the delimiter between accession numbers 
- -g/--groupings argument for the number of sequences belong to each protein class
- -m/--mode argument (split or cv) to determine which files are produced - split=training/validation/testing, cv=all
```
python datasetGet.py -f "Combined ids.txt" -d , -g 386,410,393 -m split
```

The genetic algorithm program `evolve.py` is likewise called through the command line, taking a singular argument:
- -m/--mode (split or cv) for the program to be run with training set used to fit the machine-learning models, the validation set for validation of these models, and the testing set as held-out data for evaluation of the highest-scoring model following running of the genetic algorithm. The cv option is for cross-validation in which the whole dataset is used for training and validation, with different sequences taken for each fold. 
```
python evolve.py -m cv
```


### Parameters
- `query.py` - change self.email to the user's email address for querying UniProt servers
- `chromosomeList.py` - self.origLength = number of chromosomes [models] in population
- `chromosomeList.py` - self.retLength = number of 'top' chromosomes [models] (by fitness) kept for recombination
- `chromosomeList.py` - self.vpGenNum = number of vpGen objects [6-feature sets] per model 
- `chromosomeList.py` - stopping criterion - 
- `chromosomeList.py` - 
- `chromosomeList.py` - 
