# :hatching_chick: evolve :hatched_chick:
##Genetic algorithm for heuristic optimisation of machine-learning classification of protein sequences

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

```
Files:
  datasetGet.py
  query.py
```

The genetic algorithm program `evolve.py` is likewise called through the command line, taking a singular argument:
- -m/--mode (split or cv) for the program to be run with training set used to fit the machine-learning models, the validation set for validation of these models, and the testing set as held-out data for evaluation of the highest-scoring model following running of the genetic algorithm. The cv option is for cross-validation in which the whole dataset is used for training and validation, with different sequences taken for each fold. The program must be run from the directory containing the sequence/groupings files, allowing a relative path to the python script.
```
python evolve.py -m cv
```

```
Files:
  chromosome.py
  chromosomeList.py
  evolve.py
  machineLearning.py
  mutateRecombine.py
  seqRead.py
  vpGen.py
  vpGenList.py
```

### Required packages
- [Scikit-learn] (http://scikit-learn.org/stable/)
- [NumPy] (http://www.numpy.org/)
- [Plotly] (https://plot.ly/python/getting-started/)

### Parameters
- `query.py` - change `self.email` to the user's email address for querying UniProt servers
- `chromosomeList.py` - `self.origLength` = number of chromosomes [models] in population
- `chromosomeList.py` - `self.retLength` = number of 'top' chromosomes (by fitness) kept for recombination
- `chromosomeList.py` - `self.vpGenNum` = number of vpGen objects [6-feature sets] per model 
- `chromosomeList.py` - stopping criterion - stops program if the top 10 models reach an average fitness of 1.0
``` python
while (sum(c.fitness for c in chromosomes[0:10]))/10 < 1.0:
``` 
- `chromosomeList.py` - stopping criterion - stops program if the top 1% of models all have the same ids (anchor-chain combinations) as diversity has been lost
``` python
sets = [c.ids for c in chromosomes[0:(len(chromosomes)/100)]]

if all(sets[0] == sets[i+1] for i in range((len(chromosomes)/100)-1)):
``` 
- `chromosomeList.py` - stopping criterion - stops program if the mean fitness of the top 1% of models has not changed since the prior genetic algorithm generation
```python 
meanFitness.append((sum(c.fitness for c in chromosomes[0:(len(chromosomes)/100)]))/(len(chromosomes)/100))

if self.gaPasses >= 2:

if meanFitness[self.gaPasses-1] == meanFitness[self.gaPasses-2]:
``` 
- `machineLearning.py` - if run using cross-validation, uses 3 folds 
```python 
predicted = cross_validation.cross_val_predict(clf, tCList, trainGroupings, cv=3)
```
- `machineLearning.py` - if using split mode, the highest-scoring model is evaluated on the test set using random forest classification with 5 estimators
```python 
rfc = RandomForestClassifier(n_estimators=5, random_state=1, max_features=None)
```
- `mutateRecombine.py` - `self.p` = mutation probability (1-p) e.g. p=0.98, mutation rate = 1/50

### Outputs
- `datasetGet.py` produces training_seqs.txt, validation_seqs.txt, testing_seqs.txt, training_groupings.txt, validation_groupings.txt, testing_groupings.txt or all_seqs.txt, all_groupings.txt depending on preparation mode
- `evolve.py` produces an HTML interactive _Plotly_ graph of each chromosome model of each genetic algorithm generation by fitness, with the ids and fitness as annotation. The program also produces a confusion matrix and ROC curve if using cross-validation or two confusion matrices if using held-out data
- If the HTML graphs do not automatically load - append the url (e.g. https://github.com/MatthewCarse/evolve/blob/master/Sample%20outputs/CV/6fold_cv_roc.html) to http://htmlpreview.github.io/? to view
https://htmlpreview.github.io/?https://github.com/MatthewCarse/evolve/blob/master/Sample%20outputs/CV/6fold_cv_roc.html

