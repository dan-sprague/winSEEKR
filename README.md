# SEEKR
A collection of tools for RNA/DNA sequence analysis

## Requirements

Python3.6 (Anaconda distribution recommended for other libraries)

All requirements can be installed with pip using:

pip install -r requirements.txt


### Inputs: 

#### SEEKR 

1. Fasta file of sequences of interest (any valid fasta format)
2. Fasta file of reference set of sequences (any valid fasta format)
3. Value of k

#### proteinSEEKR

1. All requirements of SEEKR
2. path to motif file in plain text

motifs MUST be in the format used in the CISBP-Rna database (http://cisbp-rna.ccbr.utoronto.ca) 


# API

## SEEKR

Create a SEEKR object by calling:
```
from SEEKR import SEEKR
SEEKRobject = SEEKR('/path/to/fasta.fa',k=5,reference='/path/to/ref.fa')
```

If you have a saved reference matrix, pass the path to the pickle file:

```
from SEEKR import SEEKR
SEEKRobject = SEEKR('/path/to/fasta.fa',k=5,reference='/path/to/saved/pickle.p')
```


### kmer_profile()
Creates a dictionary entry for each sequence-of-interest, where the key is the name of the sequence as specified in the .fasta file and the values are a vector of k-mer counts (length 4^k)

This function is called when a SEEKR object is created and exists as an attribute of the SEEKR object and can be accessed by:

```
profile = SEEKRobject.kmer_profile
```

### generate_ref(*reference*='/path/to/reference.fa')
Creates a matrix of k-mer counts that serves as a reference set to standardize k-mer counts in sequences-of-interest. Example reference sets might be GENCODE lncRNA annotations, mRNA transcripts, etc...

If a path to reference fasta is passed, function is called at instantiation and stored as an attribute

```
reference = SEEKRobject.reference
```

### correlate_kmerprofiles(*savename*=*str*)

Takes the k-mer profile of a SEEKR object and calculates the Pearson correlation between them. If a savename is passed, a csv file is saved in the current directory.

### save_kmer_profile(*savename*=*str*)
Save a pickle of the current kmer_profile attribute in the current directory

### save_ref(*savename* = *str*)
Save a pickle of the reference matrix in the current directory

## proteinSEEKR

For use in identifying likely proteins that bind a given DNA/RNA sequence

To create a proteinSEEKR object:

```
pSEEKRobject = proteinSEEKR('path/to/motifs.txt', 'path/to/fasta/file.fasta',k,reference = '/path/to/ref.fa')
```

### read_motifs()
Returns a dictionary of PWMs with filenames as keys and pandas dataframes as values

This method is called at instantiation and is stored in the attribute **pwms**

### gen_kmersdict()
Returns a dictionary with kmers ('AAAA','AAAG','AAAC','AAAT',...,'TTTT') as keys, and 0 for each value. For use in counting kmers and called in other functions

### calc_prob(p,j,pwm,kmer)
Traverses a matrix M to multiply Mi+j,c values together to produce a crude probability for finding a kmer within a motif.

p = probability and passed as None, j = window within the motif currently under consideration, pwm = matrix, kmer = current k-mer under consideration

This function is called from kmer_probabilities 

### kmer_probabilities()
Returns a dictionary of kmers as keys and a crude probability of finding that k-mer within a given PWM as values

If the value of k is larger than the size of a given motif, for example a 4 base pair motif and k = 5, each k-mer is broken down into constitutive 4-mers, and then the probability of finding each 4mer within the motif are summed to represent the likelihood of the 5mer occuring. 

### save_kmerprobs(kmerprobdict,savename)
Saves a given dictionary of k-mer probabilities (as calculated in kmer_probabilities) as a pickle file with a given savename in the current directory
 
### get_zscore_df()
Converts k-mer profile dictionary into a pandas dataframe indexed by k-mer names

### kmer_weights(probabilities)
Returns a dictionary with the PWM name as the key and probability weighted z-scores as the values. These values represent a combined measure of the enrichment of a k-mer and likelihood of binding to that k-mer

``` pSEEKRobject.kmer_weights(pSEEKRobject.kmer_probabilities())```
### get_scoredataframe(scoredict,seqnames)
Convert dictionary of k-mer weighted zscores to a pandas dataframe, name columns by sequence names

### map_proteins(path,scoredf):
Map protein names to their associated motif PWMs. Required a path to the correct .txt file providing the information and a dataframe of k-mer weighted z-scores

```pSEEKRobject.map_proteins('/path/to/RBP_info.txt',pSEEKRobject.get_scoredataframe(weighted_zscores)```

