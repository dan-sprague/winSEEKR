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
SEEKRobject = SEEKR(/path/to/fasta.fa,k=5,reference='/path/to/ref.fa')
```

If you have a saved reference matrix, simply do not specify the reference variable and the pickle file will be loaded from the current directory


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
