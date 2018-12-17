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

### kmer_profile
Creates a dictionary entry for each sequence-of-interest, where the key is the name of the sequence as specified in the .fasta file and the values are a vector of k-mer counts (length 4^k)

This function is called when a SEEKR object is created and exists as an attribute of the SEEKR object

### generate_ref *reference*=/path/to/reference.fa
Creates a matrix of k-mer counts that serves as a reference set to standardize k-mer counts in sequences-of-interest. Example reference sets might be GENCODE lncRNA annotations, mRNA transcripts, etc...


### save_kmer_profile *savename*=*str*
Save a pickle of the current kmer_profile attribute in the current directory

### save_ref *savename* = *str*
Save a pickle of the reference matrix in the current directory
