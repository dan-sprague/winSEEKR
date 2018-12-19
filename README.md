# SEEKR
A collection of tools for RNA/DNA sequence analysis

SEEKR = **SE**quence **E**valuation by **K**-mer **R**epresentation

See Jessime Kirk's & Calabrese Labs paper here: https://www.nature.com/articles/s41588-018-0207-8

## Requirements

Python3.6 (Anaconda distribution recommended for other libraries)

**cd to the directory you cloned this repo to** and install requirements using the command:

pip install -r requirements.txt

## Included Test Data

Includes
1. test_xistandrsxrepeats.fa -> represents your experimental sequences of interest
2. test_reference_gencodelncrna.fa -> represents whatever reference set of sequences you wish to use
3. pwms/ -> directory containing motifs and metadata files 

# Scripts
For standard analysis techniques using SEEKR<br/>
## seekranalysis.py -fa -ref -k -savecorr
Generates a correlation matrix based on the k-mer content of given sequences<br/><br/>

-fa : fasta file containing experimental sequences of interest<br/>
-ref : fasta file containing sequences to be used as a reference for k-mer z-score calculation<br/>
-k : value of k to use, typically between 4 and 6 inclusive<br/>
-savecorr : filename for correlation matrix to be saved in working directory, if not specified matrix will only be printed to terminal<br/>

**Usage Example**
```
python seekranalysis.py -fa ./Test_Data/test_xistandrsxrepeats.fa -ref ./Test_Data/test_reference_gencodelncrna.fa  -k 5
```

## motifscanner.py -fa -ref -k -pwm -pwmmeta --savedf
Generates a dataframe containing Position Weighted Matrix (PWM) weighted z-scores for supplied sequences.<br/><br/>

-fa : fasta file containing experimental sequences of interest<br/>
-ref : fasta file containing sequences to be used as a reference for k-mer z-score calculation<br/>
-k : value of k to use, typically between 4 and 6 inclusive<br/>
-pwms : path to **folder** containing position weight matrices (PWMs)<br/>
-pwmmeta : txt file containing the metadata for each motif<br/>
--savedf : save name for dataframe containing the results, to be saved in working directory, defaults to dataframe.csv <br/>

**Usage Example**

```
python motifscanner.py -fa ./Test_Data/test_xistandrsxrepeats.fa -ref ./Test_Data/test_reference_gencodelncrna.fa -k 4 -pwms ./Test_Data/pwms/Homo_sapiens_2018_10_23_1-29_pm/pwms_all_motifs/ -pwmmeta ./Test_Data/pwms/Homo_sapiens_2018_10_23_1-29_pm/RBP_Information_all_motifs.txt
```

# API

## SEEKR

Three parameters are required to use SEEKR:
1. Fasta file of sequences you are studying 
2. Value of k, typically between 4 and 6
3. Fasta file of sequences that comprise a reference set, for example the GENCODE lncRNA annotation fasta file

The 'reference' fasta is used to get the mean and standard deviation for each k-mer in the SEEKR algorithm


Create a SEEKR instance by calling:
```
from SEEKR import SEEKR
SEEKRobject = SEEKR(fasta_file='/path/to/fasta.fa',k=5,reference='/path/to/ref.fa')
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
Creates a matrix of k-mer counts that serves as a reference for the mean and standard deviation of each k-mer's frequency. Example reference sets might be GENCODE lncRNA annotations, mRNA transcripts, etc...

If a path to reference fasta is passed, function is called at instantiation and stored as an attribute

```
reference = SEEKRobject.reference
```

### correlate_kmerprofiles(*savename*=*str*)

Takes the k-mer profile of a SEEKR object and calculates all pairwise Pearon's correlations. If a savename is passed, a csv file is saved in the current directory.

### save_kmer_profile(*savename*=*str*)
Save a pickle of the current kmer_profile attribute in the current directory

### save_ref(*savename* = *str*)
Save a pickle of the reference matrix in the current directory

## proteinSEEKR

For use in identifying likely proteins that bind a given DNA/RNA sequence. Motifs MUST be in the format used in the CISBP-Rna database (http://cisbp-rna.ccbr.utoronto.ca) 


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

```pSEEKRobject.map_proteins('/path/to/RBP_info.txt',pSEEKRobject.get_scoredataframe(weighted_zscores))```

