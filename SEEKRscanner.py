import glob
import collections
import itertools
import pandas as pd
import os
import numpy as np

import fasta_reader as far
import numpy as np
import pickle
from itertools import groupby
from tqdm import tqdm as tqdm


from SEEKR import SEEKR


class SEEKRscanner(SEEKR):

    '''
    Scans given sequences of interest (targets) for regions of high correlation
    to a set of known functional sequences (queries)

    This is done by fragmenting the targets into regions size specified by the
    user.

    Parameters
    ----------
    queries: path to valid fasta formatted file of query sequences (example:
    Xist Repeat A,B,C,D,E etc)
    targets: path to valid fasta formatted file of target sequences (example:
    sequence with suspected cis repressive functionality)
    k: value of k for the analysis, typically k = [4,6] default = 5
    seekr_reference: reference set of sequences to standardize k-mer counts, defaults
    to included dataset of gencode lncRNA files, but optionally can
    use your own reference set

    pearsonr_reference: fasta file that specifies the sequences to be used to build
    a distribution of pearson correlation values



    Attributes
    ----------
    windowLength: Size of the sliding window over a sequence
    slide: How base pairs to skip when sliding the window over the sequence.
    threshold: How many standard deviations above mean from reference
    Inherited:
    SEEKR.kmer_profile: kmer profile of sequences of interest
    SEEKR.reference: matrix of lenght normalized k-mer counts from reference file
    SEEKR.bases: sequence alphabet (for now only ATCG)
    SEEKR.keys: list of all possible kmers for given value of k
    SEEKR.seqnames: names of sequences supplied

    '''

    def __init__(self,motif_path,fasta_file,reference,k,windowLength,slide,threshold):
        super(SEEKRscanner,self).__init__(fasta_file,k,reference)

        self.motif_path = motif_path
        self.windowLength = windowLength
        self.slide = slide
        self.treshold = threshold
