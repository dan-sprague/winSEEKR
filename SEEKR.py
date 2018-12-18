import glob
import collections
import itertools
import pandas as pd
import os
import numpy as np

import fasta_reader as far
import numpy as np
import kmerprofilelib
import pickle
from itertools import groupby
from tqdm import tqdm as tqdm

'''
@Author: Daniel Sprague
@Lab: Calabrese Lab
@Department: Department of Pharmacology, University of North Carolina at Chapel Hill

Python 3.7
Anaconda Distribution

'''


class SEEKR(object):

    '''
    Implements the core essentials of SEEKR analyis for use in more specialized
    analysis or as a standalone tool

    Parameters
    ----------

    fasta_file: fasta formatted sequences of interest
    k: value of k for the analysis, typically k = [4,6] default = 5
    reference: reference set of sequences to standardize k-mer counts, defaults
    to included dataset of gencode lncRNA files, but optionally can
    use your own reference set


    Attributes
    ----------

    kmer_profile: kmer profile of sequences of interest
    reference: matrix of lenght normalized k-mer counts from reference file
    bases: sequence alphabet (for now only ATCG)
    keys: list of all possible kmers for given value of k
    seqnames: names of sequences supplied

    '''

    def __init__(self,fasta_file,reference,k):
        self.fasta_file = fasta_file
        self.fasta = far.Reader(fasta_file)
        self.k = k
        self.bases = ['A','T','C','G']
        self.keys = [''.join(p) for p in itertools.product(self.bases,repeat=self.k)]
        if reference[:-1] is 'p':
            self.reference = pickle.load(open(reference,'rb'))
        else:
            self.reference = self.generate_ref(reference)

        self.kmer_profile = self.kmer_profile()
        self.seqnames = self.fasta.get_headers()

    '''
    generate a dictionary of kmers as keys and '0' as the value
    '''
    def gen_kmersdict(self):
        return dict(zip(self.keys,np.zeros(4**self.k)))

    '''
    Calculate length normalized, standardized k-mer counts for sequences
    of interest
    '''
    def kmer_profile(self):

        kmerprofile_dict = {}
        ref_kmers = [v for i,v in self.reference.items() if f'{self.k}mer' in i]
        seqs,headers = self.fasta.get_seqs(), self.fasta.get_headers()
        arr = (kmerprofilelib.target_norm(ref_kmers[0],seqs,self.k))
        ''' translate and log transform
        '''
        arr = np.log2(arr+np.abs(np.min(arr))+1)

        for seq_name,row in enumerate(arr):
            kmerprofile_dict[f'{headers[seq_name]}_{self.k}mers'] = row
        return kmerprofile_dict

    '''
    Generate the reference set of length normalized k-mer counts for use in
    standardization

    matrix of size nxm where n is the number of sequences in the reference
    and m is 4^k

    Default reference set of all mouse long non-coding RNAs is supplied
    '''

    def generate_ref(self,reference):
        seqs,ref,pos,arr = far.Reader(reference).get_seqs(),{},4**self.k,[]

        d = collections.OrderedDict(zip(self.keys,range(0,pos)))
        for i, seq in enumerate(seqs):
            arr.append(kmerprofilelib.count_kmers(seq,self.k,d))
        arr = np.array(arr)
        ref[f'ref_{self.k}mers'] = arr
        return ref

    '''
    Calculate pairwise pearon's correlations within a dictionary containing > 1
    k-mer profiles

    Input: Dictionary of k-mer profiles, or the kmer_profile attribute of a SEEKR object
    Output: pandas dataframe
    '''
    def correlate_kmerprofiles(self,savename=None):
        df = pd.DataFrame.from_dict(self.kmer_profile)
        corr = df.corr()
        if savename:
            corr.to_csv(f'./{savename}.csv')
        return corr

    def save_kmer_profile(self,savename='kmer_profile'):
        pickle.dump(self.kmer_profile,open(f'./{savename}.p','wb'))

    def save_ref(self,savename='ref'):
        pickle.dump(self.reference,open(f'./{savename}.p','wb'))
