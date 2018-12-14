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


class SEEKR(object):

    def __init__(self,fasta_file,k=5,reference=None):
        self.fasta_file = fasta_file
        self.fasta = far.Reader(fasta_file)
        self.k = k
        self.bases = ['A','T','C','G']
        self.keys = [''.join(p) for p in itertools.product(self.bases,repeat=self.k)]
        if reference is None:
            self.reference = pickle.load(open('./ref/ref.p','rb'))
        else:
            self.reference = self.generate_ref(reference)

        self.kmer_profile = self.kmer_profile()
        self.seqnames = self.fasta.get_headers()

    def gen_kmersdict(self):
        return dict(zip(self.keys,np.zeros(4**self.k)))

    def kmer_profile(self):

        kmerprofile_dict = {}
        ref_kmers = [v for i,v in self.reference.items() if f'{self.k}mer' in i]
        seqs = self.fasta.get_seqs()
        headers = self.fasta.get_headers()
        arr = (kmerprofilelib.target_norm(ref_kmers[0],seqs,self.k))
        for seq_name,row in enumerate(arr):
            kmerprofile_dict[f'{headers[seq_name]}_5mers'] = row
        return kmerprofile_dict

    def generate_ref(self,reference):
        seqs,ref,pos,arr = far.Reader(reference).get_seqs(),{},4**self.k,[]

        d = collections.OrderedDict(zip(self.keys,range(0,pos)))
        for i, seq in enumerate(seqs):
            arr.append(kmerprofilelib.count_kmers(seq,self.k,d))
        arr = np.array(arr)
        ref[f'ref_{self.k}mers'] = arr
        return ref

    def save_kmer_profile(self,savename='kmer_profile'):
        pickle.dump(self.kmer_profile,open(f'./{savename}.p','wb'))

    def save_ref(self,savename='ref'):
        pickle.dump(self.reference,open(f'./{savename}.p','wb'))
