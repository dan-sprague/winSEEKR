import glob
import collections
import itertools
import pandas as pd
import os
import numpy as np
import kmerprofilelib

import fasta_reader as far
import numpy as np
import pickle
from itertools import groupby
from tqdm import tqdm as tqdm


from SEEKR import SEEKR

'''
@Author: Daniel Sprague
@Lab: Calabrese Lab
@Department: Department of Pharmacology, University of North Carolina at Chapel Hill

Python 3.7
Anaconda Distribution

'''


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
    SEEKR.kmer_profile: kmer profile of queries
    SEEKR.reference: matrix of lenght normalized k-mer counts from reference file
    SEEKR.bases: sequence alphabet (for now only ATCG)
    SEEKR.keys: list of all possible kmers for given value of k
    SEEKR.fa: parsed fasta file

    '''

    def __init__(self,queries,fasta_file,reference,k,windowLength,slide,threshold):
        super(SEEKRscanner,self).__init__(queries,reference,k)

        self.fa = far.Reader(fasta_file)
        self.windowLength = windowLength
        self.slide = slide
        self.treshold = threshold

    def scan(self):
        seqs,headers = self.fa.get_seqs(), self.fa.get_headers()
        ref = [v for i,v in self.reference.items() if f'{self.k}mer' in i][0]
        for seq in seqs:
            tiles = kmerprofilelib.tile_seq(seq,self.windowLength,self.slide)
            normcount_arr = kmerprofilelib.target_norm(ref,tiles,self.k)
            normcount_arr = np.log2(normcount_arr + np.abs(np.min(normcount_arr))+1)
            R = np.zeros((len(self.kmer_profile),len(tiles)))
            for i,(query,profile) in enumerate(self.kmer_profile.items()):
                R[i] = kmerprofilelib.kmer_pearson(profile,normcount_arr)
            df = pd.DataFrame(R.T)
            df.columns = self.seqnames
        return df

    def querydist(self):
        ref = [v for i,v in self.reference.items() if f'{self.k}mer' in i][0]
        norm_ref = (ref-np.mean(ref,axis=0))/np.std(ref,axis=0)
        norm_ref_log = np.log2(norm_ref + np.abs(np.min(norm_ref))+1)
        R = np.zeros((len(self.kmer_profile),np.size(norm_ref_log,0)))

        for i,(query,profile) in enumerate(self.kmer_profile.items()):
            R[i]=kmerprofilelib.kmer_pearson(profile,norm_ref_log)
        df = pd.DataFrame(R.T)
        df.columns = self.seqnames
        return df




    # def make_plot(lncref,R,title,xtitle,savename,sd):
    #     sns.set_context(context='talk',font_scale=2)
    #     R = pd.DataFrame(R,columns=['Rval'])
    #     x,y = [None,None]
    #     lncref = pd.DataFrame(lncref,columns=['Rval'])
    #     # Points that pass a given threshold of 'significance'
    #     rpass = [(R.index[R['Rval']==point][0],point) for point in R['Rval'] if point >= lncref['Rval'].mean()+(lncref['Rval'].std()*sd)]
    #     if len(rpass) > 0: # is there anything to plot?
    #         x,y = zip(*rpass)
    #     fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,4.8),sharey=True,gridspec_kw = {'width_ratios':[1, 5]})
    #     # Plot reference KDE/Histogram next to line plot
    #     sns.distplot(lncref[lncref['Rval']<.3],vertical=True,ax=ax1,label='mm10 lncome',color='C1')
    #     ax1.grid(False)
    #     ax1.invert_xaxis()
    #     # Scatter plot of values that pass threshold
    #     ax2.scatter(x,y,label=xtitle,color='C1')
    #     # Plot all points along transcript
    #     ax2.plot(R['Rval'])
    #     ax2.set_ylim(-.1,R['Rval'].max()+.06)
    #
    #     ### plot lines correspond to standard deviations of Pearons comparisons in reference set
    #     ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval']),linestyle='--',color='orange')
    #     ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval'])*2,linestyle='--',color='orange')
    #     ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval'])*3,linestyle='--',color='orange')
    #
    #
    #     ax2.set_xlim([0, len(R)])
    #     # Below code is a for a specific experiment done earlier, remove when developing further
    #     locs = list(ax2.get_xticks())
    #     locs+=[34,77,198,207,281,292,348,991,1196,1272]
    #     ax2.set_xticks(list(np.arange(0,len(R),500))+locs)
    #     ax2.set_xticklabels(labels=[i/10 for i in np.arange(0,len(R),500)]+locs,rotation=90)
    #     #ax3 = ax2.twinx()
    #     #ax3.scatter(x=data['target'],y=data['Bit'],s=10*data['E']**2,color='r',label='nhmmer')
    #     #ax3.set(ylabel='Normalized Bit Score')
    #
    #     ax1.set(ylabel='')
    #     plt.title('')
    #     #plt.tight_layout()
    #    # fig.legend(bbox_to_anchor=(.1, .9), loc=2, borderaxespad=0,mode='expand')
    #
    #     plt.xlabel('')
    #     plt.savefig(f'{savename}_{sd}sd.pdf')
    #     plt.close()
    #     return
