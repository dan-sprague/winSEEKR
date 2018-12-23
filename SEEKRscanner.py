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

    def dseekr(sd,tiled,l,s,plot_dict,lncref,queryfiles,plotrefs,k,fa_files):
        #Tile the currently open fa file
        queryfiles_kmers = dict([(i,v) for i,v in queryfiles.items() if f'{k}mer' in i]) # extract query kmer set for appropriate value of k
        queryfiles_kmers = dict(sorted(queryfiles_kmers.items())) # "sort" the keys alphabetically
        lncref_kmers = [v for i,v in lncref.items() if f'{k}mer' in i] #extract the reference set for the approopriate value of k
        if not tiled:
            curr_tile_fa = tile_seq(fa_files,l,s)
        elif tiled:
            infile = open(fa_files)
            curr_tile_fa = [line.upper() for line in infile]
        #Calculate normalized 4,5,6mers
        curr_normcount = target_norm(lncref_kmers[0],curr_tile_fa,k) # calculate standardized and length norm'd k-mer counts
        for query in queryfiles_kmers:
            print(query)
            lncomedata = plotrefs[plot_dict[os.path.basename(query)]]
            curr_q = queryfiles_kmers[query]
            R = kmer_pearson(curr_q,curr_normcount)
            make_plot(lncomedata, R,
            f'{os.path.basename(query)}+{os.path.basename(fa_files)}',f'{os.path.basename(fa_files)}',
            f'{os.path.basename(query)}_{os.path.basename(fa_files)}_plot',sd)
            np.savetxt(f'{os.path.basename(query)}_{os.path.basename(fa_files)}_data.txt',R,delimiter=',')

    def make_plot(lncref,R,title,xtitle,savename,sd):
        sns.set_context(context='talk',font_scale=2)
        R = pd.DataFrame(R,columns=['Rval'])
        x,y = [None,None]
        lncref = pd.DataFrame(lncref,columns=['Rval'])
        # Points that pass a given threshold of 'significance'
        rpass = [(R.index[R['Rval']==point][0],point) for point in R['Rval'] if point >= lncref['Rval'].mean()+(lncref['Rval'].std()*sd)]
        if len(rpass) > 0: # is there anything to plot?
            x,y = zip(*rpass)
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,4.8),sharey=True,gridspec_kw = {'width_ratios':[1, 5]})
        # Plot reference KDE/Histogram next to line plot
        sns.distplot(lncref[lncref['Rval']<.3],vertical=True,ax=ax1,label='mm10 lncome',color='C1')
        ax1.grid(False)
        ax1.invert_xaxis()
        # Scatter plot of values that pass threshold
        ax2.scatter(x,y,label=xtitle,color='C1')
        # Plot all points along transcript
        ax2.plot(R['Rval'])
        ax2.set_ylim(-.1,R['Rval'].max()+.06)

        ### plot lines correspond to standard deviations of Pearons comparisons in reference set
        ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval']),linestyle='--',color='orange')
        ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval'])*2,linestyle='--',color='orange')
        ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval'])*3,linestyle='--',color='orange')


        ax2.set_xlim([0, len(R)])
        # Below code is a for a specific experiment done earlier, remove when developing further
        locs = list(ax2.get_xticks())
        locs+=[34,77,198,207,281,292,348,991,1196,1272]
        ax2.set_xticks(list(np.arange(0,len(R),500))+locs)
        ax2.set_xticklabels(labels=[i/10 for i in np.arange(0,len(R),500)]+locs,rotation=90)
        #ax3 = ax2.twinx()
        #ax3.scatter(x=data['target'],y=data['Bit'],s=10*data['E']**2,color='r',label='nhmmer')
        #ax3.set(ylabel='Normalized Bit Score')

        ax1.set(ylabel='')
        plt.title('')
        #plt.tight_layout()
       # fig.legend(bbox_to_anchor=(.1, .9), loc=2, borderaxespad=0,mode='expand')

        plt.xlabel('')
        plt.savefig(f'{savename}_{sd}sd.pdf')
        plt.close()
        return
