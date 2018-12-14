import numpy as np
import itertools
import math
import pandas as pd

from itertools import groupby
from scipy.stats.stats import pearsonr

import sys
import fasta_reader
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
import os
import multiprocessing

import collections

'''

# ``````````````````````````````````````````````````````````````````````````````
# This library of functions performs all necessary operations for dseekr
# lncRNA analysis

# Members:
# count_kmers, target_norm, kmer_pearson, tile_seq, make_plot
# dseekr, global_tile_seq (now redundant), global_kmer_pearson (now redundant)
# global_stats


# ``````````````````````````````````````````````````````````````````````````````
'''

'''
# This function counts length normalized kmers per kb

# Input:
#   1. DNA sequence (str)
#   2. Value of k (int)
#   3. Dictionary of kmers (OrderedDict)

# Output:
#   1. Array of kmers for sequence (ndarray)
'''

def count_kmers(fa,k,d):
    currlen = len(fa)
    vals = np.zeros(len(d))
    for i in range(currlen-k+1):
        if fa[i:i+k] in d:
            vals[d[fa[i:i+k]]]+=1
    vals = 1000*(vals/currlen)
    return vals
'''
# This function normalizes kmer counts relative to the
# mm10 gencode lncrna reference for a DNA sequence
# or list of sequences

# Dependencies:
#   1. kmerprofilelib.count_kmers

# Input:
#   1. reference array of kmer counts (ndarray)
#   2. target sequence (string or list of strings)
#   3. kmer value (int)

# Output:
#   1. Normalized kmer counts (ndarray)
'''
def target_norm(ref,target,k):
    means,sd = np.mean(ref,axis=0),np.std(ref,axis=0)
    pos,bases = 4**k,['A','T','C','G']
    keys = [''.join(p) for p in itertools.product(bases,repeat=k)]
    d = collections.OrderedDict(zip(keys,range(0,pos)))
    if len(target) == 1:
        tilenorms = np.zeros(pos)
        target_kmer_count = count_kmers(target[0].upper(),k,d)
        target_kmer_count_norm = (target_kmer_count-means)/sd
        tilenorms = target_kmer_count_norm
    else:
        tilenorms = np.zeros((len(target),pos))
        j = 0
        for seq in target:
            target_kmer_count = count_kmers(seq.upper(),k,d)
            target_kmer_count_norm = (target_kmer_count-means)/sd
            tilenorms[j] = target_kmer_count_norm
            j+=1
    return tilenorms
'''
# This function calculates the pearson correlation
# coefficient between two "kmer arrays"

# Input:
#   1. Query normalized kmers (ndarray size 1x4^k)
#   2. Target normalized kmers (ndarray size num_tilesx4^k)

# Output
#   1. One dimensional array of pearson r values (ndarray)
#      of length 4^k
'''
def kmer_pearson(query,target):
    query = np.log2(query + np.abs(np.min(query)) + 1)
    target = np.log2(target + np.abs(np.min(target)) + 1)
    R,i = np.zeros(len(target)),0
    R[R==0] = np.nan
    for row in target:
        R[i] = pearsonr(query,row)[0]
        i+=1
    return R

''' tile_seq

opens a fasta file and returns a list of strings containing tiles of length L and interval s

input: path-like string, int length of tile, int interval
output: list of strings

'''

def tile_seq(fa,length,skip):
    reader = fasta_reader.Reader(fa)
    fa = [i.upper() for i in reader.get_seqs()]
    tiles = [fa[0][i:i+length] for i in range(0,len(fa[0]),skip)]
    return tiles

''' make_plot

saves a plot of the current dseekr instance

input:
    1. reference lncome
    2. pearson R array
    3. plot title
    4. x axis title
    5. savename
    6. threshold to call significant points (standard deviations)

output:
    1. pdf file of plot

'''

def make_plot(lncref,R,title,xtitle,savename,sd):
    sns.set_context(context='talk',font_scale=2)
    R = pd.DataFrame(R,columns=['Rval'])
    x,y = [None,None]
    lncref = pd.DataFrame(lncref,columns=['Rval'])
    rpass = [(R.index[R['Rval']==point][0],point) for point in R['Rval'] if point >= lncref['Rval'].mean()+(lncref['Rval'].std()*sd)]
    if len(rpass) > 0:
        x,y = zip(*rpass)
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,4.8),sharey=True,gridspec_kw = {'width_ratios':[1, 5]})
    sns.distplot(lncref[lncref['Rval']<.3],vertical=True,ax=ax1,label='mm10 lncome',color='C1')
    ax1.grid(False)
    ax1.invert_xaxis()
    ax2.scatter(x,y,label=xtitle,color='C1')
    ax2.plot(R['Rval'])
    ax2.set_ylim(-.1,R['Rval'].max()+.06)

    ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval']),linestyle='--',color='orange')
    ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval'])*2,linestyle='--',color='orange')
    ax2.axhline(np.mean(lncref['Rval'])+np.std(lncref['Rval'])*3,linestyle='--',color='orange')


    ax2.set_xlim([0, len(R)])

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

def dseekr(sd,tiled,l,s,plot_dict,lncref,queryfiles,plotrefs,k,fa_files):
    #Tile the currently open fa file
    queryfiles_kmers = dict([(i,v) for i,v in queryfiles.items() if f'{k}mer' in i])
    queryfiles_kmers = dict(sorted(queryfiles_kmers.items()))
    lncref_kmers = [v for i,v in lncref.items() if f'{k}mer' in i]
    if not tiled:
        curr_tile_fa = tile_seq(fa_files,l,s)
    elif tiled:
        infile = open(fa_files)
        curr_tile_fa = [line.upper() for line in infile]
    #Calculate normalized 4,5,6mers
    curr_normcount = target_norm(lncref_kmers[0],curr_tile_fa,k)
    for query in queryfiles_kmers:
        print(query)
        lncomedata = plotrefs[plot_dict[os.path.basename(query)]]
        curr_q = queryfiles_kmers[query]
        R = kmer_pearson(curr_q,curr_normcount)
        make_plot(lncomedata, R,
        f'{os.path.basename(query)}+{os.path.basename(fa_files)}',f'{os.path.basename(fa_files)}',
        f'{os.path.basename(query)}_{os.path.basename(fa_files)}_plot',sd)
        np.savetxt(f'{os.path.basename(query)}_{os.path.basename(fa_files)}_data.txt',R,delimiter=',')


def global_tile_seq(fa,length,skip):
    fa = fa.upper()
    tiles = [fa[i:i+length] for i in range(0,len(fa),skip)]
    return tiles

def global_kmer_pearson(query,target):
    query = np.log2(query + np.abs(np.min(query)) + 1)
    target = np.log2(target + np.abs(np.min(target)) + 1)
    R = np.zeros(len(target))
    R[R==0] = np.nan
    i = 0
    for row in target:
        R[i] = pearsonr(query,row)[0]
        i+=1
    return R


def global_stats(sd,l,s,plot_dict,lncref,queryfiles,plotrefs,lncrnas,out,k):
    queryfiles_kmers = dict([(i,v) for i,v in queryfiles.items() if f'{k}mer' in i])
    lncref_kmers = [v for i,v in lncref.items() if f'{k}mer' in i]


    reader = fasta_reader.Reader(lncrnas)
    seqs = reader.get_seqs()

    queryfiles_kmers = dict(sorted(queryfiles_kmers.items()))

    counts = np.zeros((len(seqs),len(queryfiles_kmers)))

    #MAKE ASSIGNMENT OF QUERY TO ARRAY INDEX EXPLICIT
    #USING DICT

    queryids = dict(zip(list(queryfiles_kmers),range(len(queryfiles_kmers))))

    for i,seq in enumerate(seqs):
        curr_tile_fa = global_tile_seq(seq,l,s)
        curr_normcount = target_norm(lncref_kmers[0],curr_tile_fa,k)

        for query in queryfiles_kmers:
            lncomedata = plotrefs[query]
            thresh = np.mean(lncomedata)+np.std(lncomedata)*sd
            curr_q = queryfiles_kmers[query]
            R = global_kmer_pearson(curr_q,curr_normcount)
            sig_hits = len(R[R>thresh])
            counts[i][queryids[query]] = sig_hits/len(curr_tile_fa)

    np.savetxt(f'./{out}_{k}mers_{sd}sd.csv',counts,delimiter=',')
    return
