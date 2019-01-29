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

'''

class proteinSEEKR(SEEKR):

    '''
    Implements the SEEKR sequence analysis to predict protein binding to a
    DNA or RNA sequence

    Parameters
    ----------
    motif_path: path to motif file in matrix format. Columns are nucleotides
    and rows are positions within the motif
    fasta_file: fasta formatted sequences of interest
    k: value of k for the analysis, typically k = [4,6] default = 5
    reference: reference set of sequences to standardize k-mer counts, defaults
    to included dataset of gencode lncRNA files, but optionally can
    use your own reference set


    Attributes
    ----------
    pwms: dictionary of PWMs
    Inherited:
    SEEKR.kmer_profile: kmer profile of sequences of interest
    SEEKR.reference: matrix of lenght normalized k-mer counts from reference file
    SEEKR.bases: sequence alphabet (for now only ATCG)
    SEEKR.keys: list of all possible kmers for given value of k
    SEEKR.seqnames: names of sequences supplied

    '''

    def __init__(self,motif_path,fasta_file,reference,k):
        super(proteinSEEKR,self).__init__(fasta_file,reference,k)

        self.motif_path = motif_path
        self.pwms = self.read_motifs()
        self.kmer_profile = SEEKR.kmer_profile(self,log=False)


    '''read in motifs '''

    def read_motifs(self):
        pwms = [f for f in glob.iglob(self.motif_path+'*') if 'txt' in f]
        pwm_dict = collections.defaultdict(list)

        for file in pwms:
            try:
                '''is PWM file populated or empty'''
                pwm = pd.read_csv(file,sep='\t')
            except:
                continue
            pwm.drop('Pos',axis=1,inplace=True)
            lenmotif = pwm.shape[0]
            '''RNA to DNA alphabet'''
            pwm.rename(columns={'U':'T'},inplace=True)
            pwm_dict[os.path.basename(file)] = pwm
        return pwm_dict

    'return a dictionary of k-mers mapped to an index'
    def gen_kmersdict(self):
        return dict(zip(self.keys,np.zeros(4**self.k)))

    '''

    calculate the probability of finding a given kmer within a protein
    binding motif

    i = current index within motif
    j = window within motif (example: a 7 base pair motif has 3 possible 5mers)
    c = column of PWM, with values A,T,C, or G

    '''
    def calc_prob(self,p,frame,pwm,kmer):
        for current_position,nucleotide in enumerate(kmer):
            if p is None:
                p = self.pwms[pwm].iloc[current_position+frame][nucleotide]
            else:
                p*=self.pwms[pwm].iloc[current_position+frame][nucleotide]
        return p

    '''for all PWMs, calculate probabilities and store in dictionary'''

    def kmer_probabilities(self):
        pwm_dict = collections.defaultdict(list)

        for pwm in tqdm(self.pwms,total=len(self.pwms)):
            lenmotif = self.pwms[pwm].shape[0]

            if len(self.pwms[pwm]) < self.k:
                '''
                is the length of the motif less than the value of k?
                '''
                d = self.gen_kmersdict()
                kmers_within_kmer = [([kmer[i:i+(4)] for i in range(self.k-(4)+1)],kmer) for kmer in self.keys]
                '''
                kmers_within_kmer
                -----------------
                Find all 4-mers within the larger k-mer
                example error: 5-mer ATCGT in a 4 base pair motif

                This 5-mer does not exist within the 4 base pair motif,
                but the two 4-mers within the 5-mer,
                ATCG and TCGT can fit within the 4 base pair motif

                Calculate the probabilities of observing the 4mers
                seperately and then sum results

                No motif is < 4 base pairs, hence default of k = 4
                '''
                nkmers = lenmotif - 4 + 1
                frames = list(range(nkmers)) # possible starting positions of a k-mer in a motif

                for sub_kmers,kmer in kmers_within_kmer:
                    for sub_kmer in sub_kmers:
                        for frame in frames:
                            p = self.calc_prob(None,frame,pwm,sub_kmer)
                            d[kmer]+=p
                pwm_dict[pwm] = d.copy()
                d.clear()
            else:
                '''
                Otherwise, simply calculate the probability of the k-mer occuring
                in the motif
                '''
                d = self.gen_kmersdict()
                for kmer in self.keys:
                    nkmers = lenmotif - self.k + 1
                    frames = list(range(nkmers))
                    for frame in frames:
                        p = self.calc_prob(None,j,pwm,kmer)
                        d[kmer] += p
                pwm_dict[pwm] = d.copy()
                d.clear()
        return pwm_dict

    '''
    saves calculated probabilities dictionary to a pickle file
    '''
    def save_kmerprobs(self,kmerprobdict,savename):
        pickle.dump(kmerprobdict,open(f'{savename}.p','wb'))
        return

    '''
    get the k-mer zscores from the k-mer profile as a dataframe and re-index
    by the k-mer names (index = (AAAA,AAAT,AAAC,AAAG,...,GGGG))
    '''
    def get_zscore_df(self):
        df = pd.DataFrame.from_dict(self.kmer_profile)
        df['kmers'] = self.keys
        return df.set_index('kmers')

    '''
    weight the k-mer z-scores by the probability of finding the k-mer within
    a protein binding motif
    '''
    def kmer_weights(self,probabilities):
        z_scores = self.get_zscore_df()
        sorted_weights = np.array([kmer2weight[k] for k in z_scores.columns])
        weighted_z_scores = z_scores.values.copy() * sorted_weights
        scores_sums = weighted_z_scores.sum(axis=1)
        score_dict[pwm_path.name] = scores_sums
        return scoredict

    '''
    convert the weighted z score dictionary to a dataframe and rename columns
    to the name of the sequences
    '''
    def get_scoredataframe(self,scoredict):
        df = pd.DataFrame.from_dict(scoredict,orient='index')
        df.columns = self.seqnames
        return df

    '''
    map the proteins that recognize a given motif to the correct motif within
    the dataframe and append as a column
    '''
    def map_proteins(self,path,df):
        motifs = pd.read_csv(path,sep='\t')
        trimmed_motifs=motifs[motifs['Motif_ID']!='.']

        sorted_trimmed_motifs = trimmed_motifs.sort_values(by='Motif_ID')
        proteins_list = []
        for idx, row in df.iterrows():
            motif_name = idx[:-4]
            proteins = list(sorted_trimmed_motifs[sorted_trimmed_motifs['Motif_ID']==motif_name]['RBP_Name'])
            proteins_list.append(proteins)
        return proteins_list

    '''
    remove motifs with no associated proteins
    '''
    def prune_empty_motifs(self,df,proteins_list):
        df['Proteins']= proteins_list
        df_prune = df[df['Proteins'].str.len() > 0]
        return df_prune

    '''
    Example code for a clustered heatmap. Designed for small samples (n sequences
    ~ 10). Alternative visualizations best for larger n
    '''
    # def clusterplot(self,df,path='./plot.pdf'):
    #     cm = sns.clustermap(df,
    #               cmap='coolwarm',vmax=10,
    #               yticklabels=df[df['Proteins'].str.len() > 0]['Proteins'],
    #               figsize=(10,50),
    #               method='complete',linewidths=.5,annot=True,
    #               metric='euclidean')
    #     hm = cm.ax_heatmap.get_position()
    #     plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=10)
    #     cm.ax_heatmap.set_position([hm.x0, hm.y0, hm.width, hm.height])
    #     col = cm.ax_col_dendrogram.get_position()
    #     cm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*0.2])
    #     cm.cax.set_visible(False)
    #     cm.cax.remove()
    #     plt.savefig(path,bbox_inches='tight')
    #     return
