from SEEKR import SEEKR
import argparse
import pickle
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-fa',type=str)
parser.add_argument('-ref',type=str)
parser.add_argument('-k',type=int)
parser.add_argument('--savecorr',type=str,default='corrs')
parser.add_argument('--test',action='store_true')
args = parser.parse_args()

if args.test:
    print('Running Tests...\n')
    test = SEEKR('./Test_Data/test_kmers.fa','./Test_Data/test_reference_gencodelncrna.fa',5)
    test_raw = pickle.load(open('./Test_Data/test_rawkmers.p','rb'))
    test_profile = pickle.load(open('./Test_Data/test_kmerprofile.p','rb'))
    test_refmean = pickle.load(open('./Test_Data/test_refmean.p','rb'))
    test_refstd = pickle.load(open('./Test_Data/test_refstd.p','rb'))
    test_corr = pickle.load(open('./Test_Data/test_corr.p','rb'))

    print('Raw Kmer Counts: ',
    all(x==y for x,y in zip(test.raw_kmercounts(),test_raw)))

    print('Reference Mean: ',
    np.mean(test.reference['ref_5mers'],axis=0).all() == test_refmean.all())

    print('Reference Standard Deviation: ',
    np.std(test.reference['ref_5mers'],axis=0).all() == test_refstd.all())

    print('Kmer Profile: ',
    all(x==y for x,y in zip(test.kmer_profile,test_profile)))

    print('Correlation Matrices: ',
    test.correlate_kmerprofiles().equals(test_corr))

    exit()

SEEKRobj = SEEKR(args.fa,args.ref,args.k)


print(SEEKRobj.correlate_kmerprofiles(savename=args.savecorr))
