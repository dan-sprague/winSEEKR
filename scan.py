from SEEKRscanner import SEEKRscaner
import fasta_reader as far
import argparse
import pickle
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-fa',type=str)
parser.add_argument('-ref',type=str)
parser.add_argument('-k',type=int)
parser.add_argument('-w',default=1000,type=int)
parser.add_argument('-s',default=100,type=int)
parser.add_argument('--thresh',default=3,type=int)
parser.add_argument('--test',action='store_true')
args = parser.parse_args()


fa = far.Reader(args.fa)
headers,seqs = fa.get_headers(),fa.get_seqs()

for i,seq in enumerate(seqs):
    continue

    
