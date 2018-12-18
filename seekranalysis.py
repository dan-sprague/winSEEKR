from SEEKR import SEEKR
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-fa',type=str)
parser.add_argument('-ref',type=str)
parser.add_argument('-k',type=int)
parser.add_argument('--savecorr',type=str,default='corrs')
args = parser.parse_args()

SEEKRobj = SEEKR(args.fa,args.ref,args.k)


print(SEEKRobj.correlate_kmerprofiles(savename=args.savecorr))
