from proteinSEEKR import proteinSEEKR
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-fa',type=str)
parser.add_argument('-ref',type=str)
parser.add_argument('-k',type=int)
parser.add_argument('-pwms',type=str)
parser.add_argument('-pwmmeta',type=str)
parser.add_argument('--savedf',type=str,default='dataframe')
args = parser.parse_args()

ps = proteinSEEKR(args.pwms,args.fa,args.ref,args.k)

probs = ps.kmer_probabilities()
weights = ps.kmer_weights(probs)


df = ps.get_scoredataframe(weights)
proteins = ps.map_proteins(args.pwmmeta,df)

out_df = ps.prune_empty_motifs(df,proteins)

out_df.to_csv(f'./{args.savedf}.csv')
