import sys
import pandas as pd 


infile = sys.argv[1] # path/to/combined_counts.tsv
outfile = sys.argv[2] # path/to/combined_genecounts.tsv
#location of biomart export
biomart_Ensg2Name = sys.argv[3]



print ( "loading infile" )
df = pd.read_csv(infile, sep='\t', index_col = 0)
# Prepare infile for mapping gene symbols
df = df.T
df.reset_index(inplace = True)
convert_names_df = pd.read_csv(biomart_Ensg2Name)
gene_dict = dict(zip(convert_names_df['Gene stable ID'], convert_names_df['Gene name']))
print ( "converting gene names" ) 

df.loc[:,'index'] = df.loc[:,'index'].replace(gene_dict)
df.set_index('index', inplace = True)
# write outfile as .tsv

df.to_csv(outfile, sep = '\t')

