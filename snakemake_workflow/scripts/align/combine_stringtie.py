import sys
import os
import glob
import pandas as pd

outputs_dir = sys.argv[1]
outfile = sys.argv[2] # output filename

glob.glob(outputs_dir + '*/*/*')

for row in df_samplesheet.itertuples():
    infile = os.path.join(row.base, row.samplename, 'STAR_output','ReadsPerGene.out.tab')
    _df = pd.read_csv(infile, sep = '\t')
    _df = _df.iloc[:,:2]
    _df.columns = ['ensg_feature', row.samplename]
    _df = _df.set_index('ensg_feature').T
    if first == True:
        _df.to_csv(outfile, sep ='\t')
        first = False
    else:
        _df.to_csv(outfile, sep = '\t', header = False, mode = 'a')
