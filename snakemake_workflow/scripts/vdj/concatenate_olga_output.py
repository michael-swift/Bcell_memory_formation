import sys
import argparse

import numpy as np
import pandas as pd

######################## PATH CONFIG ################################

parser = argparse.ArgumentParser()
parser.add_argument('-input', nargs='+',
                    required=True,
                    help="path to input directory")
parser.add_argument('-output', required="True",
                    help="path to output dataframe")


######################## PARSE ARGS ################################

args = parser.parse_args()

input_paths = args.input
output = args.output

######################################################################
pdf = []
for filename in sorted(input_paths):
    try:
        pdf.append(pd.read_table(filename, header=None))
    except EmptyDataError:
        print(f"File {filename} is empty", file=sys.stderr)
pdf = pd.concat(pdf)
pdf.columns = ['cdr3_full', 'pgen_nt', 'cdr3_aa_full', 'pgen_aa']
pdf = pdf.groupby(['cdr3_full', 'pgen_nt', 'cdr3_aa_full', 'pgen_aa']).size()
pdf = pdf.reset_index()[['cdr3_full', 'pgen_nt', 'cdr3_aa_full', 'pgen_aa']]

pdf.to_csv(f'{output}', sep='\t', header=False)
