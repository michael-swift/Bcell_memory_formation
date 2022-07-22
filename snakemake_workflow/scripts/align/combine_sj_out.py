import os
import pandas as pd
import argparse


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('base_dir')
    parser.add_argument('out_file')
    parser.add_argument('infile_samplesheet')
    parser.add_argument('infile_name')
    return parser.parse_args()


def main():
    """ Combines individual STAR SJ.out.tab files into a single table. """

    args = parse_args()
    print(args.infile_samplesheet)
    df = pd.read_table(args.infile_samplesheet, sep = '\t', engine = 'python')
    samples = df.sample_index.to_list()
    files = ["{}/{}/{}".format(args.base_dir, bc, args.infile_name) for bc in samples]
    sj_cols = ['chrom', 'start','end', 'strand', 'intron_motif', 'annotated',
               'unique_mapping', 'multi_mapping', 'max_overhang']
    sjs = []
    for sample in samples:
        infile = os.path.join(args.base_dir, 
                              sample, 
                              args.infile_name)

        tmpdf = pd.read_table(infile, 
                              header=None,
                              names=sj_cols)
        
        tmpdf['sample'] = sample

        sjs.append(tmpdf)

    outdf = pd.concat(sjs)

    outdf.to_csv(args.out_file, sep='\t')
    
if __name__ == "__main__":
    main()
