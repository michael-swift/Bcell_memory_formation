import os
import pandas as pd
import argparse


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('infile_samplesheet')
    parser.add_argument('outfile')

    return parser.parse_args()


def main():
    """ Combines individual samtools depth files into a single table. """

    args = parse_args()

    df = pd.read_table(args.infile_samplesheet)

    counts = []
    for row in df.itertuples():
        infile = os.path.join(row.base, row.samplename, 'STAR_output',
                              'depth_igh_locus.txt.gz')

        try:
            tmpdf = pd.read_table(infile, header=None,
                                  names=['chrom', 'pos', row.samplename])
        except pd.errors.EmptyDataError:
            continue

        # interested in pd.Series where position is the index
        tmpdf = tmpdf.set_index('pos')[row.samplename]

        counts.append(tmpdf)

    # outer join because not all cells will have depth for all positions
    # memory consumption explodes:
    outdf = pd.concat(counts, axis=1, join='outer').fillna(0).astype(int)

    comp = None
    if args.outfile.endswith('.gz'):
        comp = 'gzip'

    outdf.to_csv(args.outfile, sep='\t', compression=comp)


if __name__ == "__main__":
    main()
