import os
import pandas as pd
import argparse


def parse_args():

    parser = argparse.ArgumentParser()

    
    parser.add_argument('base_dir')
    parser.add_argument('out_file')
    parser.add_argument('infile_samplesheet', nargs='+')
    parser.add_argument('infile_name')
    return parser.parse_args()
 

def main():
    """ Combines individual star log files into a single table. """

    args = parse_args()
    
    #df = pd.read_table(args.infile_samplesheet, 
     #                       sep='\t', 
      #                      engine='python')
    samplesheet = pd.concat([pd.read_table(
    "{}".format(x), dtype='str')
    for x in args.infile_samplesheet], ignore_index=True)

    samplesheet['sample_uid'] = samplesheet['library_uid'] + "-" + samplesheet['sample_index']
    samplesheet.set_index('sample_uid', inplace=True)
    
    samples = samplesheet.index.to_list()

    # hardcoding this which I know is stupid but YOLO
    logs = ["{}/{}/{}".format(args.base_dir, bc, args.infile_name) for bc in samples]
    num_input_reads = []
    num_uniquely_mapped_reads = []
    num_reads_multiple_loci = []
    num_reads_too_many_loci = []
    percent_reads_unmapped_too_many_mismatches = []
    percent_reads_unmapped_too_short = []
    sep='/'
    print(logs)

    for log in logs:

        
        with open(log) as f:

            for line in f:

                if "Number of input reads" in line:
                    x = line.rstrip().split()[-1]
                    num_input_reads.append(x)

                if "Uniquely mapped reads number" in line:
                    x = line.rstrip().split()[-1]
                    num_uniquely_mapped_reads.append(x)
                if "Number of reads mapped to multiple loci" in line:
                    x = line.rstrip().split()[-1]
                    num_reads_multiple_loci.append(x)

                if "Number of reads mapped to too many loci" in line:
                    x = line.rstrip().split()[-1]
                    num_reads_too_many_loci.append(x)

                if "% of reads unmapped: too many mismatches" in line:
                    x = line.rstrip().split()[-1]
                    percent_reads_unmapped_too_many_mismatches.append(x)

                if "% of reads unmapped: too short" in line:
                    x = line.rstrip().split()[-1]
                    percent_reads_unmapped_too_short.append(x)


    features = ["input",
                "uniquely_mapped",
                "multiple_loci",
                "too_many_loci",
                "percent_unmapped_too_many_mismatches",
                "percent_unmapped_too_short"]


    stats = [num_input_reads, num_uniquely_mapped_reads, num_reads_multiple_loci,
             num_reads_too_many_loci, percent_reads_unmapped_too_many_mismatches, percent_reads_unmapped_too_short]

    with open(args.out_file, 'w') as out:

        header = ["name"] + features # header
        out.write("\t".join(header) + "\n")

        for i in range(len(samples)):
            sample = samples[i]
            my_stats = [stats[x][i] for x in range(len(stats))]
            line = "\t".join([sample] + my_stats)
            out.write(line + "\n")

if __name__ == "__main__":
    main()
