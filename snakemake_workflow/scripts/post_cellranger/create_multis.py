import pandas as pd
import sys
import argparse
import shutil

#####################################################    

parser = argparse.ArgumentParser(description='Creates multi_configs from 10X template using samplesheet')
parser.add_argument('-template', 
                    required=True, 
                    nargs='+',
                    help='path to multi config template')
parser.add_argument('-outdir', 
                    default='.', 
                    help="output directory (default: working directory)")
parser.add_argument('-samplesheet', 
                    required=True)
parser.add_argument('-raw',
                    help="path to raw fastqs from illumina run",
                    required=True)
args = parser.parse_args()
SAMPLESHEET = args.samplesheet
TEMPLATE = args.template
OUTDIR = args.outdir
RAW = args.raw
samplesheet = pd.read_table(SAMPLESHEET, sep = '    ')

for sample_uid in samplesheet.sample_id:
    # copy the template file and call it multi_config_{sampleuid}.csv
    original = TEMPLATE[0]
    target = "{}{}_{}.csv".format(OUTDIR, "multi_config", sample_uid)
    print("creating ", target)
    shutil.copyfile(original, target)
    # add new line which points to the libraries with the sample_uid
    new_line = "{}_TotalSeq,{}/{},any,{}_TotalSeq,Antibody Capture".format(sample_uid, RAW, samplesheet.donor.unique()[0], sample_uid)
    file_object = open(target, "a")
    file_object.write("{}".format(new_line))
    file_object.close()
