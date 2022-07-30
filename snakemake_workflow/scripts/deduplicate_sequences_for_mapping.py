import sys
from functools import partial

from Bio import SeqIO

import gzip
import argparse

import time


####################### I/O ############################


parser = argparse.ArgumentParser(description='Deduplicates sequences')
parser.add_argument('input_path', metavar='input_path', type=str,
                    help='path to input sequence file (fasta or fastq)')
parser.add_argument('-outdir', type = str, default = '.', 
					help = 'output directory '
							'(default: working directory)')
parser.add_argument('-copy_fields', nargs='+', help = 'list of fields in sequence id to retain')
args = parser.parse_args()


## parse i/o paths
INPUT_FILENAME = args.input_path
OUTDIR = args.outdir
COPY_FILEDS = args.copy_fields

#set up output stream for parsed reads and reads that fail to parse
SAMPLE_NAME = INPUT_FILENAME.split("/")[-1].split(".fa")[0]

if SAMPLE_NAME.endswith("_cleaned"):
	SAMPLE_NAME = SAMPLE_NAME.split("_cleaned")[0]

#guess encoding of input file based on extension
_open = partial(gzip.open, mode = 'rt') if INPUT_FILENAME.endswith('.gz') else open

#guess file format based on extension
if INPUT_FILENAME.strip(".gz").endswith('.fa') or INPUT_FILENAME.strip(".gz").endswith('fasta'):
	seq_format = "fasta"
elif INPUT_FILENAME.strip(".gz").endswith('.fq') or INPUT_FILENAME.strip(".gz").endswith('fastq'):
	seq_format = "fastq"
else:
	raise IOError("Sequence input file {} cannot be read "
				 "(Allowed formats: .fasta, .fa, .fastq, .fq and compressed versions)".format(INPUT_FILENAME))
		

outfile = open("{path}/{sample}_deduplicated.fasta".format(path = OUTDIR, 
															   sample = SAMPLE_NAME), 
															   "w")
seq_id_barcode_outfile = open("{path}/{sample}_seq_ids_barcodes.tsv".format(path = OUTDIR, 
															 sample = SAMPLE_NAME), 
															 "w")
seq_dict = {}

start = time.time()

def record_id_to_dict(record):
	field_list = [x.split("=") for x in record.id.split("|") if "=" in x]
	record_id_dict = {item[0] : item[1] for item in field_list}
	return record_id_dict

with _open(INPUT_FILENAME) as f:
	records = SeqIO.parse(f, seq_format)

	total = 0

	for record in records:
		sequence = str(record.seq)
		sequence_id_dict=record_id_to_dict(record)
		
		seq_dict[sequence] = seq_dict.get(sequence,{f:[] for f in COPY_FILEDS})
		for f in COPY_FILEDS:
			seq_dict[sequence][f].append(sequence_id_dict[f])

		total += 1
		

seq_id_barcode_outfile.write("SEQ_ID\tSEQ_ID_INFO\n")

seq_id = -1
for sequence, record_info in seq_dict.items():
	seq_id += 1

	seq_id_info = "|".join(["{f}={vals}".format(f=f,vals=",".join(vals)) 
					  for f, vals in record_info.items()])

	num_duplicates=len(list(record_info.values())[0])
	new_id = "SEQ_ID={}|DUPCOUNT={}".format(seq_id, num_duplicates)
	outfile.write(">{}\n".format(new_id))
	outfile.write("{}\n".format(sequence))
	seq_id_barcode_outfile.write("{}\t{}\n".format(seq_id, seq_id_info))

outfile.close()
seq_id_barcode_outfile.close()

sys.stderr.write("############################# \n")
sys.stderr.write("# [%ds] Processed %d records.(%d unique)\n" % ((time.time() - start), total, len(seq_dict)))
sys.stderr.write("############################# \n")


