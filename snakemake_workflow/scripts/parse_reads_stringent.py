import sys
from functools import partial

from Bio import SeqIO, bgzf
import regex
from collections import OrderedDict

import gzip
import argparse

import time
import yaml

from pacbio_vdj_utils.sequence_utils import *

####################### I/O ############################


parser = argparse.ArgumentParser(description='Parses simple, non-concatemeric Pacbio hifi reads. Expects reads to be in correct orientation.')
parser.add_argument('input_path', metavar='input_path', type=str,
                    help='path to input sequence file (fasta or fastq)')
parser.add_argument('-config', type=str, required=True,
                    help='path to config file containing primer sequences')
parser.add_argument('-outdir', type = str, default = '.', 
					help = 'output directory '
							'(default: working directory)')
parser.add_argument('--amplicon', dest='amplicon', action='store_true',
					help='amplicons')
parser.add_argument('--fulllength', dest='amplicon', action='store_false', 
					help='full length transcripts')
parser.add_argument('--minibulk', dest='minibulk', action='store_true',
					help='minibulk')
parser.set_defaults(amplicon=True, minibulk=False)
args = parser.parse_args()


## parse i/o paths
INPUT_FILENAME = args.input_path
OUTDIR = args.outdir
CONFIG = args.config
AMPLICON = args.amplicon
MINIBULK = args.minibulk

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
		

failed_out = bgzf.BgzfWriter("{path}/{sample}_failed-parse.{fmt}.gz".format(path=OUTDIR, 
																	   	    sample=SAMPLE_NAME,
																	   	    fmt=seq_format), 
																		    "wb")
parsed_out = gzip.open("{path}/{sample}_parsed.tsv.gz".format(path=OUTDIR, 
															 sample=SAMPLE_NAME), 'wb')
parsed_fastqa_out = bgzf.BgzfWriter("{path}/{sample}_parsed.{fmt}.gz".format(path=OUTDIR, 
																	   sample=SAMPLE_NAME,
																	   fmt=seq_format), 'wb')

#######################  PARSER SETUP  ############################

with open(CONFIG) as configfile:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    read_config = yaml.load(configfile, Loader=yaml.FullLoader)

fwd_primer = read_config['fwd_primer']
rev_primers = read_config['rev_primers']
umi = read_config['umi_regex']
three_prime_umi = read_config['three_prime_umi_regex']

if MINIBULK:
	path_to_well_indices = read_config['well_indices']
	well_indices = fasta_to_dict(path_to_well_indices)
	well_indices = {k: rc(v) for k, v in well_indices.items()}
	well_index_lookup = {v: k for k, v in well_indices.items()}
	correct_well_index = make_corrector(well_index_lookup, max_dist=2)

if AMPLICON:
	pass
else:
	rev_primers = {k: v for k, v in rev_primers.items() if not("IG" in k)}

# reverse primers will be found in rev orientation
RevPrimerRCLookup = {v : k for k, v in rev_primers.items()}
RevPrimerLookup = {rc(v) : k for k,v in rev_primers.items()}
correct_rev_primer = make_corrector(RevPrimerLookup, max_dist=3)

# quickly check that read has plausible forward primer and reverse primer
fwd_pattern = fwd_primer + "{e<=3}" + umi # allows insertions, deletions, substitutions

rev_pattern = OR(RevPrimerRCLookup.keys()) + "{e<=3}"

if AMPLICON:
	# full_match_pattern = OR(list(I5Lookup.keys()),ENHANCEDMATCH = True) + "{e<=2}" \
	full_match_pattern = "(?e)(" + fwd_primer + "){s<=3}" \
				 + "(" + umi + ")" \
				 + "([A,C,T,G,N]*)" \
				 + OR(list(RevPrimerLookup.keys()),ENHANCEDMATCH=True) + "{e<=3}" 
	COLUMNS = ["seq_id", 				
			   "rev_comp",
				"5prime_seq", 
				"fwd_primer", 
				"umi_segment",
				"umi_seq", 
				"amplicon", 
				"rev_primer_seq", 
				"rev_primer_id", 
				"3prime_seq"]
else:
	full_match_pattern = "(?e)(" + fwd_primer + "){s<=3}" \
				 + "(" + umi + ")" \
				 + "([A,C,T,G,N]*)" \
				 + three_prime_umi \
				 + OR(list(RevPrimerLookup.keys()),ENHANCEDMATCH=True) + "{s<=3}" 

	COLUMNS = ["seq_id", 				
			   "rev_comp",
				"5prime_seq", 
				"fwd_primer", 
				"umi_segment",
				"umi_seq", 
				"amplicon",
				"3prime_umi_seq", 
				"rev_primer_seq", 
				"rev_primer_id", 
				"3prime_seq"]

if MINIBULK:
	fwd_pattern = fwd_pattern
	full_match_pattern = OR(well_index_lookup.keys(),ENHANCEDMATCH=True) + "{e<=2}" + full_match_pattern
	COLUMNS = ['well_barcode_seq', 'well_barcode_id'] + COLUMNS

fwd_pattern_obj = regex.compile(fwd_pattern)
rev_pattern_obj = regex.compile(rev_pattern)

full_pattern_obj = regex.compile(full_match_pattern)




####################### ######### ############################


header = "\t".join(COLUMNS) + "\n"

parsed_out.write(header.encode('utf-8'))

#set up counters for log
total = 0
has_primers = 0
good = 0

start = time.time()
with _open(INPUT_FILENAME) as f:
	records = SeqIO.parse(f, seq_format)

	total = 0

	for record in records:

		sequence = str(record.seq)

		total += 1
		if total%10**3 ==0:
			sys.stderr.write("[%ds] Processed %d records...(%d matches)\n" % ((time.time() - start), total, good))
		
		has_i5_fwd_primer_and_umi = fwd_pattern_obj.search(sequence)
		has_rev_primer = rev_pattern_obj.search(rc(sequence))
		#if (has_i5_fwd_primer_and_umi):
		#	print(has_i5_fwd_primer_and_umi)
		
		#if (has_rev_primer):
			#print(has_rev_primer)
		if has_i5_fwd_primer_and_umi and has_rev_primer:

			has_primers += 1

			# attempt to find a match for overall read pattern
			match_obj = full_pattern_obj.search(sequence)

			if match_obj:
				good += 1
				row=OrderedDict([ (x,"") for x in COLUMNS])
		
				row["seq_id"] = record.id
				row["rev_comp"] = "(rc)" in record.id

				match_begin = match_obj.span()[0]
				match_end = match_obj.span()[1]

				row["5prime_seq"] = sequence[:match_begin]

				if MINIBULK:
					row['well_barcode_seq'] = match_obj.groups()[0]
					row['well_barcode_id'] = correct_well_index(row['well_barcode_seq'])
					idx=1
				else:
					idx=0
				row["fwd_primer"] = match_obj.groups()[idx]
				row["umi_segment"] = match_obj.groups()[idx+1]
				row["umi_seq"] = "".join(match_obj.groups()[idx+2:idx+7])
				row["amplicon"] = match_obj.groups()[idx+7]
				if AMPLICON:
					row["rev_primer_seq"] = match_obj.groups()[idx+8]
				else:
					row['3prime_umi_seq'] = match_obj.groups()[idx+8]
					row["rev_primer_seq"] = match_obj.groups()[idx+9]
				row["3prime_seq"] = sequence[match_end:]

				row["rev_primer_id"] = correct_rev_primer(row["rev_primer_seq"])


				FASTA_RECORD_BARCODE = row["umi_seq"]
				if not AMPLICON:
					FASTA_RECORD_BARCODE += "_" + row["3prime_umi_seq"]
				# print trimmed record to fasta
				trimmed_record = record[match_obj.start(8):match_obj.end(8)]

				if MINIBULK:
					trimmed_record.id = "{}|BARCODE={}|CPRIMER={}".format(record.id,
																	 "_".join([FASTA_RECORD_BARCODE,
																	  row["rev_primer_id"]]),
																	  row["well_barcode_id"])
				else:
					trimmed_record.id = "{}|BARCODE={}|CPRIMER={}".format(record.id,
																	 FASTA_RECORD_BARCODE,
																	  row["rev_primer_id"])
				trimmed_record.name = None
				trimmed_record.description = ""

				row = "\t".join([str(x) for x in row.values()]) + "\n"
				parsed_out.write(row.encode('utf-8'))

				SeqIO.write(trimmed_record, handle=parsed_fastqa_out, format=seq_format)

				continue

		# if read failed to parse, write to failed read stream
		SeqIO.write(record, handle=failed_out, format=seq_format)

parsed_fastqa_out.close()
parsed_out.close()
failed_out.close()

sys.stderr.write("############################# \n")

sys.stderr.write("# Finished processing!\n")
sys.stderr.write("# Total reads processed: %d \n" % total)
sys.stderr.write("# Have fwd and rev primers: %d [%.2f] \n"
					 % (has_primers, float((has_primers)/total)))
sys.stderr.write("# Successfully parsed: %d [%.2f] \n"
					 % (good, float((good)/total)))

sys.stderr.write("############################# \n")
