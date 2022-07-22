import sys
from functools import partial

from Bio import SeqIO, bgzf
import regex
from collections import OrderedDict

import gzip
import argparse

import time

####################### HELPERS ############################

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def rc(seq):
	""" returns reverse complement of sequence """
	return "".join(complements.get(base, base) for base in reversed(seq))


def OR(RE_LIST,ENHANCEDMATCH=False):
	""" regex capture one or none in list of REs """
	if ENHANCEDMATCH:
		return "(?e)(" + "|".join(["(?:"+RE+")" for RE in RE_LIST]) +")"
	else:
		return "(" + "|".join(["(?:"+RE+")" for RE in RE_LIST]) +")"


def make_corrector(options, max_levenshtein_dist = 1):
	error_string = "{e<=%d}" % max_levenshtein_dist
	checkers = [(v, regex.compile("("+k+")"+error_string)) for k,v in options.items()]
	def corrector(match):
		for (v,c) in checkers:
			if c.match(match):
				return v
	return corrector


####################### SETTINGS ############################


fwd_primer = "TCGTCGGCAGCGTCAGTTGTATCAACTCAGAC"

rev_primers = { 
				"IGHA_HV" :	"GTCTCGTGGGCTCGGGGGGAAGAAGCCCTGGAC", 
				"IGHG_HV" :	"GTCTCGTGGGCTCGGGGGAAGTAGTCCTTGACCA", 
				"IGHM_HV" :	"GTCTCGTGGGCTCGGGAAGGAAGTCCTGTGCGAG", 
				"IGHE_HV" :	"GTCTCGTGGGCTCGGAAGTAGCCCGTGGCCAGG", 
				"IGHD_HV" :	"GTCTCGTGGGCTCGGTGGGTGGTACCCAGTTATCAA", 
				"mIGHG1" : "GTCTCGTGGGCTCGGCAGCAGATCCAGGGGCCAGTG"
				}

# reverse primers will be found in rev orientation
RevPrimerRCLookup = {v : k for k, v in rev_primers.items()}
RevPrimerLookup = {rc(v) : k for k,v in rev_primers.items()}
correct_rev_primer = make_corrector(RevPrimerLookup, max_levenshtein_dist = 3)

# in cleaned pacbio read configuration, i5 index is now in fwd direction
# I5Lookup = { v : k for k,v in i5_indices.items() if k in used_i5_indices}
# correct_i5_index = make_corrector(I5Lookup, max_levenshtein_dist = 2)

umi = "(.{4})(?:T)(.{4})(?:T)(.{4})(?:TCTCG)(.{2})(?:T)(.{2})(?:AT)"

# quickly check that read has plausible forward primer and reverse primer
fwd_pattern = fwd_primer + "{e<=3}" + umi 
rev_pattern = OR(RevPrimerRCLookup.keys()) + "{e<=3}"

# full_match_pattern = OR(list(I5Lookup.keys()),ENHANCEDMATCH = True) + "{e<=2}" \
full_match_pattern = "(?e)(" + fwd_primer + "){s<=3}" \
			 + "(" + umi + ")" \
			 + "([A,C,T,G,N]*)" \
			 + OR(list(RevPrimerLookup.keys()),ENHANCEDMATCH = True) + "{e<=3}" 

fwd_pattern_obj = regex.compile(fwd_pattern)
rev_pattern_obj = regex.compile(rev_pattern)

full_pattern_obj = regex.compile(full_match_pattern)

####################### I/O ############################


parser = argparse.ArgumentParser(description='Parses simple, non-concatemeric Pacbio hifi reads. Expects reads to be in correct orientation.')
parser.add_argument('input_fasta', metavar='input_path', type=str,
                    help='path to input fasta file')
parser.add_argument('-outdir', type = str, default = '.', 
					help = 'output directory '
							'(default: working directory)')
args = parser.parse_args()


## parse i/o paths
INPUT_FILENAME = args.input_fasta
OUTDIR = args.outdir

#set up output stream for parsed reads and reads that fail to parse
SAMPLE_NAME = INPUT_FILENAME.split("/")[-1].split(".fa")[0]

if SAMPLE_NAME.endswith("_simple"):
	SAMPLE_NAME = SAMPLE_NAME.split("_simple")[0]

parsed_out = gzip.open("{path}/{sample}_reattempted_passed.tsv.gz".format(path = OUTDIR, 
																	   sample = SAMPLE_NAME), 'wb')


#guess encoding of input file based on extension
_open = partial(gzip.open, mode = 'rt') if INPUT_FILENAME.endswith('.gz') else open


COLUMNS = ["seq_id", 				
		   "rev_comp",
			"5prime_seq", 
			# "i5_index_seq",
			# "i5_index_id", 
			"fwd_primer", 
			"umi_segment",
			"umi_seq", 
			"amplicon", 
			"rev_primer_seq", 
			"rev_primer_id", 
			"3prime_seq"]

header = "\t".join(COLUMNS) + "\n"

parsed_out.write(header.encode('utf-8'))

#set up counters for log
total = 0
has_primers = 0
good = 0

start = time.time()
with _open(INPUT_FILENAME) as f:
	records = SeqIO.parse(f, 'fasta')

	total = 0
	for record in records:

		sequence = str(record.seq) + str(record.seq)

		total += 1
		if total%10**3 ==0:
			sys.stderr.write("[%ds] Processed %d records...(%d matches)\n" % ((time.time() - start), total, good))
		
		has_i5_fwd_primer_and_umi = fwd_pattern_obj.search(sequence)
		has_rev_primer = rev_pattern_obj.search(rc(sequence))

		if has_i5_fwd_primer_and_umi and has_rev_primer:

			has_primers += 1

			match_obj = full_pattern_obj.search(sequence)

			if match_obj:
				good += 1
				row=OrderedDict([ (x,"") for x in COLUMNS])
		
				row["seq_id"] = record.id
				row["rev_comp"] = "(rc)" in record.id

				match_begin = match_obj.span()[0]
				match_end = match_obj.span()[1]

				row["5prime_seq"] = sequence[:match_begin]

				row["fwd_primer"] = match_obj.groups()[0]
				row["umi_segment"] = match_obj.groups()[1]
				row["umi_seq"] = ".".join(match_obj.groups()[2:7])
				row["amplicon"] = match_obj.groups()[7]
				row["rev_primer_seq"] = match_obj.groups()[8]
				row["3prime_seq"] = sequence[match_end:]


				# row["i5_index_id"] = correct_i5_index(row["i5_index_seq"])
				row["rev_primer_id"] = correct_rev_primer(row["rev_primer_seq"])

				row = "\t".join([str(x) for x in row.values()]) + "\n"
				parsed_out.write(row.encode('utf-8'))

				continue



parsed_out.close()

sys.stderr.write("############################# \n")

sys.stderr.write("# Finished processing!\n")
sys.stderr.write("# Total reads processed: %d \n" % total)
sys.stderr.write("# Have fwd and rev primers: %d [%.2f] \n"
					 % (has_primers, float((has_primers)/total)))
sys.stderr.write("# Passed parser: %d [%.2f] \n"
					 % (good, float((good)/total)))

sys.stderr.write("############################# \n")
