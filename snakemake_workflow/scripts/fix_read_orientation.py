from Bio import SeqIO, bgzf
import regex
import sys
import gzip
from functools import partial

import argparse

####################### HELPERS ############################

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def rc(seq):
	""" returns reverse complement of sequence """
	return "".join(complements.get(base, base) for base in reversed(seq))
	

####################### SETTINGS ############################


FWD_PRIMER = "CTACACGACGCTCTTCCGATCT"
#FWD_PRIMER = "TCGTCGGCAGCGTCAGTTGTATCAACTCAGAC"
FWD_PRIMER_RC = rc(FWD_PRIMER)

fwd_pattern = regex.compile('(%s){e<=3}' % (FWD_PRIMER))


############################################################

parser = argparse.ArgumentParser(description='Enforces uniform orientation in Pacbio hifi reads and removes concatemers. Read orientation is inferred based on orientation of forward primer.')
parser.add_argument('input_file', metavar='input_path', type=str,
                    help='path to input file')
parser.add_argument('-outdir', type = str, default = '.', 
					help = 'output directory '
							'(default: working directory)')
args = parser.parse_args()


## parse i/o paths
INPUT_FILENAME = args.input_file
OUTDIR = args.outdir

# infer encoding of file based on extension
_open = partial(gzip.open, mode = 'rt') if INPUT_FILENAME.endswith('.gz') else open

# infer file format based on extension
if INPUT_FILENAME.strip(".gz").endswith('.fa') or INPUT_FILENAME.strip(".gz").endswith('fasta'):
	seq_format = "fasta"
elif INPUT_FILENAME.strip(".gz").endswith('.fq') or INPUT_FILENAME.strip(".gz").endswith('fastq'):
	seq_format = "fastq"
else:
	raise IOError("Sequence input file {} cannot be read "
				 "(Allowed formats: .fasta, .fa, .fastq, .fq and compressed versions)".format(INPUT_FILENAME))
		

#set up output streams
SAMPLE_NAME = INPUT_FILENAME.split("/")[-1].split(".fa")[0]

simple_out = bgzf.BgzfWriter("{path}/{sample}_cleaned.{fmt}.gz".format(path=OUTDIR, 
																	   sample=SAMPLE_NAME,
																	   fmt=seq_format), 
							"wb")
concatemeric_out = bgzf.BgzfWriter("{path}/{sample}_concatemer.{fmt}.gz".format(path=OUTDIR, 
																	   sample=SAMPLE_NAME,
																	   fmt=seq_format), 
							"wb")
double_fwd_out = bgzf.BgzfWriter("{path}/{sample}_double_fwd_primer.{fmt}.gz".format(path=OUTDIR, 
																	   sample=SAMPLE_NAME,
																	   fmt=seq_format), 
							"wb")
no_amp_handle_out = bgzf.BgzfWriter("{path}/{sample}_no_handle.{fmt}.gz".format(path=OUTDIR, 
																	   sample=SAMPLE_NAME,
																	   fmt=seq_format), 
							"wb")


#set up counters for log file
it = 0
simple_fwd_it = 0
simple_rev_it = 0
concatemeric_fwd_it = 0
concatemeric_rev_it = 0
double_it = 0
no_handle_it = 0

#classify reads based on orientation of forward primer
with _open(INPUT_FILENAME) as f:
	for record in SeqIO.parse(f,seq_format):
		it += 1

		if it % 10**4 == 0:
			sys.stderr.write('Processed {n} reads...\n'.format(n = it))

		seq = str(record.seq)
		num_fwd_primers = len(fwd_pattern.findall(seq))
		num_fwd_primers_in_rc = len(fwd_pattern.findall(rc(seq)))

		# reads with correct chemistry should have a single fwd primer
		# if it is in forward direction just print 
		if (num_fwd_primers == 1) and (num_fwd_primers_in_rc == 0):
			SeqIO.write(record, handle = simple_out, format=seq_format)
			simple_fwd_it += 1
			continue

		# otherwise reverse complement and print
		if (num_fwd_primers == 0) and (num_fwd_primers_in_rc == 1):
			rc_record = record.reverse_complement(id=record.id + "(rc)", 
												  name = None, description = "")
			SeqIO.write(rc_record, handle = simple_out, format=seq_format)
			simple_rev_it += 1
			continue


		# print reads with forward primers in opposite orientations to separate stream
		# these are likely off-target
		if (num_fwd_primers > 0) and (num_fwd_primers_in_rc > 0):
			SeqIO.write(record, handle = double_fwd_out, format=seq_format)
			double_it += 1
			continue

		# print concatemeric reads to separate stream
		# determine direciton of concatemer
		if (num_fwd_primers > 1) and (num_fwd_primers_in_rc == 0):
			SeqIO.write(record, handle = concatemeric_out, format=seq_format)
			concatemeric_fwd_it += 1
			continue
		
		if (num_fwd_primers == 0) and (num_fwd_primers_in_rc > 1):
			rc_record = record.reverse_complement(id=record.id + "(rc)", 
												  name = None, description = "")
			SeqIO.write(rc_record, handle = concatemeric_out, format=seq_format)
			concatemeric_rev_it += 1
			continue

		# finally, record reads with no fwd primer handles
		if (num_fwd_primers == 0) and (num_fwd_primers_in_rc == 0):

			SeqIO.write(record, handle = no_amp_handle_out, format=seq_format)
			no_handle_it += 1
			continue

		sys.stderr.write("Error: oops! this should never happen. Check decision tree above.")
		sys.stderr.write("\tforward primers in the fwd and rev orientations: (%d,%d)\n" % (num_fwd_primers, num_fwd_primers_in_rc))
simple_out.close()
concatemeric_out.close()
double_fwd_out.close()
no_amp_handle_out.close()

sys.stderr.write("############################# \n")

sys.stderr.write("# Finished processing!\n")
sys.stderr.write("# Total reads processed: %d \n" % it)

sys.stderr.write("# Reads with single fwd amplification handle (in fwd or rev direction): (%d,%d) [%.2f] \n"
					 % (simple_fwd_it, simple_rev_it, float((simple_fwd_it+simple_rev_it)/it)))


sys.stderr.write("# Concatemeric reads (in fwd or rev direction): (%d,%d) [%.4f] \n"
					 % (concatemeric_fwd_it, concatemeric_rev_it, float((concatemeric_fwd_it+concatemeric_rev_it)/it)))

sys.stderr.write("# Reads with forward amplification handle in both directions: %d [%.2f] \n"
					 % (double_it, float((double_it)/it)))

sys.stderr.write("# Reads with no fwd amplification handle : %d [%.2f] \n"
					 % (no_handle_it, float((no_handle_it)/it)))

sys.stderr.write("############################# \n")
