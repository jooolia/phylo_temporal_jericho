
## I might not necessarily need to align files, but just make sure that some motifs are there that I had found in the alignment.

## Would be good to print out how many files I get rid of. 

## run file as:
#python filter_sequences_by_motifs.py amino_acid_input.faa motifs_regex filtered_output.faa non_matching.faa

import sys
file1=sys.argv[1]
file2=sys.argv[2]


from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

records_f = SeqIO.parse(open(file1,"rU"), "fasta")


from Bio import AlignIO
alignment = AlignIO.read(file1, "fasta")

print("Alignment length %i" % alignment.get_alignment_length())

# new_align = []
new_align = MultipleSeqAlignment([ ])

def check_gaps_less_7(record):
	no_term_gaps = record.seq.strip("-")
	#print no_term_gaps
	# print(len(no_term_gaps))
	# print(dir(no_term_gaps))
	no_gaps = no_term_gaps.ungap("-")
	# print no_gaps
	## instead need to look at the difference to see how many gaps....better than percent. 
	number_of_gaps = float(len(no_term_gaps))-float(len(no_gaps))
	print number_of_gaps 
	if number_of_gaps < 8.0:
		return True

## using a list comprehension to speed up! so much faster!
new_align = MultipleSeqAlignment([ record for record in alignment if check_gaps_less_7(record)])
		

AlignIO.write(new_align, file2, "fasta")
