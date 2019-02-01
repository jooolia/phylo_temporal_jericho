
## I might not necessarily need to align files, but just make sure that some motifs are there that I had found in the alignment.

## Would be good to print out how many files I get rid of. 

## run file as:
#python filter_sequences_by_motifs.py amino_acid_input.faa motifs_regex filtered_output.faa non_matching.faa

import sys
file1=sys.argv[1]
file2=sys.argv[2]


from Bio import SeqIO

records_f = SeqIO.parse(open(file1,"rU"), "fasta")


import re
records_no_stop = []
for record in records_f:
    if "X" in str(record.seq):
	    print "yes"
    else:
         print "clean seq"
         records_no_stop.append(record.id)

destination = open(file2,'w')

for item in records_no_stop:
  destination.write("%s\n" % item)

destination.close()