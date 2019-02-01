#
#Author: Julia Gustavsen
#Purpose: Count the length of all the reads on the files. 
#Usage: python length_of_fastQ_reads.py fastq_file.fastq output_file.csv
#To do:Want to get all the files together in a format so that I can use R with melt to do histograms, fastest way would be to run as bash script and then paste together in R. 

import sys

fastq_file = file( sys.argv[1], "r" )

print(sys.argv[1])

import re

#Give file an appropriate name
lengths_of_reads = open( sys.argv[2], 'w') 

#write header
lengths_of_reads.write(sys.argv[1]+ "\n")

for line in fastq_file:
    clean_read = line.strip()
    if clean_read.startswith("@"):
        sequence_read = (fastq_file.next()).strip()
 #       print sequence_read, "1st sequence read."
    #    print len(sequence_read) , "length of sequence read"
	sequence_length= str(len(sequence_read))
	lengths_of_reads.write(sequence_length + "\n")


fastq_file.close()
lengths_of_reads.close()
