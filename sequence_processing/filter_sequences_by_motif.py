
## I might not necessarily need to align files, but just make sure that some motifs are there that I had found in the alignment.

## Would be good to print out how many files I get rid of. 

## run file as:
#python filter_sequences_by_motifs.py amino_acid_input.faa motifs_regex filtered_output.faa non_matching.faa

import sys
file1=sys.argv[1]
file2=sys.argv[3]
file3=sys.argv[4]
motifs=sys.argv[2]


from Bio import SeqIO

records_f = SeqIO.parse(open(file1,"rU"), "fasta")

import re
non_matching_motif_records = {}
matching_motif_records = {}
for record in records_f:
    if re.search(motifs, str(record.seq)):
	    print "yes"
	    matching_motif_records[record.id] = record.seq
    else:
        print "no match"
        non_matching_motif_records[record.id] = record.seq



from Bio.SeqIO.FastaIO import FastaWriter
destination = open(file2,'w')
destination.close() 

length_file_1 = len(matching_motif_records)

destination = open(file2,'w')
destination.close() 
count = 0
for key, value in matching_motif_records.iteritems():
    destination = open(file2,'a')
    count = count + 1
    print count, "/", length_file_1, "good records"
    writer = FastaWriter(destination, wrap=None)    
    writer.write_header()
    ## transform into biopython record
    record = SeqIO.SeqRecord(seq= value, id = key, description = "")
    writer.write_record(record)
    writer.write_footer()            
    destination.close() 

length_file_2 = len(non_matching_motif_records)
length_original = length_file_1 + length_file_2

destination = open(file3,'w')
destination.close() 
count = 0
for key, value in non_matching_motif_records.iteritems():
    destination = open(file3,'a')
    count = count + 1
    print count, "/", length_file_2, "garbage records"
    writer = FastaWriter(destination, wrap=None)    
    writer.write_header()
    ## transform into biopython record
    record = SeqIO.SeqRecord(seq= value, id = key, description = "")
    writer.write_record(record)
    writer.write_footer()            
    destination.close() 



print "Number of original fasta records", length_original
print "Number of good records",  length_file_1 
print "Number of garbase records", length_file_2  
