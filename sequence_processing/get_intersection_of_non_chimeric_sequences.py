
# usage: python get_union_of_non_chimeric_sequences.py combined_non_chimeras.fasta combined_unique_non_chimeras.fasta


import sys
denovo_non_chimeras=sys.argv[1]
ref_non_chimeras=sys.argv[2]

## non-chimeras present in both methods
combined_non_chimeras=sys.argv[3]

import itertools
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
records_denovo= SeqIO.parse(open(denovo_non_chimeras,"rU"), "fasta")
records_ref= SeqIO.parse(open(ref_non_chimeras,"rU"), "fasta")



denovo_simple={}
for record in records_denovo:
    header=record.id
    sequence=record.seq
    denovo_simple[header] = record.seq
    
  
ref_simple={}
for record in records_ref:
    header=record.id
    sequence=record.seq
    ref_simple[header] = record.seq  
   
   
records_denovo_set=set(denovo_simple)
#print records_denovo_set

records_ref_set=set(ref_simple)

print len(records_denovo_set)

#print records_ref_set
print len(records_ref_set)


## find records that are non-chimeric in both sets
intersection_records = records_ref_set.intersection(records_denovo_set)

print len(intersection_records)

# ok use these keys to pull out the records from one of the sets of 

length_intersection=len(intersection_records)
count = 0

destination = open(combined_non_chimeras,'w')
writer = FastaWriter(destination, wrap=None) 
writer.write_header()

for record in intersection_records:
    count = count + 1
    print count, "/", length_intersection
    header=record
    sequence=denovo_simple[record]
    new = SeqIO.SeqRecord(seq = sequence,id=header, description="")      
    writer.write_record(new)

writer.write_footer()  
destination.close()

       
