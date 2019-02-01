## Python script to make a table that includes the sequence id for OTU centroids, the OTU number (as used in usearch) and the taxonomy that was assigned in motur

import sys

#file1="../results/Divided_by_primers/All_libs_combined_OTUs_picked/18s_OTUs.fasta"
#file2="../results/Divided_by_primers/All_libs_combined_OTUs_picked/18s_OTUs_numbered.fasta"
#file3="../results/Divided_by_primers/All_libs_combined_OTUs_picked/18s_OTUs.seed_v119.wang.taxonomy"
file1=sys.argv[1]
file2=sys.argv[2]
file3=sys.argv[3]

from Bio import SeqIO
original = SeqIO.parse(open(file1,"rU"), "fasta")
numbered = SeqIO.parse(open(file2,"rU"), "fasta")

import itertools

mapping_dictionary={}

for original, numbered in itertools.izip(original, numbered):
    # using Biopython eliminates a lot of issues where the fasta files have line breaks before another header line.
    # Makes it easy to just match up the files since they are in the same order. 
    # I checked that the first and last ones are matching
    mapping_dictionary[original.id]=numbered.id

# now to add the taxonomy assignment from Silva v119 database
with open(file3,'rb') as taxonomy_file:
    for line in taxonomy_file:
        newline=line.strip()
        fields = newline.split('\t')
        if fields[0] in mapping_dictionary:
            sequence_id=fields[0]
            taxonomy=fields[1]
            OTU_map= mapping_dictionary[sequence_id]
            new_fields=(OTU_map, taxonomy)
            # replace old field which had the OTU number only with a new list which has the OTU number and the taxonomy
            mapping_dictionary[sequence_id]= new_fields
            
            
    print mapping_dictionary
    

    
import csv
   
file4=sys.argv[4]   
   
writer = csv.writer(open(file4, 'wb'))
writer.writerow(["header_id", "otu_number", "silva_taxonomy"])
for key, value in mapping_dictionary.items():
   otu_number, taxonomy = value
   writer.writerow([key, otu_number, taxonomy])
   

