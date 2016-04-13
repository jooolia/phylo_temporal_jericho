
import sys
original_file=sys.argv[1]
otu_table=sys.argv[2]
project_file=sys.argv[3]

import itertools
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

total_fasta = SeqIO.parse(open(original_file,"rU"), "fasta")
project_fasta = open(project_file,'w')
project_fasta.close() 
project_fasta = open(project_file,'a')

## read in the csv file and get header names
import csv
table_normalized_otus = open(otu_table, 'rb')
reader = csv.reader(table_normalized_otus, delimiter="\t")
headers = reader.next()
print headers


writer = FastaWriter(project_fasta, wrap=None)    
writer.write_header()

for records in total_fasta:
#    print records.name
    if records.name in headers:
        writer.write_record(records)

writer.write_footer()
