# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 10:41:33 2014

@author: julia
"""


### want to look at concatenating the R1 and R2 for 16s so that I can align and then cluster them....
## run file oas: python concatenate_R1_and_R2_for_discarded_non_merged.py discarded_R1.fasta discarded_R2.fasta concatenated_outfile.fasta
 
import sys
file1=sys.argv[1]
file2=sys.argv[2]
file3=sys.argv[3]

from Bio import SeqIO
records_f = SeqIO.parse(open(file1,"rU"), "fasta")
records_r = SeqIO.parse(open(file2,"rU"), "fasta")

## just use these 2 dicts instead of printing them...make a new dict that is updated.

forward_dict={}
for keyF in records_f:
    forward_id = keyF.id[:-12]
    keyF.description=""
    forward_dict[forward_id] = keyF.seq

reverse_dict={}
for keyR in records_r:
    reverse_id = keyR.id[:-12]
    keyR.description=""
    reverse_dict[reverse_id] = keyR.seq
 
forward_set_keys = set(forward_dict.keys())
reverse_set_keys = set(reverse_dict.keys())

## find ids that are in both forward and reverse
matching_keys = forward_set_keys.intersection(reverse_set_keys)

from Bio.SeqIO.FastaIO import FastaWriter

# normal biopython writer wraps every 60 characterrs, but mothur and qiime are not happy with that. 
## overwrite old file
concat_temp=open(file3, 'w')
concat_temp.close() 

## progress bar
#bar = progressbar.ProgressBar(maxval=len(matching_keys), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
#i=0

length_matching_keys=len(matching_keys)
count = 0
import difflib


concat_temp=open(file3, 'w')
concat_temp.close() 
concat_temp=open(file3, 'a')
writer = FastaWriter(concat_temp, wrap=None)    
writer.write_header()

for key in matching_keys:
    count = count + 1
    print count, "/", length_matching_keys

    seq=difflib.SequenceMatcher(a=forward_dict[key].lower(), b=reverse_dict[key].reverse_complement().lower())
    print seq.ratio()
    if seq.ratio() < 0.40:
    # add the reverse complement of from the reverse:
        concatenated_R1_with_R2 = SeqIO.SeqRecord(seq= forward_dict[key] + reverse_dict[key].reverse_complement(), id = key, description = "")
        writer.write_record(concatenated_R1_with_R2)


writer.write_footer()            
concat_temp.close() 

