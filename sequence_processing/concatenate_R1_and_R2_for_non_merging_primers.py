
### want to look at concatenating the R1 and R2 for 16s so that I can align and then cluster them....


import sys
file1=sys.argv[1]
file2=sys.argv[2]
file3=sys.argv[3]


def file_len_fasta(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return (i + 1)/2

import itertools
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

records_r = SeqIO.parse(open(file2,"rU"), "fasta")
records_f = SeqIO.parse(open(file1,"rU"), "fasta")
destination = open(file3,'w')
destination.close() 

# normal biopython writer wraps every 60 characterrs, but mothur and qiime are not happy with that. 

## since it can be such a big file I can try to slowly write it out, instead of keeping it all in memory. 


length_file_1 = file_len_fasta(file1)

count = 0

import difflib

## first make sure the file only contains stuff from this go around. 
destination = open(file3,'w')
destination.close() 
destination = open(file3,'a')
writer = FastaWriter(destination, wrap=None)    
writer.write_header()



for (forward,reverse)in itertools.izip(records_f, records_r):
    count = count + 1
    print count, "/", length_file_1
    seq=difflib.SequenceMatcher(a=forward.seq.lower(), b=reverse.reverse_complement().lower())
    print seq.ratio()
    ## this is because of the errors I saw in the illumina sequences where the r1 and r2 when reversed were very similar, although not identical. 
    if seq.ratio() < 0.40:
        concatenated_R1_with_R2 = forward.seq + reverse.reverse_complement()
       # print concatenated_R1_with_R2.seq
        concatenated_R1_with_R2.id= forward.id
        concatenated_R1_with_R2.description = ""
      #  print concatenated_R1_with_R2
        writer.write_record(concatenated_R1_with_R2)

  
writer.write_footer()            
destination.close() 
