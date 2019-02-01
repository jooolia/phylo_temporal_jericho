### want to match the sequence id and then edit it

import os
import sys

fastq_in = sys.argv[1]
fastq_out = sys.argv[2]

drive, path = os.path.splitdrive(fastq_in)
path, filename = os.path.split(path)
file_part, extensions = os.path. splitext(filename)
print('Drive is %s Path is %s and file is %s' % (drive, path, filename))

chop_up_filename=file_part.split("_")

Library_identifier = chop_up_filename[0]+ "_" +chop_up_filename[1]
primer_type = chop_up_filename[2]

print Library_identifier,primer_type


import itertools
from Bio import SeqIO
records_f = SeqIO.parse(open(fastq_in,"rU"), "fastq")


destination = open(fastq_out,'w')

for reads in records_f:
  # print reads.id
   reads.id=Library_identifier+ "_" + reads.id
   reads.description=Library_identifier+ "_" + reads.description
 #  print reads.description
#   print reads.id, reads.seq
 #  print reads
   SeqIO.write(reads, destination, "fastq")


destination.close()
