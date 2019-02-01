## Count primers found in fastq files quickly with output as csv file
##Author: Julia Gustavsen
##Date last updated: 6 March 2014
## usage: python Counting_primers_with_python.py file_R1.fastq file_R2.fastq
#python split_adapter_trimmed_libraries_by_primer.py $trimmed_lib $lib_X_primerX.fastq
#need to use the R1 and R2 here

import csv
import os
import sys


#primer_reg_ex_csv= '/run/media/julia/OS/Users/Julia.Julia-PC/Documents/School/Thesis/ThesisProjects/JerichoAndSOGsequencing/data/Primers_reg_ex_Julia_July2013.txt'

primer_reg_ex_csv = '../../McGill_SIP_Kuwait_JAG_check/Primers_JAG_with_gp20_reg_ex_Jan2014.txt'


primer_reg_ex_dict={}
with open(primer_reg_ex_csv) as mapping_file:
    primer_mapper = csv.reader(mapping_file, delimiter='\t')
    for row in primer_mapper:
       print row
       primer_reg_ex_dict[row[0]]=row[1]
        

#read in fastq files
#fastq1 ='../data/Lib_pool98_R1_first200.fastq'
#fastq2 ='../data/Lib_pool98_R2_first200.fastq'
fastq1 = sys.argv[1]
fastq2 = sys.argv[2]

drive, path = os.path.splitdrive(fastq1)
path, filename = os.path.split(path)
file_part, extensions = os.path. splitext(filename)
#print('Drive is %s Path is %s and file is %s' % (drive, path, filename))

chop_up_filename=file_part.split("_")

Library_identifier = chop_up_filename[0]+ "_" +chop_up_filename[1]
read_orientation = chop_up_filename[2]

print Library_identifier,read_orientation


import itertools
from Bio import SeqIO
records_f = SeqIO.parse(open(fastq1,"rU"), "fastq")
records_r = SeqIO.parse(open(fastq2,"rU"), "fastq")


import re

def Find_Primer_in_Fastq(reads, primer_dict):    
    for key, value in primer_dict.iteritems():
        if re.search(value, str(reads.seq)): 
          # print key, reads.id, reads.seq
           return key, reads


for key, value in primer_reg_ex_dict.iteritems():
    filename_just_read_orientation = "../results/Divided_by_primers/%s_%s_R1.fastq" %(Library_identifier, key[:-2])
    destination = open(filename_just_read_orientation,'w')
    filename_just_read_orientation = "../results/Divided_by_primers/%s_%s_R2.fastq" %(Library_identifier, key[:-2])
    destination2 = open(filename_just_read_orientation,'w')
    destination.close()
    destination2.close()

 
def Write_sequence_matching_primer_and_mate_to_file(Find_Primer_in_Fastq, original_seq, orientation_match, mate_pair, orientation_mate):
    key, reads = Find_Primer_in_Fastq
    # filename_match = "../results/Divided_by_primers/%s_%s_%s.fastq" %(Library_identifier, key,orientation_match)
    # handle = open(filename_match, "a")
    # matching_record = SeqIO.write(reads, handle, "fastq")
    
    filename_just_read_orientation = "../results/Divided_by_primers/%s_%s_%s.fastq" %(Library_identifier, key[:-2],orientation_match)
    destination = open(filename_just_read_orientation,'a')
    matching_record = SeqIO.write(original_seq, destination, "fastq")
    
    # filename_mate = "../results/Divided_by_primers/%s_%s_%s.fastq" %(Library_identifier, key,orientation_mate)
    # handle2 = open(filename_mate, "a")
    # matching_record = SeqIO.write(mate_pair, handle2, "fastq")
    
    filename_just_read_orientation_mate = "../results/Divided_by_primers/%s_%s_%s.fastq" %(Library_identifier, key[:-2],orientation_mate)
    destination2 = open(filename_just_read_orientation_mate,'a')
    matching_record = SeqIO.write(mate_pair, destination2, "fastq")
    
    # handle.close()
    # handle2.close()
    destination.close()
    destination2.close()
    
#the only problem here is that if it matches in one direction it is not searched for in another direction. For now I can cat together all the records between primers. 
for (forward,reverse)in itertools.izip(records_f, records_r):
    #using the notation from the sequencer the headers are similar up to 12 characters from the end. So this is a way to compare the 2 records.
    if forward.id[:-12] == reverse.id[:-12]:
        rev_forward=forward.reverse_complement()
        rev_reverse=reverse.reverse_complement()
        if Find_Primer_in_Fastq(forward, primer_reg_ex_dict):
            print "Forward"
            Write_sequence_matching_primer_and_mate_to_file(Find_Primer_in_Fastq(forward, primer_reg_ex_dict),forward,"R1", reverse, "R2")
        elif Find_Primer_in_Fastq(rev_forward, primer_reg_ex_dict):
            Write_sequence_matching_primer_and_mate_to_file(Find_Primer_in_Fastq(rev_forward, primer_reg_ex_dict),forward,"R1", reverse, "R2")
            print "Forward reverse complement matching"
        elif Find_Primer_in_Fastq(reverse, primer_reg_ex_dict):    
            print "Reverse"
            Write_sequence_matching_primer_and_mate_to_file(Find_Primer_in_Fastq(reverse, primer_reg_ex_dict),reverse,"R2", forward, "R1")
        elif Find_Primer_in_Fastq(rev_reverse, primer_reg_ex_dict):
            print "Reverse reverse complement"
            Write_sequence_matching_primer_and_mate_to_file(Find_Primer_in_Fastq(rev_reverse, primer_reg_ex_dict),reverse,"R2", forward, "R1")
        else:
            print "no matches"
           
    else:
        print "sequences do not match up"
