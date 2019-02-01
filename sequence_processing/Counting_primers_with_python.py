## Count primers found in fastq files quickly with output as csv file
## usage: python Counting_primers_with_python.py file.fastq


import csv
import os
import sys

#primer_reg_ex_csv= '/run/media/julia/OS/Users/Julia.Julia-PC/Documents/School/Thesis/ThesisProjects/JerichoAndSOGsequencing/data/Primers_reg_ex_Julia_July2013.txt'

primer_reg_ex_csv = '../../McGill_SIP_Kuwait_JAG_check/Primers_JAG_with_gp20_reg_ex_Jan2014.txt'

primer_reg_ex_dict={}
with open(primer_reg_ex_csv) as mapping_file:
    primer_mapper = csv.reader(mapping_file, delimiter='\t')
    for row in primer_mapper:
       #print row
       primer_reg_ex_dict[row[0]]=row[1]
        
#primer_reg_ex_dict


#fastqfile='test2_for_regex.fastq'
fastqfile = sys.argv[1]
out_dir = sys.argv[2]

print(sys.argv[1])
    
with open(fastqfile, 'r') as content_file:
    content_fastq = content_file.read()

import re
  
def Count_Primer_in_Fastq(read_in_file, primer_dict):    
    list_counts={}
    for key, value in primer_dict.iteritems():
        match_primer=re.findall(value, read_in_file)
        count=len(match_primer)
       # print key, value, count
        list_counts[key]=count    
    return list_counts

def countfastqrecords(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines/4

Counted_primers_in_Fastq=Count_Primer_in_Fastq(content_fastq, primer_reg_ex_dict)

number_of_fastq_records = countfastqrecords(fastqfile)

Header = "Filename","Number of Fastq records", "16s_F", "16s_R","18s_F","18s_R", "gp23_F", "gp23_R", "MPL_F", "MPL_R", "CPS1", "CPS8"
Primer_count_data = fastqfile, number_of_fastq_records, Counted_primers_in_Fastq['16s_F'], Counted_primers_in_Fastq['16s_R'], Counted_primers_in_Fastq['18s_F'], Counted_primers_in_Fastq['18s_R'], Counted_primers_in_Fastq['gp23_F'],  Counted_primers_in_Fastq['gp23_R'], Counted_primers_in_Fastq['MPL_F'],  Counted_primers_in_Fastq['MPL_R'], Counted_primers_in_Fastq['CPS1'], Counted_primers_in_Fastq['CPS8']


import time
timestr = time.strftime("%Y_%m_%d")


test_out=out_dir + "primer_counts_" +timestr+ ".csv"

def addToFile(file, what):
    outfile = open( test_out, "a" )
    writer = csv.writer(outfile)
    writer.writerow(what)
    
try:
   with open(test_out):
       addToFile(test_out, Primer_count_data) 
except IOError:
   print 'Oh dear.'
   outfile = open( test_out, "w" )
   writer = csv.writer(outfile )
   writer.writerow(Header)
   writer.writerow(Primer_count_data)
   outfile.close()


