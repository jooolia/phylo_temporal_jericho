## Count primers found in fastq files quickly with output as csv file
##Author: Julia Gustavsen
##Date last updated: 6 March 2014
## usage: python count_merged_primer_files.py file_R1.fastq file_R2.fastq merged_file.fastq
#python split_adapter_trimmed_libraries_by_primer.py $trimmed_lib $lib_X_primerX.fastq
#need to use the R1 and R2 here

import csv
import os
import sys


#read in fastq files
#fastq1 ='../results/Divided_by_primers/Lib_pool11_16s_R1.fastq'
#fastq2 ='../results/Divided_by_primers/Lib_pool11_16s_R2.fastq'
#fastq_merged='../results/Divided_by_primers/Lib_pool11_16s_merged_pear.fastq.assembled.fastq'

fastq1 = sys.argv[1]
fastq2 = sys.argv[2]
fastq_merged = sys.argv[3]


drive, path = os.path.splitdrive(fastq1)
path, filename = os.path.split(path)
file_part, extensions = os.path. splitext(filename)
#print('Drive is %s Path is %s and file is %s' % (drive, path, filename))
#print file_part
chop_up_filename=file_part.split("_")

library_identifier = chop_up_filename[1]
primer_name = chop_up_filename[2]


def countfastqrecords(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines/4


def percent_merged(fastq1, fastq2, fastq_merged):
    count_reads1 = countfastqrecords(fastq1)
    count_reads2 = countfastqrecords(fastq2)
    count_merged = countfastqrecords(fastq_merged)
    percent_merged_via_reads1 = (count_merged/float(count_reads1))*100
    percent_merged_via_reads2 = (count_merged/float(count_reads2))*100
    #these two equations should give the same percentage
    assert percent_merged_via_reads1 == percent_merged_via_reads2
    return percent_merged_via_reads1

number_of_fastq_records_read1=countfastqrecords(fastq1)
number_of_fastq_records_read2=countfastqrecords(fastq2)
number_of_merged_records=countfastqrecords(fastq_merged)
percent_merged_fastq_file=percent_merged(fastq1,fastq2, fastq_merged)
print number_of_fastq_records_read1
print percent_merged_fastq_file


Header = "Filename","Library_pool", "Primer_name","Number_of_Fastq_records_in_read1","Number_of_Fastq_records_in_read2", "Number_merged","Percent_merged"
Merged_count_data = fastq1,library_identifier, primer_name, number_of_fastq_records_read1, number_of_fastq_records_read2, number_of_merged_records, percent_merged_fastq_file


import time
timestr = time.strftime("%Y_%m_%d")

merge_out="../results/Merge_counts/merge_counts.csv"

def addToFile(file, what):
    outfile = open(merge_out, "a" )
    writer = csv.writer(outfile)
    writer.writerow(what)
    
try:
   with open(merge_out):
       addToFile(merge_out, Merged_count_data) 
except IOError:
   print 'Oh dear.'
   outfile = open(merge_out, "w" )
   writer = csv.writer(outfile)
   writer.writerow(Header)
   writer.writerow(Merged_count_data)
   outfile.close()