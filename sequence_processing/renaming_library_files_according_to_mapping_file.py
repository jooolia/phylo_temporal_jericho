#!/usr/bin/env python

import csv
import os

mapping_csv= '../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv'

Rename_dict={}

with open(mapping_csv) as mapping_file:
    csv_mapper = csv.reader(mapping_file)
    for row in csv_mapper:
#        print row[6]
        Rename_dict[row[0]]=row[6]

suffix = '.fastq'
Read1 = '_R1'
Read2 = '_R2'

 # List the files in the current directory
os.chdir('../data')

import re

for filename in os.listdir('.'):
     root, extension = os.path.splitext(filename)
     print root
    # some names have a trailing 001 that needs to be removed
     simple_root = re.sub(r"\_001$", "", root)
     print simple_root
     print extension
     if Read1 in simple_root: 
          print simple_root
          stripped_root = simple_root[:-3] #removes the R1 
          if stripped_root in Rename_dict:
            print "yay"
            os.rename(filename, 'Lib_pool'+ Rename_dict[stripped_root]+ Read1 + suffix)
            print "Successfully renamed", filename, "to", 'Lib_pool'+ Rename_dict[stripped_root]+ Read1 + suffix
        
     if Read2 in simple_root: 
          print simple_root
          stripped_root = simple_root[:-3] #removes the R2
          if stripped_root in Rename_dict:
              os.rename(filename, 'Lib_pool'+ Rename_dict[stripped_root]+ Read2 + suffix)
              print "Successfully renamed.", filename, "to", 'Lib_pool'+ Rename_dict[stripped_root]+ Read2 + suffix

os.chdir('../scripts')
