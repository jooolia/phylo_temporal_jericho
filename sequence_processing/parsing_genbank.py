#!/usr/bin/python
## Using python scripts adapted from Adina Chuang Howe (@teeniedeenie) 
## to download metadeta from Genbank. 
## Modified by: Julia Gustavsen (j.gustavsen@gmail.com)
## Date created: 11 September 2014
## Purpose: Divide retrieved genbank files up into a file that is easily readable by 
## making into 

import sys 

file=sys.argv[1]

genbank_record=open(file, 'r')

n = 0
l = []

for line in genbank_record:
   # print line
    if line.startswith("SOURCE"):
        newline=line.strip()
        source_divided=(newline.split(None, 1))
        print source_divided[1]     
    if "ORGANISM" in line:
        newline=line.strip()
        organism_divided=(newline.split(None, 1))
        print organism_divided[1]
    if "country" in line:
        newline=line.strip()
        country_divided=(newline.split("=", 1))
        print country_divided[1]
    if "isolate" in line:
        newline=line.strip()
        isolate_divided=(newline.split("=", 1))
        print isolate_divided[1]
    if "isolation_source" in line:
        newline=line.strip()
        isolation_source_divided=(newline.split("=", 1))
        print isolation_source_divided[1]
    if "collection_date" in line:
        newline=line.strip()
        collection_date_divided=(newline.split("=", 1))
        print collection_date_divided[1]
    # This is to deal with records that do not have these specific fields. 
    # I should make these into functions so that it is more efficient
    # The string is divided into 2 because the 2nd item is called when I go to print out the file. 
    try:
        country_divided 
    except NameError:
        country_divided = "" , ""
    try:
        collection_date_divided
    except NameError:
        collection_date_divided = "" , ""  
    try:
        isolate_divided 
    except NameError:
        isolate_divided = "" , ""

import csv


Header = "Filename","Source", "Organism","Isolate","Isolation_source", "Location","Date"
Genbank_data = file,source_divided[1], organism_divided[1], country_divided[1], isolate_divided[1], isolation_source_divided[1],collection_date_divided[1] 

import time
timestr = time.strftime("%Y_%m_%d")

genbank_out="Genbank_data_"+timestr+".csv"

def addToFile(file, what):
    outfile = open(genbank_out, "a" )
    writer = csv.writer(outfile)
    writer.writerow(what)
    
try:
   with open(genbank_out):
       addToFile(genbank_out, Genbank_data) 
except IOError:
   print 'Oh dear.'
   outfile = open(genbank_out, "w" )
   writer = csv.writer(outfile)
   writer.writerow(Header)
   writer.writerow(Genbank_data)
   outfile.close()
