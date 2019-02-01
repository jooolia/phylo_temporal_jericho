#!/usr/bin/python

## Using python scripts adapted from Adina Chuang Howe (@teeniedeenie) 
## to download metadeta from Genbank. 
## Modified by: Julia Gustavsen (j.gustavsen@gmail.com)
## Purpose: Fetch genbank data for accession numbers from relevant papers. 
import urllib2
import os
import sys
import time

if len(sys.argv) != 3:
    print "USAGE: fetch_genome.py <genome_id_list> <out_dir>"
    sys.exit(1)

url_template = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=%s&rettype=fasta&retmode=text"
# can change rettype to fasta to get just the fasta files &rettype=fasta&retmode=text

if not os.path.exists(sys.argv[2]):
    os.makedirs(sys.argv[2])

for id in open(sys.argv[1]):
    id = id.strip()
    if id == "":
        continue

    sys.stdout.write("Fetching %s..." % id)
    sys.stdout.flush()
    gbk_out_file = os.path.join(sys.argv[2], id + ".fasta")
    if os.path.exists(gbk_out_file):
        print "already fetched"

    open(gbk_out_file, "w").write(urllib2.urlopen(url_template % id).read())
    print "Done"
    time.sleep(1.0/3)

 