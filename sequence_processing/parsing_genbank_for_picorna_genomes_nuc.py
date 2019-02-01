#!/usr/bin/python
## Using python scripts adapted from Adina Chuang Howe (@teeniedeenie) 
## to download metadeta from Genbank. 
## Modified by: Julia Gustavsen (j.gustavsen@gmail.com)
## Date created: 11 September 2014
## Purpose: Divide retrieved genbank files up into a file that is easily readable by 
## making into 


import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#genome=SeqIO.read(sys.argv[1], 'genbank')


import re

def findWholeWord(w):
    return re.compile(r'\b({0})\b'.format(w), flags=re.IGNORECASE).search


n = 0
l = []
gp20_list = []
target='RNA-dependent RNA polymerase'
alt_target='RdRp'
for record in list(SeqIO.parse(sys.argv[1], 'genbank')):
    #print record
    print len(record)
    org = record.annotations["source"]
    gi_number = record.annotations["gi"]
    for feat in record.features:
        if feat.type == "CDS":
            #print feat.qualifiers
            try:
                #if findWholeWord('gp20')(feat.qualifiers['note'][0]) or 'gp20' in feat.qualifiers['product'][0]:
                if target in feat.qualifiers['note'][0] or target in feat.qualifiers['product'][0] or target in feat.qualifiers['gene'][0] or alt_target in feat.qualifiers['gene'][0]:
                    print "yes"
                    print feat
                    gp20_translation = feat.qualifiers['translation'][0]
                    gp20_nucleotide = record.seq
		            ## want to get the specific nucleotide location
                    #print gp20_nucleotide
                    DNA_pol_nuc=gp20_nucleotide[feat.location.start.position:feat.location.end.position]
                    database_id = feat.qualifiers['db_xref'][0]
                    note = feat.qualifiers['note'][0]
                    product = feat.qualifiers['product'][0]
                    info_about_organism_and_gene = org, gi_number,  product, DNA_pol_nuc
                    gp20_list.append(info_about_organism_and_gene)
            except KeyError:
                print('This feature does not have a note')
                try:
                    if target in feat.qualifiers['product'][0] or target in feat.qualifiers['gene'][0] or alt_target in feat.qualifiers['gene'][0]:
                        print "yes"
                        print feat
                        gp20_translation = feat.qualifiers['translation'][0]
                        gp20_nucleotide = record.seq
                        DNA_pol_nuc=gp20_nucleotide[feat.location.start.position:feat.location.end.position]
                        #print gp20_nucleotide
                        #database_id = feat.qualifiers['db_xref'][0]
                        #note = "no note"
                        product = feat.qualifiers['product'][0]
                        info_about_organism_and_gene = org, gi_number, product, DNA_pol_nuc
                        gp20_list.append(info_about_organism_and_gene)
                except KeyError:
                    print('This feature does not have a product')
                    pass


print gp20_list
import csv


Header = "Organism", "gi_number", "feature_product","RdRp_nucleotide"
# Genbank_data = file,source_divided[1], organism_divided[1], country_divided[1], isolate_divided[1], isolation_source_divided[1],collection_date_divided[1] 

import time
timestr = time.strftime("%Y_%m_%d")

genbank_out=sys.argv[2]

outfile = open(genbank_out, "w" )
writer = csv.writer(outfile, delimiter='\t',)
writer.writerow(Header)
writer.writerows(gp20_list)
outfile.close()

