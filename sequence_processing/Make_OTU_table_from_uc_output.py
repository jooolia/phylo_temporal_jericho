
#Purpose: Python script to sort the OTU table so that it is only centroids and the number of reads in each centroid. Will then feed this to R. 
#Author: Julia Gustavsen
#Created: 19 August 2013
#Updated: 27 April 2014 
#usage; python Make_OTU_table_from_uc_output.py uclust_by_primer.uc Output_OTU_table.csv

import sys
#sys.arg[0] is script name


Uclust_uc = sys.argv[1]
#Uclust_uc ="../results/Divided_by_primers/All_libs_combined_OTUs_picked/Total_16s_R1_otus_97.00_identity.uc"

OTU_table = sys.argv[2] 

Usearch_global_uc = sys.argv[3]
#OTU_table = "../results/Divided_by_primers/All_libs_combined_OTUs_picked/test_OTU_table.csv"

#Parse uclust output from the Miseq sequences into something that can be parsed into an OTU table. I have renamed all the reads to start with the pool name. 
uclust = open(Uclust_uc, 'r')
OTUDictionary = {}
for line in uclust:
    fields = line.split()
    otu_id = fields[1]
    sections = fields[8].split("_") #split up the description so that I can get out the library pool name from file name
    name_and_size = fields[8].split(";") #split description so I can get at the number of occurrences
    OTU_name = name_and_size[0]
    OTU_abun = name_and_size[1].split("=") #get rid of the "size" part of this.
    pool_name = sections[0]
    centroid_id = [OTU_name, OTU_abun[1]]
    if fields[0] == 'S': # add query to dictionary as new sequence
        OTUDictionary[otu_id] = [centroid_id] 
    elif fields[0] == 'H': # add hit to match list
        OTUDictionary[otu_id].append(centroid_id)

#Test to see if I can just add these extra reads to this dictionary easily
#but have to see if it matches what is there already. 
uclust = open(Usearch_global_uc, 'r')
for line in uclust:
    fields = line.split()
    otu_id = fields[1]
    seed_id = fields[9]
    sections = fields[8].split("_") #split up the description so that I can get out the library pool name from file name
    name_and_size = fields[8].split(";") #split description so I can get at the number of occurrences
    OTU_name = name_and_size[0]
    OTU_abun = name_and_size[1].split("=") #get rid of the "size" part of this.
    pool_name = sections[0]
    centroid_id = [OTU_name, OTU_abun[1]]
    if seed_id in OTUDictionary.values():

        
    if fields[0] == 'S': # add query to dictionary as new sequence
        OTUDictionary[otu_id] = [centroid_id] 
    elif fields[0] == 'H': # add hit to match list
        OTUDictionary[otu_id].append(centroid_id)


        


OTUTable = {}

for key, value in OTUDictionary.items():
    otu_id = key
    for matching_reads_in_OTU in value:
        read_name = matching_reads_in_OTU[0]
        split_name = matching_reads_in_OTU[0].split("_")  #split up the description so that I can get out the library pool name from file name
        pool = split_name[0]
        OTU_abun = int(matching_reads_in_OTU[1]) #change this to an integer so that I can add the "abundances" together
        try:
            OTUTable[(otu_id),(pool)] += OTU_abun
        except:
            OTUTable[(otu_id),(pool)] = OTU_abun
    
print OTUTable

# Print the OTU table to a CSV file. 
    
import csv

def getUniqueValues(seq):
    "Return sorted list of unique values in sequence"
    values = list(set(seq))
    values.sort()
    return values

def dataArray(data2d, rowIterField=0, rowLabel='', defaultVal=''):
    # get all unique unit and test labels
    rowLabels = getUniqueValues(key[rowIterField] for key in data2d)
    colLabels = getUniqueValues(key[1-rowIterField] for key in data2d)

    # create key-tuple maker
    if rowIterField==0:
        key = lambda row,col: (row, col)
    else:
        key = lambda row,col: (col, row)

    # header row
    yield [rowLabel] + colLabels
    for row in rowLabels:
        # data rows
        yield [row] + [data2d.get(key(row,col), defaultVal) for col in colLabels]

def main():
    with open(OTU_table, 'wb') as outf:
        outcsv = csv.writer(outf)
        outcsv.writerows(dataArray(OTUTable, 0, 'OTU_id'))

if __name__=="__main__":
    main()

