## Prints out CSV with long format of length of sequences in the file. Could slim down the file name later. 

from Bio import SeqIO
import sys

fastq_file=sys.argv[1]
length_dir=sys.argv[2]
sizes = [len(rec) for rec in SeqIO.parse(fastq_file, "fastq")]
print len(sizes), min(sizes), max(sizes)
#sizes


# print to file
import csv

Header = "Filename","Length"


import time
timestr = time.strftime("%Y_%m_%d")

length_out=length_dir+ "/Length_data.csv"

def addToFile(file, what):
    outfile = open(length_out, "a" )
    writer = csv.writer(outfile)
    writer.writerow(what)

# actually print it out in long format for easier work in R

print range(len(sizes))

#print sizes[2]
for one_record in range(len(sizes)):
	#print one_record
	Fastq_length_data = fastq_file,sizes[one_record]
	print Fastq_length_data
	try:
            with open(length_out):
                addToFile(length_out, Fastq_length_data) 
        except IOError:
            print 'Oh dear.'
            outfile = open(length_out, "w" )
            writer = csv.writer(outfile)
            writer.writerow(Header)
            writer.writerow(Fastq_length_data)
            outfile.close()
	