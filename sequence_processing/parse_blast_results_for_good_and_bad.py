## file for parsing blast
## need to see those that are kept and those that are not kept. 

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys

blast_xml=sys.argv[1]
original_fasta=sys.argv[2]
good_records=sys.argv[3]
bad_records=sys.argv[4]

## so want to print the records that are found in the blast xml file to a new file
## first start by getting id of the records and then compare to the fasta file record
blast_records = NCBIXML.parse(open(blast_xml))

good_ids=[]
bad_ids=[]
for record in blast_records:
    #print record
    #print dir(record)
    print "record id", record.query
    if record.alignments:
        # to print the "best" matches bit-score
        print record.alignments[0].hsps[0].expect
        # to print the "best" matches bit-score
        print record.alignments[0].hsps[0].score
        for alignment in record.alignments:
            identities_list=[]
            for something in alignment.hsps:
                #print dir(something)
                identities_list.append(something.identities)
                print something.identities
            # if any of the identities in the high-scoring pair
            # is low I will get rid of that sequence    
            if any(i >= 20 for i in identities_list):
                print "yes"
                good_ids.append(record.query)
            else:
                 print "identity too low"
                 bad_ids.append(record.query)

                  ## but only want id added once....
    elif not record.alignments:
        print "no hit"
        bad_ids.append(record.query)
    else:
        print "something weird"

#print "List of good ids", good_ids
#print "List of bad ids", bad_ids


original_fasta_file = SeqIO.parse(open(original_fasta,"rU"), "fasta")

destination_good = open(good_records,'w')
destination_bad = open(bad_records,'w')


for reads in original_fasta_file:
    print reads.id
    if reads.id in good_ids:
        print "hoorah"
        SeqIO.write(reads, destination_good, "fasta")
    elif reads.id in bad_ids:
        print "we did not find it"
        SeqIO.write(reads, destination_bad, "fasta")

destination_good.close()
destination_bad.close()