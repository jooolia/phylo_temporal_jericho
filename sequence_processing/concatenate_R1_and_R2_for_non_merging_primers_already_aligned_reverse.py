
### want to look at concatenating the R1 and R2 for 16s so that I can align and then cluster them....

import sys
file1=sys.argv[1]
file2=sys.argv[2]
file3=sys.argv[3]


def file_len_fasta(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return (i + 1)/2

import itertools
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

records_r = SeqIO.parse(open(file2,"rU"), "fasta")
records_f = SeqIO.parse(open(file1,"rU"), "fasta")


length_file_1 = file_len_fasta(file1)

count = 0

import difflib

## first make sure the file only contains stuff from this go around. 
destination = open(file3,'w')
destination.close() 
destination = open(file3,'a')
writer = FastaWriter(destination, wrap=None)    
writer.write_header()


regular_forward_dict={}
forward_list=[]
for record in records_f:
    #print record.id
    forward_id = record.id.split("_2_N")[0] # pulls out part that is the beginning and should be 
    #print forward_id
    forward_list.append(forward_id)
    regular_forward_dict[record.id] = record.seq


regular_reverse_dict={}
reverse_list=[]
for record in records_r:
    #print record.id
    reverse_id = record.id.split("_1_N")[0]
    #print reverse_id
    reverse_list.append(reverse_id)
    regular_reverse_dict[record.id]=record.seq

records_r.close()
records_f.close()

 
forward_set = set(forward_list)
reverse_set = set(reverse_list)


# ## find ids that are in both forward and reverse
matching_keys = forward_set.intersection(reverse_set)


#print matching_keys

matching_set =set(matching_keys)

def add_to_overall_dict_with_set(sequence_dictionary, matching_set, new_overall_dict):
     length_keys=len(sequence_dictionary.keys())
     count = 0
     for key, value in sequence_dictionary.iteritems():
         count = count + 1
         #print count, "/", length_keys
         for key_short in matching_set:
             if key.startswith(key_short):
                 print "yay"
               #  .setdefault workks with there is no key existing. 
                 new_overall_dict.setdefault(key_short,[]).extend([key, value])
     return new_overall_dict

overall_dict={}

matching_forward = add_to_overall_dict_with_set(regular_forward_dict, matching_set, overall_dict)

#print matching_forward

matching_forward_and_reverse = add_to_overall_dict_with_set(regular_reverse_dict, matching_set, matching_forward)

#print matching_forward_and_reverse
print "Now getting ready to write that big file!!! Concatenation here we come!"

 ## So now dictionary structure is short_key: [forwardid, forward_seq, reverseid, reverse_seq]
concat_temp=open(file3, 'w')
concat_temp.close() 
concat_temp=open(file3, 'a')
writer = FastaWriter(concat_temp, wrap=None)   
writer.write_header()
 ### maybe could make it faster by having the reverse and forward as a series of lists in a dictionary. Short_key: [keyF, keyR, seqF, seqR]
count=0
length_overall_keys=len(matching_forward_and_reverse.keys())

import difflib

for key, value in matching_forward_and_reverse.iteritems():
      count = count + 1
      print count, "/", length_overall_keys
      forward_key = value[0]
      print forward_key
      reverse_key = value[2]
      print reverse_key
      forward_seq = value[1]
      print forward_seq
      reverse_seq = value[3]
      print reverse_seq
      ## I have had a problem with seuqences in the R1 and R2 being almost identical when translated. I don't know how this would happen, 
      ## but it is seriously messing with my data and I start to lose a lot of sequences because I have to get rid of those 10 and greater
      ## and this is quite problematic. Thus I have used this section to compare R1 and R2 to each other and if they are less than 40% similar
      ## they are kept. I based the 40% on looking at some test cases and blasting the resulting sequences. 
      seq=difflib.SequenceMatcher(a=forward_seq.lower(), b=reverse_seq.lower())
      print seq.ratio()
      if seq.ratio() < 0.40:
           print "the plus side"
           R1_comes_first = SeqIO.SeqRecord(seq= forward_seq + reverse_seq, id = forward_key + "joined_with" + reverse_key, description = "")
           writer.write_record(R1_comes_first)           
      else:
        "sequences show a weird similarity"


writer.write_footer() 
concat_temp.close()

