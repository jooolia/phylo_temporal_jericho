#!/bin/bash
## Author: Julia Gustavsen
## Modified by: Julia Gustavsen

#source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit

Sequence_directory=../results/Divided_by_primers/All_libs_combined

Total_18s=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.ng.fasta 
Total_16s=$Sequence_directory/Total_16s_concat_annotated.good.good.filter.ng.fasta

Total_gp23=$Sequence_directory/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng.fasta

MPL_reads=$Sequence_directory/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng.fasta

seqtk=/home/labop/Data/software/seqtk/seqtk
#seqtk=/home/julia/seqtk/seqtk

### Databases
gp23_db=../data/gp23_sequence_references/gp23_nuc
MPL_db=../data/MPL_sequence_references/MPL_nuc


echo =========================================================================
echo "remove contaminants"
echo =========================================================================


echo "blast against the gene"

#$Total_18s $Total_16s 

for files in $MPL_reads $Total_gp23 
do 

basename $files
short_name=`basename -s .fasta $files`
if [[ $( echo $files) == *MPL* ]]
 then 
 database=$MPL_db
elif [[ $( echo $files) == *gp23* ]] 
then
database=$gp23_db
else echo "weird file!"
  fi


## skip temporarily
tblastx -query $files -db $database -evalue 1e-3 -outfmt 6 -out $Sequence_directory/$short_name"_blast.result" -num_threads 10 -max_target_seqs 1

echo "Seqtk for pulling out the good hits"
## was too slow! Just want to take the ids that are in the blast.result! More efficient for time and space. 

#python parse_blast_results_for_good_and_bad.py $Sequence_directory/$short_name"_blast.xml" $files  $Sequence_directory/$short_name"_good_by_blast.fasta" $Sequence_directory/$short_name"_bad_by_blast.fasta" 

## weird that the blast has given more than 1 target. :(
echo "cutting the sequences up"	
cut -f 1 $Sequence_directory/$short_name"_blast.result" > keep.list
#Extract sequences with names in file name.lst, one sequence name per line:

echo "taking the unique hits"
uniq keep.list > keep_unique.list

echo "now using subseq in the unique hits"

$seqtk subseq $files keep_unique.list > $Sequence_directory/$short_name"_good_by_blast.fasta"


done

#deactivate