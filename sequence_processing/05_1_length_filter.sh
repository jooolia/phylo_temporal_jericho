
## Author: Julia Gustavsen
## Trying to filter sequences before OTU picking so that I can get better OTUs. I am using Biopython to look for the motif and then print out the good sequences. 
## This will only be used for the viral sequences. 
## To-do: make note of where these motifs come from in the reference alignment
set -o nounset -o pipefail -o errexit

## not sure that I need this now that I am filtering at the nucleotide level. But using this as a space wehre I filter by length

#Use the python script for a number of files 
Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated

#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32
usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32



Sequences_gp23=$Translated_seqs_dir/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 
Sequences_MPL=$Translated_seqs_dir/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 



for files in $Sequences_gp23 $Sequences_MPL
do
$usearch8 -sortbylength $files -fastaout temp.fasta -minseqlength 66
mv temp.fasta $files
done