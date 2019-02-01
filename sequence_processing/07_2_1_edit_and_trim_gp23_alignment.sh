

set -o nounset -o pipefail -o errexit


Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated

#trimal=/home/julia/trimal/source/trimal
trimal=/home/labop/trimal/source/trimal

#clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64
clustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64

#hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 
hmm_build=hmmbuild

#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32
usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32


gp23_alignment=../results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.fasta
Sequences_gp23=$Translated_seqs_dir/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 


###### First profile 
head -n 50000 $Sequences_gp23"_clustal_aligned_all_profile.fasta" > $Sequences_gp23"_clustal_aligned_all_profile_first50000.fasta"
tail -n 50000 $Sequences_gp23"_clustal_aligned_all_profile.fasta" > $Sequences_gp23"_clustal_aligned_all_profile_last50000.fasta"

$trimal -keepheader -in $Sequences_gp23"_clustal_aligned_all_profile.fasta" -out $Sequences_gp23"_clustal_aligned_all_profile_trimmed.fasta" -selectcols { 1-53,656-664 } -fasta

#aliview ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profile_first50000.fasta


### Frofile B
head -n 50000 $Sequences_gp23"_clustal_aligned_all_profileB.fasta" > $Sequences_gp23"_clustal_aligned_all_profileB_first50000.fasta"
tail -n 50000 $Sequences_gp23"_clustal_aligned_all_profileB.fasta" > $Sequences_gp23"_clustal_aligned_all_profileB_last50000.fasta"

#aliview ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profileB_first50000.fasta

$trimal -keepheader -in $Sequences_gp23"_clustal_aligned_all_profileB.fasta" -out $Sequences_gp23"_clustal_aligned_all_profileB_trimmed.fasta" -selectcols { 1-25,465-466 } -fasta


### Profile C
head -n 50000 $Sequences_gp23"_clustal_aligned_all_profileC.fasta" > $Sequences_gp23"_clustal_aligned_all_profileC_first50000.fasta"
tail -n 50000 $Sequences_gp23"_clustal_aligned_all_profileC.fasta" > $Sequences_gp23"_clustal_aligned_all_profileC_last50000.fasta"

#aliview ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profileC_first50000.fasta

$trimal -keepheader -in $Sequences_gp23"_clustal_aligned_all_profileC.fasta" -out $Sequences_gp23"_clustal_aligned_all_profileC_trimmed.fasta" -selectcols { 1-20,432-495 } -fasta


### Profile D
head -n 50000 $Sequences_gp23"_clustal_aligned_all_profileD.fasta" > $Sequences_gp23"_clustal_aligned_all_profileD_first50000.fasta"
tail -n 50000 $Sequences_gp23"_clustal_aligned_all_profileD.fasta" > $Sequences_gp23"_clustal_aligned_all_profileD_last50000.fasta"

#aliview ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profileD_first50000.fasta

$trimal -keepheader -in $Sequences_gp23"_clustal_aligned_all_profileD.fasta" -out $Sequences_gp23"_clustal_aligned_all_profileD_trimmed.fasta" -selectcols { 1-5,344-344 } -fasta

### Profile E
head -n 50000 $Sequences_gp23"_clustal_aligned_all_profileE.fasta" > $Sequences_gp23"_clustal_aligned_all_profileE_first50000.fasta"
tail -n 50000 $Sequences_gp23"_clustal_aligned_all_profileE.fasta" > $Sequences_gp23"_clustal_aligned_all_profileE_last50000.fasta"

#aliview ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profileE_first50000.fasta


$trimal -keepheader -in $Sequences_gp23"_clustal_aligned_all_profileE.fasta" -out $Sequences_gp23"_clustal_aligned_all_profileE_trimmed.fasta" -selectcols { 1-40,516-581 } -fasta


cat $Sequences_gp23"_clustal_aligned_all_profile_trimmed.fasta" $Sequences_gp23"_clustal_aligned_all_profileB_trimmed.fasta" $Sequences_gp23"_clustal_aligned_all_profileC_trimmed.fasta" $Sequences_gp23"_clustal_aligned_all_profileD_trimmed.fasta" $Sequences_gp23"_clustal_aligned_all_profileE_trimmed.fasta" > $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed.fasta"


sed 's/-//g' $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed.fasta" > $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_sed.fasta"

### look at this ...
head -n 50000 $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_sed.fasta" > $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_first50000.fasta"
tail -n 50000 $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_sed.fasta" > $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_last50000.fasta"

## make onto 1 line. 
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_sed.fasta" > $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_sed_1L.fasta"


$usearch8 -sortbylength $Sequences_gp23"_clustal_aligned_all_profileALL_trimmed_sed_1L.fasta" -fastaout $Sequences_gp23"_clustal_aligned_all_profile_trimmed_length.fasta" -minseqlength 145

mv $Sequences_gp23"_clustal_aligned_all_profile_trimmed_length.fasta" $Translated_seqs_dir/Total_gp23_pear.assembled_annotated_quality_length_filt.fasta
## maybe have to trim all separately and then concat and then length trim...