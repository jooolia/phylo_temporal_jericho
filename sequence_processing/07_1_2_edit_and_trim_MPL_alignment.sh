source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit


Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated

#trimal=/home/julia/trimal/source/trimal
trimal=/home/labop/trimal/source/trimal

#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32
usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32

Sequences_MPL=$Translated_seqs_dir/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 

head -n 50000 $Sequences_MPL"_clustal_aligned_all_profile.fasta" > $Sequences_MPL"_clustal_aligned_all_profile_first50000.fasta"
tail -n 50000 $Sequences_MPL"_clustal_aligned_all_profile.fasta" > $Sequences_MPL"_clustal_aligned_all_profile_last50000.fasta"

#aliview ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profile_first50000.fasta


$trimal -keepheader -in $Sequences_MPL"_clustal_aligned_all_profile.fasta" -out $Sequences_MPL"_clustal_aligned_all_profile_trimmed.fasta" -selectcols { 1-71,426-478 } -fasta

head -n 50000 $Sequences_MPL"_clustal_aligned_all_profile_trimmed.fasta" > $Sequences_MPL"_clustal_aligned_all_profile_trimmed_first50000.fasta"
tail -n 50000 $Sequences_MPL"_clustal_aligned_all_profile_trimmed.fasta" > $Sequences_MPL"_clustal_aligned_all_profile_trimmed_last50000.fasta"

## do not get rid of trailing - on the header. maybe this is causing problem... replace gaps except lines that start with 
sed -i '/^>/!s/-//g' $Sequences_MPL"_clustal_aligned_all_profile_trimmed.fasta"


$usearch8 -sortbylength $Sequences_MPL"_clustal_aligned_all_profile_trimmed.fasta" -fastaout $Sequences_MPL"_clustal_aligned_all_profile_trimmed_length_filt.fasta" -minseqlength 99

## only want seqs that have gdd

motifs_MPL='GDD'

python filter_sequences_by_motif.py $Sequences_MPL"_clustal_aligned_all_profile_trimmed_length_filt.fasta" $motifs_MPL $Translated_seqs_dir/temp_mpl.fasta $Translated_seqs_dir/Total_MPL_garbage.fasta

## make onto 1 line. 
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $Translated_seqs_dir/temp_mpl.fasta > $Translated_seqs_dir/temp_mpl_1L.fasta

mv $Translated_seqs_dir/temp_mpl_1L.fasta $Translated_seqs_dir/Total_MPL_pear.assembled_lines_annotated_quality_length_filt.fasta


#sed -i 's/\;/_/g;s/\+//g' $Sequences_MPL"_clustal_aligned_all_profile_trimmed_length_filt.fasta"

set +ue
deactivate
