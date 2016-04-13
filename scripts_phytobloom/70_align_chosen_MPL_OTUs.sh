## get top 10 OTUs, pull out reads, derep, align and then edit those alignments


source /home/labop/Data/Julia/Julia_python/bin/activate


set -o nounset -o pipefail -o errexit


## Programs
## Programs

usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32
#usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32

clustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64
#clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64
#hmm_build=hmmbuild
hmm_build=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild 
hmm_align=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmalign
esl_reformat=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-reformat
seqtk=/home/julia_g/seqtk/seqtk
#seqtk=seqtk

#mothur1_3=/home/labop/mothur/mothur
mothur1_3=/home/julia_g/mothur-1.36.1/mothur

for file in MPL
# gp23 
do

echo $file
if [ "$file" = "gp23" ]
then
Identity=0.95
Percent_ID=`echo "$Identity* 100" | bc`
echo $Percent_ID

# reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profile_trimmed_length_filt.fasta
# miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_annotated_quality_length_filt_otus_95.00.fasta
# miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_lines_annotated_quality_length_filt_otus_95.00_sed.fasta


# normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_gp23.tsv
normalized_otu_table_top_10=../data/OTU_table_Jericho_time_series_gp23_normalized_top_10.tsv
# Times_series=../results/Total_gp23_filtered_otus_95.00_for_Time_series.fasta
# Times_series_top_10=../results/Total_gp23_filtered_otus_95.00_for_Time_series_top_10.fasta
# old_alignment=../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta 

elif [ "$file" = "MPL" ]
then
Identity=0.95
Percent_ID=`echo "$Identity* 100" | bc`
echo $Percent_ID

# reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profile_trimmed_length_filt.fasta
# miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95.00.fasta
# miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95.00_sed.fasta


# normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_MPL.tsv
normalized_otu_table_top_10=../data/OTU_table_Jericho_time_series_MPL_normalized_top_10.tsv
# Times_series=../results/Total_MPL_filtered_otus_95.00_for_Time_series.fasta
# Times_series_top_10=../results/Total_MPL_filtered_otus_95.00_for_Time_series_top_10.fasta
 old_alignment=../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta 

fi

# sed '/^>/s/;/\./g;/^>/s/=/\./g ' $miseq_amplicons > $miseq_amplicons_sed 

# python filter_fasta_based_on_normalzed_data.py $miseq_amplicons_sed $normalized_otu_table $Times_series

#cp $reads ../results/current_reads.fasta

#$hmm_build --informat afa ../results/current.hmm $old_alignment 

#head $normalized_otu_table_top_10 -n 1 | sed s/\"//g > top_10_OTUs.txt
	for OTU in OTU_1 OTU_544 OTU_283 OTU_1588
	do

	$hmm_build --informat afa ../results/"derep_"$file"_"$OTU"like_reads_aligned_edit.hmm" ../results/"derep_"$file"_"$OTU"like_reads_aligned_edit.fasta" 



	## what if could pull in reads only for project:
	cut -f 3 -d , ../results/$file"_VC_number_with_pool.csv" |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/$file"_pool_numbers.txt"

	cut -f 3 -d , ../results/$file"_VC_number_with_pool_with_all_times.csv" |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/$file"_pool_numbers_all_times.txt"

	echo "pool17redo_" >> ../results/$file"_pool_numbers_all_times.txt"
	echo "pool16redo_" >> ../results/$file"_pool_numbers_all_times.txt"

	cat ../results/$file"_pool_numbers.txt"
	cat ../results/$file"_pool_numbers_all_times.txt"

	## change to a single line...
	## make onto 1 line. 
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/$file"_"$OTU"like_reads.fasta" > ../results/$file"_"$OTU"like_reads_1L.fasta"

	## search for the reads within those vcs that are within the time series. 
	grep -A 1  --no-group-separator -f ../results/$file"_pool_numbers.txt" ../results/$file"_"$OTU"like_reads_1L.fasta" >  ../results/$file"_"$OTU"like_reads_only_Jericho.fasta"

	grep -A 1  --no-group-separator -f ../results/$file"_pool_numbers_all_times.txt" ../results/$file"_"$OTU"like_reads_1L.fasta" >  ../results/$file"_"$OTU"like_reads_all_times.fasta"
		
	$hmm_align -o ../results/$file"_"$OTU"like_reads_aligned.stock" ../results/"derep_"$file"_"$OTU"like_reads_aligned_edit.hmm" ../results/$file"_"$OTU"like_reads_only_Jericho.fasta" 

	$esl_reformat --informat stockholm -o ../results/$file"_"$OTU"like_reads_aligned.fasta" afa ../results/$file"_"$OTU"like_reads_aligned.stock" 


	$hmm_align -o ../results/$file"_"$OTU"like_reads_aligned_all_times.stock" ../results/"derep_"$file"_"$OTU"like_reads_aligned_edit.hmm" ../results/$file"_"$OTU"like_reads_all_times.fasta" 

	$esl_reformat --informat stockholm -o ../results/$file"_"$OTU"like_reads_aligned_all_times.fasta" afa ../results/$file"_"$OTU"like_reads_aligned_all_times.stock" 

	#echo "done with aligning dereplicated seqs for " $file
	done
done