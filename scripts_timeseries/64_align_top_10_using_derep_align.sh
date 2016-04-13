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
hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 
hmm_align=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmalign
esl_reformat=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-reformat
seqtk=/home/julia/seqtk/seqtk
#seqtk=seqtk

#mothur1_3=/home/labop/mothur/mothur
mothur1_3=/home/julia/mothur/mothur

for file in MPL gp23 
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

head $normalized_otu_table_top_10 -n 1 | sed s/\"//g > top_10_OTUs.txt
	for OTU in `cat top_10_OTUs.txt`
	do

	$hmm_build --informat afa ../results/"derep_"$file"_"$OTU"like_reads_aligned_edit.hmm" ../results/"derep_"$file"_"$OTU"like_reads_aligned_edit.fasta" 



	## what if could pull in reads only for project:
	cut -f 3 -d , ../results/$file"_VC_number_with_pool.csv" |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/$file"_pool_numbers.txt"

	echo "pool17redo_" >> ../results/$file"_pool_numbers.txt"
	echo "pool16redo_" >> ../results/$file"_pool_numbers.txt"

	cat ../results/$file"_pool_numbers.txt"

	## change to a single line...
	## make onto 1 line. 
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/$file"_"$OTU"like_reads.fasta" > ../results/$file"_"$OTU"like_reads_1L.fasta"

	## search for the reads within those vcs that are within the time series. 
	grep -A 1  --no-group-separator -f ../results/$file"_pool_numbers.txt" ../results/$file"_"$OTU"like_reads_1L.fasta" >  ../results/$file"_"$OTU"like_reads_only_Jericho.fasta"
		
	$hmm_align -o ../results/$file"_"$OTU"like_reads_aligned.stock" ../results/"derep_"$file"_"$OTU"like_reads_aligned_edit.hmm" ../results/$file"_"$OTU"like_reads_only_Jericho.fasta" 

	$esl_reformat --informat stockholm -o ../results/$file"_"$OTU"like_reads_aligned.fasta" afa ../results/$file"_"$OTU"like_reads_aligned.stock" 


	#echo "done with aligning dereplicated seqs for " $file
	done
done


for file in 16s 18s 
do

echo $file
if [ "$file" = "16s" ]
then 
Identity=0.97
Percent_ID=`echo "$Identity* 100" | bc`
echo $Percent_ID

# reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Total_16s_forwardgood.good.pcr.ng.fasta
# miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ngnon_chimeras_ref97.00.fasta
# miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_R1_annotated_trimmed_paired.good.good.ngnon_chimeras_ref97.00_sed.fasta
# normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv
normalized_otu_table_top_10=../data/OTU_table_Jericho_time_series_16s_R1_normalized_top_10.tsv
# Times_series=../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.fasta
# Times_series_top_10=../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series_top_10.fasta
old_alignment=../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.fasta

elif [ "$file" = "18s" ]
then
Identity=0.97
Percent_ID=`echo "$Identity* 100" | bc`
echo $Percent_ID

# reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng.fasta
# miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.fasta
# miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.ng_otus_97.00_sed.fasta
# normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_18S.tsv
normalized_otu_table_top_10=../data/OTU_table_Jericho_time_series_18s_normalized_top_10.tsv
# Times_series=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.fasta
# Times_series_top_10=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series_top_10.fasta
old_alignment=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.fasta
fi


# sed '/^>/s/;/\./g;/^>/s/=/\./g ' $miseq_amplicons > $miseq_amplicons_sed 

# python filter_fasta_based_on_normalzed_data.py $miseq_amplicons_sed $normalized_otu_table $Times_series

#cp $reads ../results/current_reads.fasta

#$hmm_build --informat afa ../results/current.hmm $old_alignment 

head $normalized_otu_table_top_10 -n 1 | sed s/\"//g > top_10_OTUs.txt
	for OTU in `cat top_10_OTUs.txt`
	do
	## or could even try to pull all o

	echo $OTU > ../results/$OTU".list"

	# 	$seqtk subseq $Times_series ../results/$OTU".list" > ../results/$file"_"$OTU"like.fasta"

	# 	## sometimes filenames too long

	# 	$usearch8 -usearch_global ../results/current_reads.fasta -db ../results/$file"_"$OTU"like.fasta" -strand both -id $Identity -matched ../results/$file"_"$OTU"like_reads.fasta"

	# 	## derep to quickly aling

	# 	$usearch8 -derep_fulllength ../results/$file"_"$OTU"like_reads.fasta" -sizeout -fastaout ../results/"derep_"$file"_"$OTU"like_reads.fasta"

	# 	## align using previous alignment

	# $mothur1_3 "#align.seqs(fasta=../results/"derep_"$file"_"$OTU"like_reads.fasta", reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"	

	# $mothur1_3 "#filter.seqs(fasta=../results/"derep_"$file"_"$OTU"like_reads.align", vertical=T)"

	cut -f 3 -d , ../results/$file"_VC_number_with_pool.csv" |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/$file"_pool_numbers.txt"

	echo "pool17redo_" >> ../results/$file"_pool_numbers.txt"
	echo "pool16redo_" >> ../results/$file"_pool_numbers.txt"

	cat ../results/$file"_pool_numbers.txt"

	## search for the reads within those vcs that are within the time series. 
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/$file"_"$OTU"like_reads.fasta" > ../results/$file"_"$OTU"like_reads_1L.fasta"

	## search for the reads within those vcs that are within the time series. 
	grep -A 1  --no-group-separator -f ../results/$file"_pool_numbers.txt" ../results/$file"_"$OTU"like_reads_1L.fasta" >  ../results/$file"_"$OTU"like_reads_only_Jericho.fasta"

	head ../results/$file"_"$OTU"like_reads_only_Jericho.fasta"
	head $old_alignment
	# ../results/derep_16s_OTU_1.size.748217.like_reads.filter.fasta
		
 	# $mothur1_3 "#align.seqs(fasta=../results/$file"_"$OTU"like_reads_only_Jericho.fasta", reference=../results/"derep_"$file"_"$OTU"like_reads.filter_edit.fasta", flip=t, processors=8); summary.seqs(fasta=current)"	
 	$mothur1_3 "#align.seqs(fasta=../results/$file"_"$OTU"like_reads_only_Jericho.fasta", reference=$old_alignment, flip=t, processors=8)"	

 	$mothur1_3 "#filter.seqs(fasta=../results/$file"_"$OTU"like_reads_only_Jericho.align", vertical=T)"

	# ../results/derep_16s_OTU_1.size.748217.like_reads.filter.fasta
	mv ../results/$file"_"$OTU"like_reads_only_Jericho.filter.fasta" ../results/$file"_"$OTU"like_reads_aligned.fasta"	

	#echo "done with aligning dereplicated seqs for " $file
	done
done


python send_email_python.py $0

set +u +e
deactivate