
## get top 10 OTUs, pull out reads, derep, align and then edit those alignments


source /home/labop/Data/Julia/Julia_python/bin/activate


set -o nounset -o pipefail -o errexit


## Programs
## Programs

#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32
usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32

clustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64
#clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64
hmm_build=hmmbuild
#hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 

#seqtk=/home/julia/seqtk/seqtk
seqtk=seqtk



reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profile_trimmed_length_filt.fasta
miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_annotated_quality_length_filt_otus_95.00.fasta
miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_lines_annotated_quality_length_filt_otus_95.00_sed.fasta


normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_gp23.tsv
normalized_otu_table_top_10=../data/OTU_table_Jericho_time_series_gp23_normalized_top_10.tsv
Times_series=../results/Total_gp23_filtered_otus_95.00_for_Time_series.fasta
Times_series_top_10=../results/Total_gp23_filtered_otus_95.00_for_Time_series_top_10.fasta
old_alignment=../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta 

# sed '/^>/s/;/\./g;/^>/s/=/\./g ' $miseq_amplicons > $miseq_amplicons_sed 

# python filter_fasta_based_on_normalzed_data.py $miseq_amplicons_sed $normalized_otu_table $Times_series

head ../results/proportional_gp23_groupI.csv -n 1 | sed s/\"//g |sed 's/\,/\n/g' > ../results/proportional_gp23_groupI_OTU_names.txt
cat ../results/proportional_gp23_groupI_OTU_names.txt
#	echo $OTU > ../results/$OTU".list"

$seqtk subseq $miseq_amplicons_sed ../results/proportional_gp23_groupI_OTU_names.txt > ../results/gp23_groupI_OTU_like.fasta


makeblastdb -in ../results/MPL_groupA_OTU_like.fasta -out ../results/MPL_groupA_OTU_like -dbtype prot -hash_index

#blastx -query $MPL_reads -db ../results/MPL_groupA_OTU_like -evalue 1e-3 -outfmt 6 -out ../results/blast_mpl_A_to_nucleotide.txt -num_threads 2 -max_target_seqs 1

#$usearch8 -usearch_global $MPL_reads -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.95 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta

## Need to subsample by libraries 
for pool in `cat ../results/MPL_pool_numbers.txt `
do
echo $pool
number_reads=`grep -c $pool $MPL_reads`
echo $number_reads
if (("$number_reads" > "5000"))
then
echo "yay"
grep $pool $MPL_reads -A 1 --no-group-separator | $seqtk sample -s 100 - 5000 >> ../results/MPL_reads_pools_above_5000_reads.fasta
fi
	# ## so if above certain number I want to take a random subset of the reads. 
done





$usearch8 -usearch_local ../results/MPL_reads_pools_above_5000_reads.fasta -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.70 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta

	## derep to quickly aling

$usearch8 -derep_fulllength ../results/MPL_groupA_OTU_like_nucleotide.fasta -sizeout -fastaout ../results/derep_MPL_groupA_OTU_like_nucleotide.fasta -topn 20000


$mothur1_3 "#count.seqs(name=../results/derep_MPL_groupA_OTU_like_nucleotide.fasta, processors=8); 
align.seqs(fasta=../results/derep_MPL_groupA_OTU_like_nucleotide.fasta, reference=../results/derep_MPL_groupA_OTU_like_nucleotide_edit.align, flip=t, processors=8);
summary.seqs(fasta=current)"

## have a pretty decent alignment of nucleotide reads. 

###### here!!!!

## just jericho
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/MPL_groupA_OTU_like_nucleotide.fasta > ../results/MPL_nuc_reads_1L.fasta

	## search for the reads within those vcs that are within the time series. 
	grep -A 1  --no-group-separator -f ../results/MPL_pool_numbers.txt ../results/MPL_nuc_reads_1L.fasta >  ../results/MPL_nuc_reads_1Llike_reads_only_Jericho.fasta


 $mothur1_3 "#count.seqs(name=../results/MPL_nuc_reads_1Llike_reads_only_Jericho.fasta, processors=8); 
align.seqs(fasta=../results/MPL_nuc_reads_1Llike_reads_only_Jericho.fasta, reference=../results/derep_MPL_groupA_OTU_like_nucleotide_edit.align, flip=t, processors=8);
summary.seqs(fasta=current)"


$mothur1_3 "#screen.seqs(fasta=../results/MPL_nuc_reads_1Llike_reads_only_Jericho.align, minlength=300, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"


python filter_seqs_with_high_percent_gaps.py ../results/MPL_nuc_reads_1Llike_reads_only_Jericho.good.align ../results/MPL_nuc_reads_1Llike_reads_only_Jericho.good.gappy_removed.fasta


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/MPL_nuc_reads_1Llike_reads_only_Jericho.good.gappy_removed.fasta | sed 's/_//g' | sed 's/test/_/g' > ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta 	

head ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta	

entropy-analysis ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta --quick



head $normalized_otu_table_top_10 -n 1 | sed s/\"//g > top_10_OTUs.txt
	for OTU in `cat top_10_OTUs.txt`
	do

	sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/$file"_"$OTU"like_reads_aligned.fasta" | sed 's/_//g' | sed 's/test/_/g' > ../results/$file"_"$OTU"like_reads_aligned_new_names.fasta" 	

	head ../results/$file"_"$OTU"like_reads_aligned_new_names.fasta"	

	entropy-analysis ../results/$file"_"$OTU"like_reads_aligned_new_names.fasta" --quick

	oligotype ../results/$file"_"$OTU"like_reads_aligned_new_names.fasta" ../results/$file"_"$OTU"like_reads_aligned_new_names.fasta-ENTROPY" -c 1 -M 50 -o ../results/"oligotyping"$file"_"$OTU"like_reads"/ --gen-html 

		#echo "done with aligning dereplicated seqs for " $file
	done
done

python send_email_python.py $0

set +u +e
deactivate




