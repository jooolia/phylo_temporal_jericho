

##source /home/labop/Data/Julia/Julia_python/bin/activate


#set -o nounset -o pipefail -o errexit


## Programs

usearch8=/home/julia_g/usearch/usearch8.1.1861_i86linux32
#usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32


#hmm_build=hmmbuild
hmm_build=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild 

seqtk=/home/julia_g/seqtk/seqtk
#seqtk=seqtk

#mothur1_3=/home/labop/mothur/mothur
mothur1_3=/home/julia_g/mothur-1.36.1/mothur/mothur

Database_dir=../../Taxonomic_databases


reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Total_16s_forwardgood.good.pcr.ng.fasta
## seems to be some confusions with this....
miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ngnon_chimeras_ref97.00.fasta
miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_R1_annotated_trimmed_paired.good.good.ngnon_chimeras_ref97.00_sed.fasta
normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv
normalized_otu_table_top_10=../data/OTU_table_Jericho_time_series_16s_R1_normalized_top_10.tsv
Times_series=../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.fasta
Times_series_top_10=../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series_top_10.fasta
old_alignment=../../JerichoTimeSeries/results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.fasta

cyano_otus=../../JerichoTimeSeries/results/normalized_16s_cyano_otus.csv


head $cyano_otus -n 1 | sed s/\"//g |sed 's/\,/\n/g' > ../results/cyano_OTU_names.list
cat ../results/cyano_OTU_names.list
#	echo $OTU > ../results/$OTU".list"

sed s/\.size\.*\.//g  $miseq_amplicons_sed > ../results/amplicons_no_size_16s.fasta

$seqtk subseq ../results/amplicons_no_size_16s.fasta ../results/cyano_OTU_names.list > ../results/cyano_OTU_like.fasta

cut -f 3 -d , ../results/S16_VC_number_with_pool.csv |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/16s_pool_numbers.txt

## Need to subsample by libraries 
for pool in `cat ../results/16s_pool_numbers.txt `
do
echo $pool
number_reads=`grep -c $pool $reads`
echo $number_reads
if (("$number_reads" > "19000"))
then
echo "yay"
grep $pool $reads -A 1 --no-group-separator | $seqtk sample -s 100 - 19000 >> ../results/16s_reads_pools_above_5000_reads.fasta
fi
	# ## so if above certain number I want to take a random subset of the reads. 
done


#blastx -query $MPL_reads -db ../results/MPL_groupA_OTU_like -evalue 1e-3 -outfmt 6 -out ../results/blast_mpl_A_to_nucleotide.txt -num_threads 2 -max_target_seqs 1

#$usearch8 -usearch_global $MPL_reads -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.95 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta
#cp $reads ../results/current.fasta

$usearch8 -usearch_global ../results/16s_reads_pools_above_5000_reads.fasta -db ../results/cyano_OTU_like.fasta -strand both -id 0.97 -matched ../results/cyano_OTU_like_nucleotide.fasta

	## derep to quickly aling


	## derep to quickly aling

$usearch8 -derep_fulllength ../results/cyano_OTU_like_nucleotide.fasta -sizeout -fastaout ../results/derep_cyano_OTU_like_nucleotide.fasta -topn 4000

	## align using previous alignment

$mothur1_3 "#count.seqs(name=../results/derep_cyano_OTU_like_nucleotide.fasta); align.seqs(fasta=../results/derep_cyano_OTU_like_nucleotide.fasta, reference=$Database_dir/silva.seed_v119.align, flip=t, processors=8); summary.seqs(fasta=current)"	

$mothur1_3 "#filter.seqs(fasta=../results/derep_cyano_OTU_like_nucleotide.align, vertical=T)"

	## align using previous alignment




## just jericho
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/cyano_OTU_like_nucleotide.fasta > ../results/cyano_OTU_like_nucleotide_1L.fasta

	## search for the reads within those vcs that are within the time series. 
grep -A 1  --no-group-separator -f ../results/16s_pool_numbers.txt ../results/cyano_OTU_like_nucleotide_1L.fasta >  ../results/cyano_OTU_like_nucleotide_1L_just_jericho.fasta


sed -i 's/\./_/g' ../results/cyano_OTU_like_nucleotide_1L_just_jericho.fasta


$mothur1_3 "# align.seqs(fasta=../results/cyano_OTU_like_nucleotide_1L_just_jericho.fasta, reference=$Database_dir/silva.seed_v119.align, flip=t, processors=8);summary.seqs(fasta=current)"

#$mothur1_3 "# align.seqs(fasta=../results/cyano_OTU_like_nucleotide_1L_just_jericho.fasta, reference=../results/derep_cyano_OTU_like_nucleotide.filter_edit.fasta, flip=t, processors=8);summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/cyano_OTU_like_nucleotide_1L_just_jericho.align, minlength=150, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"




sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good.align | sed 's/_//g' | sed 's/test/_/g' > ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta 	

head ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta

entropy-analysis ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta  --quick


 oligotype ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -c 5 -M 50 -o ../results/oligotyping_cyano_otus/ --gen-html --skip-blast-search

## 98 97 21
 oligotype ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta ../results/cyano_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -C 38,39,61,66,81,112,121,122,130,161,97,21 -M 50 -o ../results/oligotyping_cyano_otus/ --gen-html --skip-blast-search
### 111, 254, 302 base locations of interest. 

######################################################
########### Do over whole time series
#####################################################

cut -f 3 -d , ../results/S16_VC_number_with_pool_with_all_times.csv |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/16s_pool_numbers_all_times.txt

## Need to subsample by libraries 
for pool in `cat ../results/16s_pool_numbers_all_times.txt `
do
echo $pool
number_reads=`grep -c $pool $reads`
echo $number_reads
if (("$number_reads" > "19000"))
then
echo "yay"
grep $pool $reads -A 1 --no-group-separator | $seqtk sample -s 100 - 19000 >> ../results/16s_reads_all_times_pools_above_5000_reads.fasta
fi
	# ## so if above certain number I want to take a random subset of the reads. 
done


#blastx -query $MPL_reads -db ../results/MPL_groupA_OTU_like -evalue 1e-3 -outfmt 6 -out ../results/blast_mpl_A_to_nucleotide.txt -num_threads 2 -max_target_seqs 1

#$usearch8 -usearch_global $MPL_reads -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.95 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta
#cp $reads ../results/current.fasta

$usearch8 -usearch_global ../results/16s_reads_all_times_pools_above_5000_reads.fasta -db ../results/cyano_OTU_like.fasta -strand both -id 0.97 -matched ../results/cyano_OTU_like_nucleotide_all_times.fasta


# $usearch8 -derep_fulllength ../results/cyano_OTU_like_nucleotide.fasta -sizeout -fastaout ../results/derep_cyano_OTU_like_nucleotide_all_times.fasta -topn 4000

# 	## align using previous alignment

# $mothur1_3 "#count.seqs(name=../results/derep_cyano_OTU_like_nucleotide.fasta); align.seqs(fasta=../results/derep_cyano_OTU_like_nucleotide.fasta, reference=../../JerichoTimeSeries/results/cyano_OTU_like.align, flip=t, processors=8); summary.seqs(fasta=current)"	

# $mothur1_3 "#filter.seqs(fasta=../results/derep_cyano_OTU_like_nucleotide.align, vertical=T)"

# 	## align using previous alignment


## just jericho
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/cyano_OTU_like_nucleotide_all_times.fasta > ../results/cyano_OTU_like_nucleotide_all_times_1L.fasta

	## search for the reads within those vcs that are within the time series. 
grep -A 1  --no-group-separator -f ../results/16s_pool_numbers_all_times.txt ../results/cyano_OTU_like_nucleotide_all_times_1L.fasta >  ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.fasta


sed -i 's/\./_/g' ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.fasta


$mothur1_3 "# align.seqs(fasta=../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.fasta, reference=$Database_dir/silva.seed_v119.align, flip=t, processors=8);summary.seqs(fasta=current)"

#$mothur1_3 "# align.seqs(fasta=../results/cyano_OTU_like_nucleotide_1L_just_jericho.fasta, reference=../results/derep_cyano_OTU_like_nucleotide.filter_edit.fasta, flip=t, processors=8);summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.align, minlength=150, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"

sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good.align | sed 's/_//g' | sed 's/test/_/g' > ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta 	

head ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta

entropy-analysis ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta  --quick


 oligotype ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -c 5 -M 50 -o ../results/oligotyping_cyano_otus_all_times/ --gen-html --skip-blast-search

## 98 97 21
 oligotype ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta ../results/cyano_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -C 38,39,61,66,81,112,121,122,130,161,97,21,212,211 -M 50 -o ../results/oligotyping_cyano_otus_all_times/ --gen-html --skip-blast-search
### 111, 254, 302 base locations of interest. 






python send_email_python.py $0

set +u +e
deactivate