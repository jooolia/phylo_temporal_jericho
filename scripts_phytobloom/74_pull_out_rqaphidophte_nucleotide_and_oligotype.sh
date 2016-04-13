

#source /home/labop/Data/Julia/Julia_python/bin/activate


set -o nounset -o pipefail 
#-o errexit

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

reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng.fasta
miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.fasta
miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.ng_otus_97.00_sed.fasta
normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_18S.tsv
normalized_otu_table_top_10=../data/OTU_table_Jericho_time_series_18s_normalized_top_10.tsv
Times_series=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.fasta
Times_series_top_10=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series_top_10.fasta
old_aligment=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.fasta
raphidophyte_otus=../../JerichoTimeSeries/results/normalized_18s_raphidophyte_otus.csv


head $raphidophyte_otus -n 1 | sed s/\"//g |sed 's/\,/\n/g' > ../results/raphidophyte_OTU_names.list
cat ../results/raphidophyte_OTU_names.list
#	echo $OTU > ../results/$OTU".list"

sed s/\.size\.*\.//g  $miseq_amplicons_sed > ../results/amplicons_no_size_18s.fasta

$seqtk subseq ../results/amplicons_no_size_18s.fasta ../results/raphidophyte_OTU_names.list > ../results/raphidophyte_OTU_like.fasta


## what if could pull in reads only for project:
cut -f 3 -d , ../results/S18_VC_number_with_pool.csv |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/S18_pool_numbers.txt

#echo "pool17redo_" >> ../results/$file"_pool_numbers.txt"
#echo "pool16redo_" >> ../results/$file"_pool_numbers.txt"

cat ../results/S18_pool_numbers.txt

rm ../results/18s_reads_pools_above_10000_reads.fasta

## Need to subsample by libraries 
for pool in `cat ../results/S18_pool_numbers.txt`
do
echo $pool
number_reads=`grep -c $pool $reads`
echo $number_reads
if (("$number_reads" > "10000"))
then
echo "yay"
grep $pool $reads -A 1 --no-group-separator | $seqtk sample -s 100 - 10000 >> ../results/18s_reads_pools_above_10000_reads.fasta
fi

	# ## so if above certain number I want to take a random subset of the reads. 
done
	

#blastx -query $MPL_reads -db ../results/MPL_groupA_OTU_like -evalue 1e-3 -outfmt 6 -out ../results/blast_mpl_A_to_nucleotide.txt -num_threads 2 -max_target_seqs 1

#$usearch8 -usearch_global $MPL_reads -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.95 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta
#cp $reads ../results/current.fasta
echo "usearching!"

$usearch8 -usearch_global ../results/18s_reads_pools_above_10000_reads.fasta -db ../results/raphidophyte_OTU_like.fasta -strand both -id 0.97 -matched ../results/raphidophyte_OTU_like_nucleotide.fasta


$usearch8 -derep_fulllength ../results/raphidophyte_OTU_like_nucleotide.fasta -sizeout -fastaout ../results/derep_raphidophyte_OTU_like_nucleotide.fasta -topn 4000

	## align using previous alignment

$mothur1_3 "#count.seqs(name=../results/derep_raphidophyte_OTU_like_nucleotide.fasta); align.seqs(fasta=../results/derep_raphidophyte_OTU_like_nucleotide.fasta, reference=../../JerichoTimeSeries/results/derep_raphidophyte_OTU_like_nucleotide.filter_edit.fasta, flip=t, processors=8); summary.seqs(fasta=current)"	

$mothur1_3 "#filter.seqs(fasta=../results/derep_raphidophyte_OTU_like_nucleotide.align, vertical=T)"


## just jericho
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/raphidophyte_OTU_like_nucleotide.fasta > ../results/raphidophyte_OTU_like_nucleotide_1L.fasta

	## search for the reads within those vcs that are within the time series. 
grep -A 1  --no-group-separator -f ../results/S18_pool_numbers.txt ../results/raphidophyte_OTU_like_nucleotide_1L.fasta >  ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.fasta


$mothur1_3 "#count.seqs(name=../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.fasta, processors=8); 
align.seqs(fasta=../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.fasta, reference=../results/derep_raphidophyte_OTU_like_nucleotide.filter.fasta, flip=t, processors=8);
summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.align, minlength=150, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good.align | sed 's/_//g' | sed 's/test/_/g' > ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta 	

head ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta 		

entropy-analysis ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta  --quick

oligotype ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta  ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/oligotyping_raphidophyte_OTUs/ --gen-html --skip-blast-search


### 111, 254, 302 base locations of interest. 
oligotype ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta  ../results/raphidophyte_OTU_like_nucleotide_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -C 218,88,109,47,213,148,108 -M 50 -o ../results/oligotyping_raphidophyte_OTUs/ --gen-html --skip-blast-search



######################################
### do also with all times

cut -f 3 -d , ../results/S18_VC_number_with_pool_with_all_times.csv |tail -n +2 | sed 's/^/pool/' | sed 's/\"//g' | sed 's/pool\([0-9]\+\)/pool\1_/g' > ../results/S18_pool_numbers_all_times.txt

#echo "pool17redo_" >> ../results/$file"_pool_numbers.txt"
#echo "pool16redo_" >> ../results/$file"_pool_numbers.txt"

cat ../results/S18_pool_numbers_all_times.txt


## Need to subsample by libraries 
for pool in `cat ../results/S18_pool_numbers_all_times.txt `
do
echo $pool
number_reads=`grep -c $pool $reads`
echo $number_reads
if (("$number_reads" > "10000"))
then
echo "yay"
grep $pool $reads -A 1 --no-group-separator | $seqtk sample -s 100 - 10000 >> ../results/18s_all_times_reads_pools_above_10000_reads.fasta
fi
	# ## so if above certain number I want to take a random subset of the reads. 
done


#blastx -query $MPL_reads -db ../results/MPL_groupA_OTU_like -evalue 1e-3 -outfmt 6 -out ../results/blast_mpl_A_to_nucleotide.txt -num_threads 2 -max_target_seqs 1

#$usearch8 -usearch_global $MPL_reads -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.95 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta
#cp $reads ../results/current.fasta

$usearch8 -usearch_global ../results/18s_all_times_reads_pools_above_10000_reads.fasta -db ../results/raphidophyte_OTU_like.fasta -strand both -id 0.97 -matched ../results/raphidophyte_OTU_like_nucleotide_all_times.fasta


# $usearch8 -derep_fulllength ../results/raphidophyte_OTU_like_nucleotide_all_times.fasta -sizeout -fastaout ../results/derep_raphidophyte_OTU_like_nucleotide.fasta -topn 4000

	## align using previous alignment

# $mothur1_3 "#count.seqs(name=../results/derep_raphidophyte_OTU_like_nucleotide.fasta); align.seqs(fasta=../results/derep_raphidophyte_OTU_like_nucleotide.fasta, reference=../../JerichoTimeSeries/results/derep_raphidophyte_OTU_like_nucleotide.filter_edit.fasta, flip=t, processors=8); summary.seqs(fasta=current)"	

# $mothur1_3 "#filter.seqs(fasta=../results/derep_raphidophyte_OTU_like_nucleotide.align, vertical=T)"


## just jericho
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/raphidophyte_OTU_like_nucleotide_all_times.fasta > ../results/raphidophyte_OTU_like_nucleotide_all_times_1L.fasta

	## search for the reads within those vcs that are within the time series. 
grep -A 1  --no-group-separator -f ../results/S18_pool_numbers_all_times.txt ../results/raphidophyte_OTU_like_nucleotide_all_times_1L.fasta >  ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.fasta


$mothur1_3 "#count.seqs(name=../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.fasta, processors=8); 
align.seqs(fasta=../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.fasta, reference=../results/derep_raphidophyte_OTU_like_nucleotide.filter.fasta, flip=t, processors=8);
summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.align, minlength=150, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good.align | sed 's/_//g' | sed 's/test/_/g' > ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta 	

head ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta 		

entropy-analysis ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta  --quick

oligotype ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta  ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/oligotyping_raphidophyte_all_times_OTUs/ --gen-html --skip-blast-search


### 111, 254, 302 base locations of interest. 
oligotype ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta  ../results/raphidophyte_OTU_like_nucleotide_all_times_1L_just_jericho.good_edit_new_names.fasta-ENTROPY -C 218,88,109,47,213,148,108 -M 50 -o ../results/oligotyping_raphidophyte_all_times_OTUs/ --gen-html --skip-blast-search


#python send_email_python.py $0

set +u +e
deactivate