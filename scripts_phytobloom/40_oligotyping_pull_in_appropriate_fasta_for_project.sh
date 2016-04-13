
## Idea is to oligotype the heterosigma and the dinoflagellate otus that dominate in the bloom and see if there are some dynamics that I'm able to pull out. 

## 1 Pull out appropriate reads from the data, so that would be at 97% using usearch for the specific otus

## find OTUS of interest


## pull out reads using usearch
## use oligotyping to separate reads into types
## observe types over time. 


## Programs
Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated


hmm_build=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild 
hmm_align=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmalign
esl_reformat=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/esl-reformat
seqtk=/home/julia_g/seqtk/seqtk
usearch8=/home/julia_g/usearch/usearch8.1.1861_i86linux32
#usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32


#hmm_build=hmmbuild
hmm_build=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild 

seqtk=/home/julia_g/seqtk/seqtk
#seqtk=seqtk

#mothur1_3=/home/labop/mothur/mothur
mothur1_3=/home/julia_g/mothur-1.36.1/mothur/mothur


S18_reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng.fasta

S18_miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.fasta

S18_miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.ng_otus_97.00_sed.fasta

normalized_otu_table=../data/OTU_table_bloom_time_series_normalized_18S.tsv

normalized_otu_table_top_100=../data/OTU_table_bloom_time_series_18s_normalized_top_100.tsv

Times_series_18s=../results/Total_18s_filtered_otus_97.00_for_bloom_Time_series.fasta

Times_series_18s_top_100=../results/Total_18s_filtered_otus_97.00_for_bloom_Time_series_top_100.fasta

Database_dir=../../Taxonomic_databases

sed '/^>/s/;/\./g;/^>/s/=/\./g ' $S18_miseq_amplicons > $S18_miseq_amplicons_sed 

python filter_fasta_based_on_normalzed_data.py $S18_miseq_amplicons_sed $normalized_otu_table $Times_series_18s


#python filter_fasta_based_on_normalzed_data.py $S18_miseq_amplicons_sed $normalized_otu_table_top_100 $Times_series_18s_top_100


## which OTU is dominant on 1252? OTU 2 which is a Heterosigma classification. 

## pull out this OTU
for pool in `cat ../results/S18_pool_numbers.txt `
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





echo "OTU_2.size.159373." > ../results/haka.list

$seqtk subseq $Times_series_18s ../results/haka.list > ../results/heterosigma_otu.fasta

$usearch8 -usearch_global ../results/18s_reads_pools_above_10000_reads.fasta -db ../results/heterosigma_otu.fasta -strand both -id 0.97 -matched ../results/heterosigma_reads.fasta
		
# ##re-running some parts of this

## what about dereplicating and then using that aligment for aligning all the reads?? quick try...

$usearch8 -derep_fulllength ../results/heterosigma_reads.fasta -sizeout -fastaout derep_haka.fa

#$mothur1_3 "#count.seqs(name=derep_haka.fa); align.seqs(fasta=derep_haka.fa, reference=../../Taxonomic_databases/silva.seed_v119.align, flip=t, processors=8); summary.seqs(fasta=current)"	


#$mothur1_3 "#filter.seqs(fasta=derep_haka.align, vertical=T)"


## edit the dereplicated read alignment. 
#aliview derep_haka.filter.fasta
## aliview derep_haka.filter_edit.fasta

$mothur1_3 "#count.seqs(name=../results/heterosigma_reads.fasta); align.seqs(fasta=../results/heterosigma_reads.fasta, reference=derep_haka.align, flip=t, processors=8); summary.seqs(fasta=current)"		
 
#aliview ../results/heterosigma_reads.align  
#aliview ../results/heterosigma_reads_edit.align

#$mothur1_3 "#filter.seqs(fasta=../results/heterosigma_reads.align, vertical=T)"

#aliview ../results/heterosigma_reads.filter.fasta
## aliview ../results/heterosigma_reads.filter_edit.fasta


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/heterosigma_reads.align  > temp1.fasta
sed 's/_//g' temp1.fasta > temp2.fasta
sed 's/test/_/g' temp2.fasta > ../results/heterosigma_reads_new_names.fasta 



entropy-analysis ../results/heterosigma_reads_new_names.fasta --quick

oligotype ../results/heterosigma_reads_new_names.fasta ../results/heterosigma_reads_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/haka_only_bloom/ --gen-html --skip-blast-search



## dinoflagellate
echo "OTU_3.size.86173." > ../results/dino.list

$seqtk subseq $Times_series_18s ../results/dino.list > ../results/dino_otu.fasta

$usearch8 -usearch_global ../results/18s_reads_pools_above_10000_reads.fasta -db ../results/dino_otu.fasta -strand both -id 0.97 -matched ../results/dino_reads.fasta
		
# ##re-running some parts of this

## what about dereplicating and then using that aligment for aligning all the reads?? quick try...

$usearch8 -derep_fulllength ../results/dino_reads.fasta -sizeout -fastaout derep_dino.fa

$mothur1_3 "#count.seqs(name=derep_haka.fa); align.seqs(seed=clear,fasta=derep_dino.fa, reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"	


$mothur1_3 "#filter.seqs(fasta=derep_dino.align, vertical=T)"


## edit the dereplicated read alignment. 
aliview derep_dino.filter.fasta
## aliview derep_dino.filter_edit.fasta

$mothur1_3 "#align.seqs(fasta=../results/dino_reads.fasta, reference=derep_dino.filter.fasta, flip=t, processors=8); summary.seqs(fasta=current)"		
 
#aliview ../results/dino_reads.align  
#aliview ../results/dino_reads_edit.align

#$mothur1_3 "#filter.seqs(fasta=../results/dino_reads_edit.align, vertical=T)"

#aliview ../results/heterosigma_reads.filter.fasta
## aliview ../results/heterosigma_reads.filter_edit.fasta

## need to format the names
sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/dino_reads.align > temp1.fasta
sed 's/_//g' temp1.fasta > temp2.fasta
sed 's/test/_/g' temp2.fasta > ../results/dino_reads_new_names.fasta 


## Use library sequencing info to identify which pool corresponds to which VC
#Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=83)
cut -f 3 -d , ../results/S18_VC_number_with_pool.csv |tail -n +2 | sed 's/^/pool/' >../results/S18_pool_numbers.txt

grep -A 1 -f  ../results/S18_pool_numbers.txt ../results/dino_reads_new_names.fasta >  ../results/dino_reads_new_names_only_bloom.fasta

## weirdly there were some extra empty columns at the end of the reads. Opened in aliview and got rid of that

### fix

entropy-analysis ../results/dino_reads_new_names_only_bloom.fasta --quick

oligotype ../results/dino_reads_new_names_only_bloom.fasta ../results/dino_reads_new_names_only_bloom.fasta-ENTROPY -c 3 -M 50 -o ../results/dino_only_bloom/ --gen-html --skip-blast-search



#########################################
### The second most abundant oligotype for dinos
########################################

echo "OTU_48.size.6381." > ../results/dino2.list


$seqtk subseq $Times_series_18s ../results/dino2.list > ../results/dino2_otu.fasta

$usearch8 -usearch_global ../results/18s_reads_pools_above_10000_reads.fasta -db ../results/dino2_otu.fasta -strand both -id 0.97 -matched ../results/dino2_reads.fasta
		
# ##re-running some parts of this

## what about dereplicating and then using that aligment for aligning all the reads?? quick try...

$usearch8 -derep_fulllength ../results/dino2_reads.fasta -sizeout -fastaout derep_dino2.fa

$mothur1_3 "# align.seqs(seed=clear,fasta=derep_dino2.fa, reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"	


$mothur1_3 "#filter.seqs(fasta=derep_dino2.align, vertical=T)"


## edit the dereplicated read alignment. 
#Saliview derep_dino2.filter.fasta
## aliview derep_haka.filter_edit.fasta

$mothur1_3 "#count.seqs(name=../results/dino2_reads.fasta); align.seqs(fasta=../results/dino2_reads.fasta, reference=derep_dino2.filter.fasta, flip=t, processors=8); summary.seqs(fasta=current)"		
 
#aliview ../results/dino2_reads.align  
#aliview ../results/dino_reads_edit.align

#$mothur1_3 "#filter.seqs(fasta=../results/dino_reads_edit.align, vertical=T)"

#aliview ../results/heterosigma_reads.filter.fasta
## aliview ../results/heterosigma_reads.filter_edit.fasta

## need to format the names
sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/dino2_reads.align > temp1.fasta
sed 's/_//g' temp1.fasta > temp2.fasta
sed 's/test/_/g' temp2.fasta > ../results/dino2_reads_new_names.fasta 


entropy-analysis ../results/dino2_reads_new_names.fasta  --quick

oligotype ../results/dino2_reads_new_names.fasta ../results/dino2_reads_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/dino2_only_bloom/ --gen-html --skip-blast-search




########################################
################ do for all time points.


## dinoflagellate
echo "OTU_3.size.86173." > ../results/dino.list

$seqtk subseq $Times_series_18s ../results/dino.list > ../results/dino_otu.fasta

$usearch8 -usearch_global ../results/18s_all_times_reads_pools_above_10000_reads.fasta -db ../results/dino_otu.fasta -strand both -id 0.97 -matched ../results/dino_reads_all_times.fasta
		
# ##re-running some parts of this

## what about dereplicating and then using that aligment for aligning all the reads?? quick try...

# $usearch8 -derep_fulllength ../results/dino_reads.fasta -sizeout -fastaout derep_dino.fa

# $mothur1_3 "#count.seqs(name=derep_haka.fa); align.seqs(fasta=derep_dino.fa, reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"	


# $mothur1_3 "#filter.seqs(fasta=derep_dino.align, vertical=T)


$mothur1_3 "#align.seqs(seed=clear,fasta=../results/dino_reads_all_times.fasta, reference=../results/dino_reads_edit.align, flip=t, processors=8); summary.seqs(fasta=current)"		
 
$mothur1_3 "#filter.seqs(fasta=../results/dino_reads_all_times.align, vertical=T)"

#aliview ../results/heterosigma_reads.filter.fasta
## aliview ../results/heterosigma_reads.filter_edit.fasta

## need to format the names
sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/dino_reads_all_times.filter.fasta > temp1.fasta
sed 's/_//g' temp1.fasta > temp2.fasta
sed 's/test/_/g' temp2.fasta > ../results/dino_reads_new_names_all_times.fasta 


## Use library sequencing info to identify which pool corresponds to which VC
#Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=83)
# cut -f 3 -d , ../results/S18_VC_number_with_pool.csv |tail -n +2 | sed 's/^/pool/' >../results/S18_pool_numbers.txt

grep -A 1 -f  ../results/S18_pool_numbers_all_times.txt ../results/dino_reads_new_names_all_times.fasta >  ../results/dino_reads_new_names_Jericho_all_times.fasta

## weirdly there were some extra empty columns at the end of the reads. Opened in aliview and got rid of that

### fix

entropy-analysis ../results/dino_reads_new_names_Jericho_all_times.fasta --quick

oligotype ../results/dino_reads_new_names_Jericho_all_times.fasta ../results/dino_reads_new_names_Jericho_all_times.fasta-ENTROPY -c 3 -M 50 -o ../results/dino_all_times/ --gen-html --skip-blast-search

# -C 76,85,214
#########################################
### The second most abundant oligotype for dinos
########################################

echo "OTU_24948.size.1." > ../results/dino2.list


$seqtk subseq $Times_series_18s ../results/dino2.list > ../results/dino2_otu.fasta

$usearch8 -usearch_global ../results/18s_all_times_reads_pools_above_10000_reads.fasta -db ../results/dino2_otu.fasta -strand both -id 0.97 -matched ../results/dino2_reads_all_times.fasta
		
# ##re-running some parts of this

## what about dereplicating and then using that aligment for aligning all the reads?? quick try...

 $usearch8 -derep_fulllength ../results/dino2_reads.fasta -sizeout -fastaout ../results/derep_dino2.fa

 # $mothur1_3 "# align.seqs(fasta=derep_dino2.fa, reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"	


# $mothur1_3 "#filter.seqs(fasta=derep_dino2.align, vertical=T)"


## edit the dereplicated read alignment. 
# aliview derep_dino2.filter.fasta
## aliview derep_haka.filter_edit.fasta

#$mothur1_3 "#align.seqs(seed=clear,fasta=../results/dino2_reads_all_times.fasta, reference=../results/derep_dino2_edit.align, flip=t, processors=8); summary.seqs(fasta=current)"		
 
$mothur1_3 "#align.seqs(seed=clear,fasta=../results/dino2_reads_all_times.fasta, reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"



#aliview ../results/dino2_reads.align  
#aliview ../results/dino_reads_edit.align

#$mothur1_3 "#filter.seqs(fasta=../results/dino_reads_edit.align, vertical=T)"

#aliview ../results/heterosigma_reads.filter.fasta
## aliview ../results/heterosigma_reads.filter_edit.fasta

## need to format the names
sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/dino2_reads_all_times.align > temp1.fasta
sed 's/_//g' temp1.fasta > temp2.fasta
sed 's/test/_/g' temp2.fasta > ../results/dino2_reads_all_times_new_names.fasta 


entropy-analysis ../results/dino2_reads_all_times_new_names.fasta   --quick

# oligotype ../results/dino2_reads_all_times_new_names.fasta  ../results/dino2_reads_all_times_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/dino2_all_times/ --gen-html --skip-blast-search

oligotype ../results/dino2_reads_all_times_new_names.fasta  ../results/dino2_reads_all_times_new_names.fasta-ENTROPY -C 70,84,102,79,87,215,214 -M 50 -o ../results/dino2_all_times/ --gen-html --skip-blast-search


