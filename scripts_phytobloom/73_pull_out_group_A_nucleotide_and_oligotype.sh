

set -o nounset  -o errexit
#-o pipefail
usearch8=/home/julia_g/usearch/usearch8.1.1861_i86linux32
#usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32


#hmm_build=hmmbuild
hmm_build=/home/julia_g/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmbuild 

seqtk=/home/julia_g/seqtk/seqtk
#seqtk=seqtk

#mothur1_3=/home/labop/mothur/mothur
mothur1_3=/home/julia_g/mothur-1.36.1/mothur/mothur
## good MPL
MPL_reads=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.fasta

miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95.00_sed.fasta
## make list of group A fasta
head ../../JerichoTimeSeries/results/proportional_MPL_groupA.csv -n 1 | sed s/\"//g |sed 's/\,/\n/g' > ../results/proportional_MPL_groupA_OTU_names.txt
cat ../results/proportional_MPL_groupA_OTU_names.txt
#	echo $OTU > ../results/$OTU".list"

$seqtk subseq $miseq_amplicons_sed ../results/proportional_MPL_groupA_OTU_names.txt > ../results/MPL_groupA_OTU_like.fasta


makeblastdb -in ../results/MPL_groupA_OTU_like.fasta -out ../results/MPL_groupA_OTU_like -dbtype prot -hash_index

#blastx -query $MPL_reads -db ../results/MPL_groupA_OTU_like -evalue 1e-3 -outfmt 6 -out ../results/blast_mpl_A_to_nucleotide.txt -num_threads 2 -max_target_seqs 1

#$usearch8 -usearch_global $MPL_reads -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.95 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta

#rm ../results/MPL_reads_pools_above_5000_reads.fasta

## Need to subsample by libraries 
for pool in `cat ../results/MPL_pool_numbers.txt `
do
echo $pool
number_reads=`grep -c $pool $MPL_reads`
echo $number_reads
if (("$number_reads" > "3300"))
then
echo "yay"
grep $pool $MPL_reads -A 1 --no-group-separator | $seqtk sample -s 100 - 3300 >> ../results/MPL_reads_pools_above_5000_reads.fasta
fi
	# ## so if above certain number I want to take a random subset of the reads. 
done


#### look for Group A

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


oligotype ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta  ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta-ENTROPY -c 1 -M 50 -o ../results/oligotyping_MPL_group_A_from_nuc_like_reads/ --gen-html --skip-blast-search

oligotype ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta  ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta-ENTROPY -C 5,121,122,266,299,311,314,269 -M 50 -o ../results/oligotyping_MPL_group_A_from_nuc_like_reads/ --gen-html --skip-blast-search
### 111, 254, 302 base locations of interest.
## or should it be 54? 



#####################################################################################
############ Try with OTUs over time that seem important ############################
### Of note seem to be 3, 5 ,4 and 288

echo OTU_3 > OTU3.list

$seqtk subseq $miseq_amplicons_sed OTU3.list > ../results/MPL_OTU_3.fasta

echo OTU_4 > OTU4.list

$seqtk subseq $miseq_amplicons_sed OTU4.list > ../results/MPL_OTU_4.fasta


echo OTU_5 > OTU5.list

$seqtk subseq $miseq_amplicons_sed OTU5.list > ../results/MPL_OTU_5.fasta

echo OTU_288 > OTU288.list

$seqtk subseq $miseq_amplicons_sed OTU288.list > ../results/MPL_OTU_288.fasta

makeblastdb -in ../results/MPL_OTU_3.fasta -out ../results/MPL_OTU_3 -dbtype prot -hash_index
makeblastdb -in ../results/MPL_OTU_4.fasta -out ../results/MPL_OTU_4 -dbtype prot -hash_index
makeblastdb -in ../results/MPL_OTU_5.fasta -out ../results/MPL_OTU_5 -dbtype prot -hash_index
makeblastdb -in ../results/MPL_OTU_288.fasta -out ../results/MPL_OTU_288 -dbtype prot -hash_index

#blastx -query $MPL_reads -db ../results/MPL_groupA_OTU_like -evalue 1e-3 -outfmt 6 -out ../results/blast_mpl_A_to_nucleotide.txt -num_threads 2 -max_target_seqs 1

#$usearch8 -usearch_global $MPL_reads -db ../results/MPL_groupA_OTU_like.fasta -strand both -id 0.95 -matched ../results/MPL_groupA_OTU_like_nucleotide.fasta

## create empty file
rm ../results/MPL_reads_pools_above_5000_reads.fasta
## Need to subsample by libraries 
for pool in `cat ../results/MPL_pool_numbers_all_times.txt `
do
echo $pool
number_reads=`grep -c $pool $MPL_reads`
echo $number_reads
if (("$number_reads" > "2500"))
then
echo "yay"
grep $pool $MPL_reads -A 1 --no-group-separator | $seqtk sample -s 100 - 2500 >> ../results/MPL_reads_pools_above_5000_reads.fasta
fi
	# ## so if above certain number I want to take a random subset of the reads. 
done


#### look for OTUs

$usearch8 -usearch_local ../results/MPL_reads_pools_above_5000_reads.fasta -db ../results/MPL_OTU_3.fasta -strand both -id 0.70 -matched ../results/MPL_OTU_3_nucleotide.fasta

$usearch8 -usearch_local ../results/MPL_reads_pools_above_5000_reads.fasta -db ../results/MPL_OTU_4.fasta -strand both -id 0.70 -matched ../results/MPL_OTU_4_nucleotide.fasta
$usearch8 -usearch_local ../results/MPL_reads_pools_above_5000_reads.fasta -db ../results/MPL_OTU_5.fasta -strand both -id 0.70 -matched ../results/MPL_OTU_5_nucleotide.fasta
$usearch8 -usearch_local ../results/MPL_reads_pools_above_5000_reads.fasta -db ../results/MPL_OTU_288.fasta -strand both -id 0.70 -matched ../results/MPL_OTU_288_nucleotide.fasta

# 	## derep to quickly aling

# $usearch8 -derep_fulllength ../results/MPL_groupA_OTU_like_nucleotide.fasta -sizeout -fastaout ../results/derep_MPL_groupA_OTU_like_nucleotide.fasta -topn 20000


# $mothur1_3 "#count.seqs(name=../results/derep_MPL_groupA_OTU_like_nucleotide.fasta, processors=8); 
# align.seqs(fasta=../results/derep_MPL_groupA_OTU_like_nucleotide.fasta, reference=../results/derep_MPL_groupA_OTU_like_nucleotide_edit.align, flip=t, processors=8);
# summary.seqs(fasta=current)"

## have a pretty decent alignment of nucleotide reads. 

###### here!!!!

##  so that they are on one line
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/MPL_OTU_3_nucleotide.fasta > ../results/MPL_OTU_3_nuc_1L.fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/MPL_OTU_4_nucleotide.fasta > ../results/MPL_OTU_4_nuc_1L.fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/MPL_OTU_5_nucleotide.fasta > ../results/MPL_OTU_5_nuc_1L.fasta
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ../results/MPL_OTU_288_nucleotide.fasta > ../results/MPL_OTU_288_nuc_1L.fasta
	# ## search for the reads within those vcs that are within the time series. 
	# grep -A 1  --no-group-separator -f ../results/MPL_pool_numbers.txt ../results/MPL_nuc_reads_1L.fasta >  ../results/MPL_nuc_reads_1Llike_reads_only_Jericho.fasta

### now do a bit more separately
##### OTU 3

## if I could separate this and keep it separate the would be great!
$usearch8 -derep_fulllength ../results/MPL_OTU_3_nuc_1L.fasta -sizeout -fastaout ../results/derep_MPL_OTU_3_nuc_1L.fasta -topn 2000

mafft ../results/derep_MPL_OTU_3_nuc_1L.fasta > ../results/derep_MPL_OTU_3_nuc.align

$mothur1_3 "#count.seqs(name=../results/MPL_OTU_3_nuc_1L.fasta, processors=8); 
align.seqs(fasta=../results/MPL_OTU_3_nuc_1L.fasta, reference=../results/derep_MPL_OTU_3_nuc_edit.align, flip=t, processors=8);
summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/MPL_OTU_3_nuc_1L.align, minlength=300, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"

python filter_seqs_with_high_percent_gaps.py ../results/MPL_OTU_3_nuc_1L.good.align ../results/MPL_OTU_3_nuc_1L.good.gappy_removed.fasta


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/MPL_OTU_3_nuc_1L.good.gappy_removed.fasta | sed 's/_//g' | sed 's/test/_/g' > ../results/MPL_OTU_3_nuc_edit_new_names.fasta 	

head ../results/MPL_OTU_3_nuc_edit_new_names.fasta	

entropy-analysis ../results/MPL_OTU_3_nuc_edit_new_names.fasta --quick

oligotype ../results/MPL_OTU_3_nuc_edit_new_names.fasta  ../results/MPL_OTU_3_nuc_edit_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/oligotyping_MPL_OTU_3_nuc_all_times/ --gen-html --skip-blast-search

oligotype ../results/MPL_OTU_3_nuc_edit_new_names.fasta  ../results/MPL_OTU_3_nuc_edit_new_names.fasta-ENTROPY -C 53,118,221,59,416,245,276,443,272 -M 50 -o ../results/oligotyping_MPL_OTU_3_nuc_all_times/ --gen-html --skip-blast-search


########### OTU 4 

$usearch8 -derep_fulllength ../results/MPL_OTU_4_nuc_1L.fasta -sizeout -fastaout ../results/derep_MPL_OTU_4_nuc_1L.fasta -topn 200

mafft ../results/derep_MPL_OTU_4_nuc_1L.fasta > ../results/derep_MPL_OTU_4_nuc.align

$mothur1_3 "#count.seqs(name=../results/MPL_OTU_4_nuc_1L.fasta, processors=8); 
align.seqs(fasta=../results/MPL_OTU_4_nuc_1L.fasta, reference=../results/derep_MPL_OTU_4_nuc_edit.align, flip=t, processors=8);
summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/MPL_OTU_4_nuc_1L.align, minlength=300, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"


python filter_seqs_with_high_percent_gaps.py ../results/MPL_OTU_4_nuc_1L.good.align ../results/MPL_OTU_4_nuc_1L.fasta.good.gappy_removed.fasta


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/MPL_OTU_4_nuc_1L.fasta.good.gappy_removed.fasta | sed 's/_//g' | sed 's/test/_/g' > ../results/MPL_OTU_4_nuc_edit_new_names.fasta 	

head ../results/MPL_OTU_4_nuc_edit_new_names.fasta	

entropy-analysis ../results/MPL_OTU_4_nuc_edit_new_names.fasta --quick


oligotype ../results/MPL_OTU_4_nuc_edit_new_names.fasta  ../results/MPL_OTU_4_nuc_edit_new_names.fasta-ENTROPY -c 1 -M 50 -o ../results/oligotyping_MPL_OTU_4_nuc_all_times/ --gen-html --skip-blast-search

#### OTu 5


$usearch8 -derep_fulllength ../results/MPL_OTU_5_nuc_1L.fasta -sizeout -fastaout ../results/derep_MPL_OTU_5_nuc_1L.fasta -topn 200

mafft ../results/derep_MPL_OTU_5_nuc_1L.fasta > ../results/derep_MPL_OTU_5_nuc.align

$mothur1_3 "#count.seqs(name=../results/MPL_OTU_5_nuc_1L.fasta, processors=8); 
align.seqs(fasta=../results/MPL_OTU_5_nuc_1L.fasta, reference=../results/derep_MPL_OTU_5_nuc_edit.align, flip=t, processors=8);
summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/MPL_OTU_5_nuc_1L.align, minlength=300, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"


python filter_seqs_with_high_percent_gaps.py ../results/MPL_OTU_5_nuc_1L.good_edit.align ../results/MPL_OTU_5_nuc_1L.fasta.good.gappy_removed.fasta


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/MPL_OTU_5_nuc_1L.fasta.good.gappy_removed.fasta | sed 's/_//g' | sed 's/test/_/g' > ../results/MPL_OTU_5_nuc_edit_new_names.fasta 	

head ../results/MPL_OTU_5_nuc_edit_new_names.fasta	

entropy-analysis ../results/MPL_OTU_5_nuc_edit_new_names.fasta --quick


oligotype ../results/MPL_OTU_5_nuc_edit_new_names.fasta  ../results/MPL_OTU_5_nuc_edit_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/oligotyping_MPL_OTU_5_nuc_all_times/ --gen-html --skip-blast-search

oligotype ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta  ../results/MPL_nuc_reads_1Llike_reads_only_Jericho_edit_new_names.fasta-ENTROPY -C 25, 90, 193
 -M 50 -o ../results/oligotyping_MPL_group_A_from_nuc_like_reads/ --gen-html --skip-blast-search



 ### OTU 288

$usearch8 -derep_fulllength ../results/MPL_OTU_288_nuc_1L.fasta -sizeout -fastaout ../results/derep_MPL_OTU_288_nuc_1L.fasta -topn 200

mafft ../results/derep_MPL_OTU_288_nuc_1L.fasta > ../results/derep_MPL_OTU_288_nuc.align

$mothur1_3 "#count.seqs(name=../results/MPL_OTU_288_nuc_1L.fasta, processors=8); 
align.seqs(fasta=../results/MPL_OTU_288_nuc_1L.fasta, reference=../results/derep_MPL_OTU_288_nuc_edit.align, flip=t, processors=8);
summary.seqs(fasta=current)"

$mothur1_3 "#screen.seqs(fasta=../results/MPL_OTU_288_nuc_1L.align, minlength=300, maxhomop=8, processors=8,maxambig=0); summary.seqs(fasta=current)"

python filter_seqs_with_high_percent_gaps.py ../results/MPL_OTU_288_nuc_1L.good_edit.align ../results/MPL_OTU_288_nuc_1L.fasta.good.gappy_removed.fasta


sed 's/Lib_pool\([0-9]\+\)/Lib-pool\1test/g' ../results/MPL_OTU_288_nuc_1L.fasta.good.gappy_removed.fasta | sed 's/_//g' | sed 's/test/_/g' > ../results/MPL_OTU_288_nuc_edit_new_names.fasta 	

head ../results/MPL_OTU_288_nuc_edit_new_names.fasta	

entropy-analysis ../results/MPL_OTU_288_nuc_edit_new_names.fasta --quick

oligotype ../results/MPL_OTU_288_nuc_edit_new_names.fasta  ../results/MPL_OTU_288_nuc_edit_new_names.fasta-ENTROPY -c 3 -M 50 -o ../results/oligotyping_MPL_OTU_288_nuc_all_times/ --gen-html --skip-blast-search


oligotype ../results/MPL_OTU_288_nuc_edit_new_names.fasta  ../results/MPL_OTU_288_nuc_edit_new_names.fasta-ENTROPY -C 17,132,321,317,200,167,140,272,171,164,81 -M 50 -o ../results/oligotyping_MPL_OTU_288_nuc_all_times/ --gen-html --skip-blast-search