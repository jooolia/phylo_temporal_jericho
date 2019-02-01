
##Fix naming!!! to long and indicipherable!
source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit
## could think about analyzing the chloroplast data separately. 

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32

python_script_dir=/home/labop/drive5_py
#python_script_dir=/home/julia/usearch/

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

vsearch=/home/labop/vsearch/vsearch-1.0.6-linux-x86_64

## note change resuls_small_test to the regular "results" when analysis is further along
Sequence_directory=../results/Divided_by_primers/All_libs_combined

Database_dir=../../Taxonomic_databases

Translated_seqs_dir=$Sequence_directory/Translated

### OTU picking with uclust

Picked_OTUs=$Translated_seqs_dir/Picked_OTUs

mkdir -p $Picked_OTUs

Gold_data=/Data/blastdb/gold_v_microbiomeutil-r20110519.fa
gp23_ref=/Data/blastdb/G23dbP.fasta
MPL_ref=/Data/RDRPdbP



########################## 16s and 18s ########################
## Sequences from script 02 aftder trimming and filtering with mothur

echo =========================================================================
echo " Time for 16s and 18s"
echo =========================================================================


Identity=0.97
Percent_ID=`echo "$Identity* 100" | bc`


Sequences_16s=$Sequence_directory/Total_16s_concat_annotated.good.good.ng.fasta
Sequences_16s_R1=$Sequence_directory/Total_16s_forwardgood.good.pcr.ng.fasta
#Sequences_16s_R1=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.good.good.filter.ng.fasta ## confusion here. Which one..

#Sequences_16s_R1=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.good.good.ng.fasta
Sequences_16s_R2=../results/Divided_by_primers/All_libs_combined/Total_16s_reversegood.good.pcr.ng.fasta
#Sequences_16s_R2=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.good.good.ng.fasta
Sequences_18s=../results/Divided_by_primers/All_libs_combined/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng.fasta

for file in $Sequences_18s $Sequences_16s_R1 $Sequences_16s_R2 $Sequences_16s
do 

	echo =========================================================================
	echo " Dereplicate and chimera check sequences"
	echo =========================================================================

	basename_reads=$(basename $file .fasta)
	echo $basename_reads


	## classify and the remove...
	#$mothur1_3 "#classify.seqs(fasta=$file, reference=$Database_dir/silva.nr_v119.align, taxonomy=$Database_dir/silva.nr_v119.tax, iters=1000, processors=8)"

	
	echo =========================================================================
	echo " OTU picking at 97%"
	echo =========================================================================

	$vsearch -derep_fulllength $file -sizeout -output derep.fa --minseqlength 100

	$vsearch -sortbysize derep.fa -sizeout -output derep.fas -minsize 1
	# # #         #~ # Make sure that usearch is not stopping to early as the defaulf for max reject is 25. 
	count_lines_for_max_reject=`grep -c '^>' derep.fas`
	echo "number for max reject"
	echo $count_lines_for_max_reject

	# ## does not support -threads
	$usearch8 -cluster_otus derep.fas -strand both -id $Identity -otus $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -maxrejects $count_lines_for_max_reject -uparseout $Picked_OTUs/$basename_reads"_"$Percent_ID"results.txt" -relabel OTU_ -sizein -sizeout

	#$usearch8 -uchime_denovo $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -chimeras $Picked_OTUs/$basename_reads"chimeras.fasta" -nonchimeras $Picked_OTUs/$basename_reads"non_chimeras_denovo"$Percent_ID".fasta" -uchimeout $Picked_OTUs/$basename_reads"centroids_"$Percent_ID".uchime" -abskew 2

	$vsearch -uchime_ref $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -db $Gold_data -strand plus -uchimeout $Picked_OTUs/$basename_reads.uchime -chimeras $Picked_OTUs/$basename_reads"chimeras_ref.fasta" -nonchimeras $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta"


	# echo =========================================================================
	# echo "Take the intersect of the two chimera methods"
	# echo =========================================================================


	 #python get_intersection_of_non_chimeric_sequences.py $Picked_OTUs/$basename_reads"non_chimeras_denovo"$Percent_ID".fasta" $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta" $Picked_OTUs/$basename_reads"combined_unique_non_chimeras_"$Percent_ID".fasta"



	 echo =========================================================================
	 echo "Classify sequences using Silva database"
	 echo =========================================================================

	 $mothur1_3 "#classify.seqs(fasta=$Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta", reference=$Database_dir/silva.nr_v119.align, taxonomy=$Database_dir/silva.nr_v119.tax, iters=1000, processors=8)"

	# ## gives file with .good.filter.ng_outs_97.00.nr_v119.taxonomy
	# ## remove the chloroplast reads

    ## this had reduced then number of reads that I got, so will keep so that my reads/sample are better and then remove the total number!!  
	$mothur1_3 "#remove.lineage(taxonomy=$Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".nr_v119.wang.taxonomy", taxon=Chloroplast-Mitochondria-unknown-Archaea, fasta=$Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta")"

	#python $python_script_dir/fasta_number.py $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta" OTU_ > $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_numbered.fasta" 

	$usearch8 -usearch_global $file -db $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".pick.fasta" -strand plus -id $Identity -uc $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_global.uc" 

	# ## on server moved my own script to my own directory
	python $python_script_dir/uc2otutab_JAG.py $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_global.uc" > $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"global_OTU_table.tsv"


	# ## generates $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".pick.fasta

	# ## not sure if I need to remove the lineage from the numbered otus too
	# python Mapping_taxonomy_to_OTUs.py $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta"\
	#  $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_numbered.fasta"\
	#  $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".nr_v119.wang.taxonomy"\
	#  $Picked_OTUs/$basename_reads"OTUs_"$Percent_ID"taxonomy_and_otus_mapped.csv"


	# 		rm derep.fa
	# 		rm derep.fas

	python send_email_python.py $0

done

set +ue
deactivate

