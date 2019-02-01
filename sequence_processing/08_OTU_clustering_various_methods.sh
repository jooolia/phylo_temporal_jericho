
##Fix naming!!! to long and indicipherable!
source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit
## could think about analyzing the chloroplast data separately. 

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

#usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32
usearch8=/home/labop/usearchv8.0.1623
#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32

python_script_dir=/home/labop/drive5_py
#python_script_dir=/home/julia/usearch/

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

vsearch=/home/labop/vsearch/vsearch-1.0.6-linux-x86_64

## note change resuls_small_test to the regular "results" when analysis is further along
Sequence_directory=../results/Divided_by_primers/All_libs_combined

Database_dir=../../Taxonomic_databases

gp23_ref=/Data/blastdb/G23db.fasta
MPL_ref=/Data/Julia/qiime/PicornaLikedb.fasta


Translated_seqs_dir=$Sequence_directory/Translated

Sequences_gp23=../results/Divided_by_primers/All_libs_combined/Translated/Total_gp23_pear.assembled_annotated_quality_length_filt.fasta 
Sequences_MPL=../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_length_filt.fasta 
#Sequences_MPL=$Translated_seqs_dir/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 


### OTU picking with uclust

Picked_OTUs=$Translated_seqs_dir/Picked_OTUs

mkdir -p $Picked_OTUs

Gold_data=/Data/blastdb/gold_v_microbiomeutil-r20110519.fa
gp23_ref=/Data/blastdb/G23dbP.fasta
MPL_ref=/Data/RDRPdbP


echo =========================================================================
echo " Picking OTUs for MPL and gp23"
echo =========================================================================


for file in $Sequences_MPL $Sequences_gp23  
do
echo $file
if [[ $( echo $file) == *gp23* ]]
then
Identity=0.95
Percent_ID=`echo "$Identity* 100" | bc`
database=$gp23_ref
echo $Percent_ID

elif [[ $( echo $file) == *MPL* ]]
then
Identity=0.95
Percent_ID=`echo "$Identity* 100" | bc`
database=$MPL_ref
echo $Percent_ID
fi
basename_reads=$(basename $file .fasta)
echo $basename_reads

		$usearch8 -derep_fulllength $file -sizeout -fastaout derep.fa

		$usearch8 -sortbysize derep.fa -fastaout derep.fas -minsize 1
		# Make sure that usearch is not stopping to early as the defaulf for max reject is 25. 
		count_lines_for_max_reject=`grep -c '^>' derep.fas`
		echo "number for max reject"
		echo $count_lines_for_max_reject

		
		$usearch8 -cluster_otus derep.fas -strand both -id $Identity -relabel OTU_ -otus $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -maxrejects $count_lines_for_max_reject -uparseout $Picked_OTUs/$basename_reads"_"$Percent_ID"results.txt"
## not working with protein sequences


        #$usearch8 -uchime_denovo $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -chimeras $Picked_OTUs/$basename_reads"chimeras_denovo.fasta" -nonchimeras $Picked_OTUs/$basename_reads"non_chimeras_denovo"$Percent_ID".fasta" -uchimeout $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".uchime" -abskew 2

        #$vsearch -uchime_ref $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -db $database -strand plus -uchimeout $Picked_OTUs/$basename_reads.uchime -chimeras $Picked_OTUs/$basename_reads"chimeras_ref.fasta" -nonchimeras $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta"




# echo =========================================================================
# echo "Take the intersect of the two chimera methods"
# echo =========================================================================



# python get_intersection_of_non_chimeric_sequences.py $Picked_OTUs/$basename_reads"non_chimeras_denovo"$Percent_ID".fasta" $Picked_OTUs/$basename_reads"non_chimeras_ref"$Percent_ID".fasta" $Picked_OTUs/$basename_reads"combined_unique_non_chimeras_"$Percent_ID".fasta"


		## original sequences to compare to:

		#python $python_script_dir/fasta_number.py $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" OTU_ > $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_numbered.fasta" 




		$usearch8 -usearch_global $file -db $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -strand plus -id $Identity -uc $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_global.uc" 
		#$vsearch -usearch_global $file -db $Picked_OTUs/$basename_reads"_otus_"$Percent_ID".fasta" -strand plus -id $Identity -uc $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_global.uc"

		## on server moved my own script to my own directory
		python $python_script_dir/uc2otutab_JAG.py $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"_global.uc" > $Picked_OTUs/$basename_reads"_otus_"$Percent_ID"global_OTU_table.tsv"

	done


set +ue
deactivate


