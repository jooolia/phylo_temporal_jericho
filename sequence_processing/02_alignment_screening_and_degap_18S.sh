
set -o nounset -o pipefail -o errexit
 set +u 
##Fix naming!!! to long and indicipherable!
source /home/labop/Data/Julia/Julia_python/bin/activate
set -u
## note change resuls_small_test to the regular "results" when analysis is further along

## Change here for different directories
Sequence_directory=../results/Divided_by_primers/All_libs_combined
#Sequence_directory=../results_subsampled_1000/Divided_by_primers/All_libs_combined


Database_dir=../../Taxonomic_databases

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

## testing out the new vsearch for things requiring large amounts of memory
vsearch=/home/labop/vsearch/vsearch-1.0.6-linux-x86_64

python_script_dir=/home/labop/drive5_py
#python_script_dir=/home/julia/usearch/

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

Silva_gold=$Database_dir/silva.gold.align
#Gold_data=$Database_dir/gold_v_microbiomeutil-r20110519.fa

Gold_data=/Data/blastdb/gold_v_microbiomeutil-r20110519.fa
############ 18s #################

#commented out once so that I could use uchime. 


	## first convert to fasta, then subsample, then move:
	#zipped_assembled_18s_fastq=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.fastq.gz
	assembled_18s_fastq=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.fastq
	assembled_18s_fasta=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.fasta

	 #zcat $zipped_assembled_18s_fastq > $assembled_18s_fastq

	# convert to fasta
	$usearch7 -fastq_filter $assembled_18s_fastq  -fastaout $assembled_18s_fasta -fastq_ascii 33


## after first checking where the sequences hit I used this command (was down lower, but now I have moved it up)
# using the custom cut database will save memory and time! 
#mothur "#pcr.seqs(fasta=$Database_dir/silva.seed_v119.align, start=37000, end=43116, keepdots=F, processors=8)"

#mv $Database_dir/silva.seed_v119.pcr.align $Database_dir/silva.seed_v119_cut_for_18s




# ##re-running some parts of this
$mothur1_3 "#count.seqs(name=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.fasta); align.seqs(fasta=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.fasta, reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"


$mothur1_3 "#screen.seqs(fasta=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.align, optimize=start, end=5546, maxhomop=8, processors=8); filter.seqs(fasta=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.good.align, vertical=T)"


$mothur1_3 "#summary.seqs(fasta=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.fasta,processors=8)"


$mothur1_3 "#pcr.seqs(fasta=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.fasta, start=71, end=671, processors=8)"

$mothur1_3 "#degap.seqs(fasta=$Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.fasta)"

## gives $Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng.fasta ###



#rm 	$assembled_18s_fastq
#rm	$assembled_18s_fasta
gzip -f $Sequence_directory/Total_18s_pear.assembled_lines_annotated_quality_trim.align


python send_email_python.py $0

set +u +e
deactivate

