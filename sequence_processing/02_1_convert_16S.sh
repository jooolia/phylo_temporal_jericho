

################ need to fix all of this! Now do chimera checking after OTU picking....

source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit


## Change here for different directories
#Sequence_directory=../results_subsampled_1000/Divided_by_primers/All_libs_combined
Sequence_directory=../results/Divided_by_primers/All_libs_combined

Database_dir=../../Taxonomic_databases

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

vsearch=/home/labop/vsearch/vsearch-1.0.6-linux-x86_64

fastx_q_to_a=/home/labop/fastx_bin/fastq_to_fasta

python_script_dir=/home/labop/drive5_py
#python_script_dir=/home/julia/usearch/

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

#Gold_data=$Database_dir/gold_v_microbiomeutil-r20110519.fa
Gold_data=/Data/blastdb/gold_v_microbiomeutil-r20110519.fa

####### Try with 16s #######################################

### Organize sequences first so that I can process them all into an alignment ###

## now should have all of the sequence here. 
	zipped_16s_fastq_R1=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.fastq.gz
	R1_16s_fastq=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.fastq
	R1_16s_fasta=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.fasta
   # zcat $zipped_16s_fastq_R1 > $R1_16s_fastq

	# convert to fasta
#$usearch7 -fastq_filter $R1_16s_fastq  -fastaout $R1_16s_fasta -fastq_ascii 33
#$fastx_q_to_a -i $R1_16s_fastq -o $R1_16s_fasta 

## too big for usearch

sed -n '1~4s/^@/>/p;2~4p' $R1_16s_fastq> $R1_16s_fasta


	zipped_16s_fastq_R2=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.fastq.gz
	R2_16s_fastq=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.fastq
	R2_16s_fasta=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.fasta
    #zcat $zipped_16s_fastq_R2 > $R2_16s_fastq

	# convert to fasta
#$usearch7 -fastq_filter $R2_16s_fastq  -fastaout $R2_16s_fasta -fastq_ascii 33
sed -n '1~4s/^@/>/p;2~4p' $R2_16s_fastq  > $R2_16s_fasta


python send_email_python.py $0

set +u +e
deactivate