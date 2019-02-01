

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

python_script_dir=/home/labop/drive5_py
#python_script_dir=/home/julia/usearch/

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

#Gold_data=$Database_dir/gold_v_microbiomeutil-r20110519.fa
Gold_data=/Data/blastdb/gold_v_microbiomeutil-r20110519.fa



####################################################
#### Try using mothur to do alignments and screen alignments ####


## then will cut the database and do all

## after examining this ut the database down to a length that is useful. 
#$mothur1_3 "#pcr.seqs(fasta=$Database_dir/silva.seed_v119.align, start=6000, end=32000, keepdots=F, processors=8)"

#mv $Database_dir/silva.seed_v119.pcr.align $Database_dir/silva.seed_v119_cut_for_16s


 echo " Align Concatenated R1 and R2 16s"


### should probably increase the sequence length to 350. 
#$mothur1_3 "#count.seqs(name=$Sequence_directory/Total_16s_concat_annotated.fasta, processors=8);screen.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.fasta, minlength=350, maxhomop=8, processors=8); summary.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.fasta); align.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.fasta, reference=$Database_dir/silva.seed_v119.align, flip=t, processors=8); summary.seqs(fasta=current)" 


## Also do the R1 and R2 so I can compare to the concatenated
echo " Align R1 16s"
R1_16s_fasta=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.fasta

$mothur1_3 "#count.seqs(name=$R1_16s_fasta, processors=8);screen.seqs(fasta=$R1_16s_fasta, minlength=150, maxhomop=8, processors=8); summary.seqs(fasta=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.good.fasta); align.seqs(fasta=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.good.fasta, reference=$Database_dir/silva.seed_v119.align, flip=t, processors=8); summary.seqs(fasta=current)"



 echo " Align R2 16s"
R2_16s_fasta=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.fasta

$mothur1_3 "#count.seqs(name=$R2_16s_fasta, processors=8);screen.seqs(fasta=$R2_16s_fasta, minlength=150, maxhomop=8, processors=8); summary.seqs(fasta=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.good.fasta); align.seqs(fasta=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.good.fasta, reference=$Database_dir/silva.seed_v119.align, flip=t, processors=8); summary.seqs(fasta=current)"

python send_email_python.py $0

set +u +e
deactivate