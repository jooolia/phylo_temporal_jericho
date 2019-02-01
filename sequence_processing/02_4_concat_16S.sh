

################ need to fix all of this! Now do chimera checking after OTU picking....

source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit


## Change here for different directories
#Sequence_directory=../results_subsampled_1000/Divided_by_primers/All_libs_combined
Sequence_directory=../results/Divided_by_primers/All_libs_combined

Database_dir=../../Taxonomic_databases

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32
usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32

vsearch=/home/labop/vsearch/vsearch-1.0.6-linux-x86_64

fastx_q_to_a=/home/labop/fastx_bin/fastq_to_fasta

python_script_dir=/home/labop/drive5_py
#python_script_dir=/home/julia/usearch/

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

#Gold_data=$Database_dir/gold_v_microbiomeutil-r20110519.fa
Gold_data=/Data/blastdb/gold_v_microbiomeutil-r20110519.fa

Sequences_16s_R1=$Sequence_directory/Total_16s_R1_annotated_trimmed_paired.good.good.ng.fasta
Sequences_16s_R2=$Sequence_directory/Total_16s_R2_annotated_trimmed_paired.good.good.ng.fasta



Sequences_16s_R1_F=../results/Divided_by_primers/All_libs_combined/Total_16s_R1_annotated_trimmed_paired_forwardgood.good.pcr.ng.fasta
Sequences_16s_R1_R=../results/Divided_by_primers/All_libs_combined/Total_16s_R1_annotated_trimmed_paired_reversegood.good.pcr.ng.fasta
Sequences_16s_R2_F=../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired_forwardgood.good.pcr.ng.fasta
Sequences_16s_R2_R=../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired_reversegood.good.pcr.ng.fasta

echo "concatenate based on the primer forward or reverse"
cat $Sequences_16s_R1_F $Sequences_16s_R2_F > ../results/Divided_by_primers/All_libs_combined/Total_16s_forwardgood.good.pcr.ng.fasta
cat $Sequences_16s_R1_R $Sequences_16s_R2_R > ../results/Divided_by_primers/All_libs_combined/Total_16s_reversegood.good.pcr.ng.fasta

echo "concatenate based on the read R1 or R2"
cat $Sequences_16s_R1_F $Sequences_16s_R1_R > ../results/Divided_by_primers/All_libs_combined/Total_16s_R1good.good.pcr.ng.fasta
cat $Sequences_16s_R2_F $Sequences_16s_R2_R > ../results/Divided_by_primers/All_libs_combined/Total_16s_R2good.good.pcr.ng.fasta


#python concatenate_R1_and_R2_for_non_merging_primers.py ../results/Divided_by_primers/All_libs_combined/Total_16s_R1good.good.pcr.ng.fasta ../results/Divided_by_primers/All_libs_combined/Total_16s_R2good.good.pcr.ng.fasta $Sequence_directory/Total_16s_concat_annotated.fasta

## the alignment process would have already flipped the sequences so it is not necessary to flip the sequence for the concatenation. these sequences are also trimmed so that there should be no overlapping problems in the middle. 


echo "concatenate based on the read R1 forward and R2 reverse"
python concatenate_R1_and_R2_for_non_merging_primers_already_aligned_forward.py $Sequences_16s_R1_F $Sequences_16s_R2_R $Sequence_directory/Total_16s_concat_annotated_A.fasta

echo "concatenate based on the read R2 forward and R1 reverse"
python concatenate_R1_and_R2_for_non_merging_primers_already_aligned_reverse.py $Sequences_16s_R2_F $Sequences_16s_R1_R $Sequence_directory/Total_16s_concat_annotated_B.fasta

cat $Sequence_directory/Total_16s_concat_annotated_A.fasta $Sequence_directory/Total_16s_concat_annotated_B.fasta >  $Sequence_directory/Total_16s_concat_annotated.fasta



python send_email_python.py $0

set +u +e
deactivate