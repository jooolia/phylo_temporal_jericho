
source /home/labop/Data/Julia/Julia_python/bin/activate

################ need to fix all of this! Now do chimera checking after OTU picking....
set -o nounset -o pipefail -o errexit


## Change here for different directories
#Sequence_directory=../results_subsampled_1000/Divided_by_primers/All_libs_combined
Sequence_directory=../results/Divided_by_primers/All_libs_combined

Database_dir=../../Taxonomic_databases

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

## testing out the new vsearch for things requiring large amounts of memory
vsearch=/home/labop/vsearch/vsearch-1.0.6-linux-x86_64

python_script_dir=/home/labop/drive5_py
#python_script_dir=/home/julia/usearch/

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

#Gold_data=$Database_dir/gold_v_microbiomeutil-r20110519.fa
Gold_data=/Data/blastdb/gold_v_microbiomeutil-r20110519.fa



### should probably increase the sequence length to 350. 
$mothur1_3 "#count.seqs(name=$Sequence_directory/Total_16s_concat_annotated.fasta, processors=8);screen.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.fasta, minlength=300, maxhomop=8, processors=8); summary.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.fasta); align.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.fasta, reference=$Database_dir/silva.seed_v119.align, flip=t, processors=8); summary.seqs(fasta=current)"


echo " All together 16s"

#$mothur1_3 "#screen.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.align, start=6388, end=26970, maxhomop=8,maxambig=0, processors=8, minlength=400);summary.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.good.align, processors=8)"

$mothur1_3 "#screen.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.align, optimize=start, optimize=end, maxhomop=8,maxambig=0, processors=8, minlength=300);summary.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.good.align, processors=8)"


## remove gaps so that I can use with uparse
$mothur1_3 "#degap.seqs(fasta=$Sequence_directory/Total_16s_concat_annotated.good.good.align)"

rm $Sequence_directory/Total_16s_concat_annotated.good.good.align
gzip -f $Sequence_directory/Total_16s_concat_annotated.good.align
