source /home/labop/Data/Julia/Julia_python/bin/activate

### align with mafft a set of nucleotides in the hopes of making database for mothur filtering of gp20
## trying to do this to get rid of so many problematic sequences that have been found in the sequences...

set -o nounset -o pipefail -o errexit

## or could I find if from my searches already where I downloaded genbank files??

Sequence_directory=../results/Divided_by_primers/All_libs_combined


mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

################ 16s ###########################
  for reads in 1 2
do


  ################### R1        ####################
  #16s_zipped=$Sequence_directory/Total_16s_R$reads"_annotated_trimmed_paired.fastq.gz"
  #16s_fastq=../results/Divided_by_primers/All_libs_combined/Total_16s_R$reads"_annotated_trimmed_paired.fastq"
S16_fasta=../results/Divided_by_primers/All_libs_combined/Total_16s_R$reads"_annotated_trimmed_paired.fasta"


  
  #rm $16s_fastq

basename $S16_fasta
short_name=`basename -s .fasta $S16_fasta`

gunzip $Sequence_directory/$short_name".good.align"

## starting with $Sequence_directory/Total_16s_R2_annotated_trimmed_paired.good.align
#../results/Divided_by_primers/All_libs_combined/Total_16s_R1_annotated_trimmed_paired.good.align

#### 
echo "########## Screen seqs to divide into forward and reverse primers for 16s then, could concatenate"


echo "16s forward"

start_pos=6443
end_pos=15697
direction="forward"

$mothur1_3 "#screen.seqs(fasta=$Sequence_directory/$short_name".good.align", start=$start_pos, end=$end_pos, maxhomop=8,maxambig=0, processors=8, minlength=190);summary.seqs(fasta=$Sequence_directory/$short_name".good.good.align", processors=8)"


## remove gaps so that I can use with uparse
#$mothur1_3 "#degap.seqs(fasta=$Sequence_directory/$short_name".good.good.align")"
mv $Sequence_directory/$short_name".good.good.align" $Sequence_directory/$short_name"_"$direction"good.good.align"
# mv $Sequence_directory/$short_name".good.good.ng.fasta" $Sequence_directory/$short_name"_"$direction"_no_gaps.fasta"

echo " 16s Reverse"

start_pos=21910
end_pos=27163
direction="reverse"

$mothur1_3 "#screen.seqs(fasta=$Sequence_directory/$short_name".good.align", start=$start_pos, end=$end_pos, maxhomop=8,maxambig=0, processors=8, minlength=190);summary.seqs(fasta=$Sequence_directory/$short_name".good.good.align", processors=8)"


# #; filter.seqs(fasta=$Sequence_directory/$short_name".good.align", vertical=T);summary.seqs(fasta=current) "

# #aliview $Sequence_directory/$short_name".good.align"

#$mothur1_3 "#degap.seqs(fasta=$Sequence_directory/$short_name".good.good.align")"
mv $Sequence_directory/$short_name".good.good.align" $Sequence_directory/$short_name"_"$direction"good.good.align"
# mv $Sequence_directory/$short_name".good.good.ng.fasta" $Sequence_directory/$short_name"_"$direction"_no_gaps.fasta"

gzip $Sequence_directory/$short_name".good.align"
python send_email_python.py $Sequence_directory/$short_name".good.align"


done

set +u +e
deactivate