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


echo "16s forward"


direction="forward"


$mothur1_3 "#pcr.seqs(fasta=$Sequence_directory/$short_name"_"$direction"good.good.align", start=6430, end=15697, processors=8);degap.seqs(fasta=$Sequence_directory/$short_name"_"$direction"good.good.pcr.align", processors=8)"


echo " 16s Reverse"
direction="reverse"

$mothur1_3 "#pcr.seqs(fasta=$Sequence_directory/$short_name"_"$direction"good.good.align", start=21791, end=28463, processors=8);degap.seqs(fasta=$Sequence_directory/$short_name"_"$direction"good.good.pcr.align",processors=8)"

#mv $Sequence_directory/$short_name".good.good.align" $Sequence_directory/$short_name"_"$direction"good.good.align"
#mv $Sequence_directory/$short_name".good.good.ng.fasta" $Sequence_directory/$short_name"_"$direction"_no_gaps.fasta"


python send_email_python.py $Sequence_directory/$short_name".good.good.ng.fasta"


done

set +u +e
deactivate