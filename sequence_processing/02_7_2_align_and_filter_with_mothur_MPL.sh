source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit


Sequence_directory=../results/Divided_by_primers/All_libs_combined


usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur


  ## first convert to fasta, then subsample, then move:
  zipped_assembled_MPL_fastq=$Sequence_directory/Total_MPL_pear.assembled_lines_annotated_quality_trim.fastq.gz
  assembled_MPL_fastq=$Sequence_directory/Total_MPL_pear.assembled_lines_annotated_quality_trim.fastq
  assembled_MPL_fasta=$Sequence_directory/Total_MPL_pear.assembled_lines_annotated_quality_trim.fasta

  ##### change to zcat for speed!!!

   # unzip fastq
  # gzip -d $zipped_assembled_MPL_fastq

#zcat $zipped_assembled_MPL_fastq > $assembled_MPL_fastq

## could I use zcat here??
  # convert to fasta
  $usearch7 -fastq_filter $assembled_MPL_fastq  -fastaout $assembled_MPL_fasta -fastq_ascii 33 -fastq_minlen 300


files=$assembled_MPL_fasta
basename $files
short_name=`basename -s .fasta $files`

$mothur1_3 "#count.seqs(name=$files, processors=8); 
align.seqs(fasta=$files, reference=../results/MPL_environmental_refs_nuc_derep_MAFFT_edited.fasta, flip=t, processors=8);
summary.seqs(fasta=current)"


python send_email_python.py "check MPL alignmetn"
read -p "Please edit 02_7_2_align_and_filter_with_mothur_MPL.sh to reflect appropriate screening"


$mothur1_3 "#screen.seqs(fasta=$Sequence_directory/$short_name".align", start=100, end=500, minlength=300, maxhomop=8, processors=8,maxambig=0);
 summary.seqs(fasta=$Sequence_directory/$short_name".good.align");
  filter.seqs(fasta=$Sequence_directory/$short_name".good.align", vertical=T);
  summary.seqs(fasta=current) "

#aliview $Sequence_directory/$short_name".good.only_1000"$direction".fasta"

$mothur1_3 "#degap.seqs(fasta=$Sequence_directory/$short_name".good.align")"

## gives ../results/Divided_by_primers/All_libs_combined/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng.fasta

set +u +e
deactivate