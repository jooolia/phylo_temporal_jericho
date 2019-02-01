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

gp23_nt_database=/Data/blastdb/G23db.fasta

  ## first convert to fasta, then subsample, then move:
  zipped_assembled_gp23_fastq=$Sequence_directory/Total_gp23_pear.assembled_lines_annotated_quality_trim.fastq.gz
  assembled_gp23_fastq=$Sequence_directory/Total_gp23_pear.assembled_lines_annotated_quality_trim.fastq
  assembled_gp23_fasta=$Sequence_directory/Total_gp23_pear.assembled_lines_annotated_quality_trim.fasta


  ##### change to zcat for speed!!!

   # unzip fastq
  # gzip -d $zipped_assembled_MPL_fastq

#zcat $zipped_assembled_MPL_fastq > $assembled_MPL_fastq

## too big for usearch


### use andre's alignment!!!
gp23_nuc_ref=../../JerichoAndSOGsequencing/data/Alignments/g23-1400namefixalignRefineAgain_edited_GOS_removed_back_translated.fasta


sed -n '1~4s/^@/>/p;2~4p' $assembled_gp23_fastq> $assembled_gp23_fasta

files=$assembled_gp23_fasta
basename $files
short_name=`basename -s .fasta $files`

$mothur1_3 "#count.seqs(name=$files, processors=8); 
align.seqs(fasta=$files, reference=$gp23_nuc_ref, flip=t, processors=8);
summary.seqs(fasta=current)"


#                 Start   End     NBases  Ambigs  Polymer NumSeqs
# Minimum:        0       0       0       0       1       1
# 2.5%-tile:      6       624     3       0       1       174730
# 25%-tile:       6       885     327     0       4       1747299
# Median:         14      988     415     0       6       3494597
# 75%-tile:       14      1034    448     0       6       5241895
# 97.5%-tile:     989     1034    490     0       6       6814463
# Maximum:        1043    1043    556     0       46      6989192
# Mean:   58.4121 936.876 382.666 0       5.18096
# # of Seqs:      6989192


##now
#                 Start   End     NBases  Ambigs  Polymer NumSeqs
# Minimum:        1       1       1       0       1       1
# 2.5%-tile:      1       4       3       0       1       267
# 25%-tile:       4       642     299     0       4       2665
# Median:         4       657     402     0       6       5329
# 75%-tile:       295     657     447     0       6       7993
# 97.5%-tile:     655     657     486     0       6       10391
# Maximum:        657     657     501     0       8       10657
# Mean:   125.742 589.787 329.362 0       4.9041
# of Seqs:      10657


python send_email_python.py "check gp23 alignment"
read -p "Please edit 02_7_3_align_and_filter_with_mothur_gp23.sh to reflect appropriate screening"


$mothur1_3 "#screen.seqs(fasta=$Sequence_directory/$short_name".align", start=10, end=650, minlength=300, maxhomop=8, processors=8,maxambig=0);
 summary.seqs(fasta=$Sequence_directory/$short_name".good.align");
  filter.seqs(fasta=$Sequence_directory/$short_name".good.align", vertical=T);
  summary.seqs(fasta=current) "

#aliview $Sequence_directory/$short_name".good.only_1000"$direction".fasta"

$mothur1_3 "#degap.seqs(fasta=$Sequence_directory/$short_name".good.align")"

## gives ../results/Divided_by_primers/All_libs_combined/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng.fasta

set +e +u
deactivate