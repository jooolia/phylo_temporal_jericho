source /home/labop/Data/Julia/Julia_python/bin/activate

set -o nounset -o pipefail -o errexit
## Translation of the chimera checked sequences

## Do not need for 18s and 16s sequences

seqtk=/home/labop/Data/software/seqtk/seqtk
#seqtk=/home/julia/seqtk/seqtk

Sequence_directory=../results/Divided_by_primers/All_libs_combined

# gp23_reads=$Sequence_directory/Total_gp23_pear.assembled_lines_annotated_quality_trim.fasta
# # basename .fasta gives :Total_gp23_pear.assembled_lines_annotated_quality_trim
# # will give Total_gp23_pear.assembled_lines_annotated_quality_trim.faa in the translated seqs dir
# MPL_reads=$Sequence_directory/Total_MPL_pear.assembled_lines_annotated_quality_trim.fasta

## did another filtering step with mothur
gp23_reads=$Sequence_directory/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.fasta
# basename .fasta gives :Total_gp23_pear.assembled_lines_annotated_quality_trim
# will give Total_gp23_pear.assembled_lines_annotated_quality_trim.faa in the translated seqs dir

MPL_reads=$Sequence_directory/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.fasta


#FragGeneScan=/home/julia/FragGeneScan1.18/
FragGeneScan=/Data/software/FragGeneScan1.16

Translated_seqs_dir=$Sequence_directory/Translated

mkdir -p $Translated_seqs_dir



for files in $MPL_reads $gp23_reads 
do 
basename_reads=$(basename $files .fasta)
echo $basename_reads
$FragGeneScan/run_FragGeneScan.pl -genome=$files -out $Translated_seqs_dir/$basename_reads -complete=0 -train=illumina_10
echo "Python for getting seqs without stop codons"
python keep_seqs_without_stop_codons.py $Translated_seqs_dir/$basename_reads".faa" keep_no_stop.txt
echo "Seqtk for parsing from original file"
$seqtk subseq $Translated_seqs_dir/$basename_reads".faa" keep_no_stop.txt > test.faa
mv test.faa $Translated_seqs_dir/$basename_reads".faa"

done

set +u +e
deactivate