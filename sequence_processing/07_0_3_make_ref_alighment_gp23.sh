## migth be useful to work on sequence alignment with all of the data? 
## FOR PROTEINS
## set up some directories for keeping downloaded references

## ../data/gp23_sequence_references etc
### Get references from Genbank
### Nucleotide alignments for viruses
source /home/labop/Data/Julia/Julia_python/bin/activate


set -o nounset -o pipefail -o errexit

## Programs
Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated
clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64

hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 

usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32

## NCBI protein alignment
NCBI_gp23_ref=../data/gp23_sequence_references/NCBI_ref_pfam07068_gp23.fasta
NCBI_gp23_ref_renamed=../data/gp23_sequence_references/NCBI_ref_pfam07068_gp23_renamed.fasta
### Do alignments ###
Andre_ref=../../JerichoAndSOGsequencing/data/Alignments/g23-1400namefixalignRefineAgain_edited.fas

gp23_ref=../data/gp23_sequence_references

## put the aligned and manipulated data in ../results/alignments/

align_dir=../results/alignments

mkdir -p $align_dir


echo "=========================================="
echo " Parse protein references"
echo "=========================================="


for file in $gp23_ref/*.gb
do
base_file_name=$(basename -s ".gb" $file)
python parsing_genbank_for_gp23_prot.py $file $gp23_ref/$base_file_name"_prot.tsv" 
done

echo "=========================================="
echo " Make sure names are human readable"
echo "=========================================="

for file in $gp23_ref/*_prot.tsv
do
echo $file
base_file_name=$(basename -s ".tsv" $file)
awk 'BEGIN {FS="\t"} NR <= 1 {next} {print ">"$1"_gi_"$2"\n"$4}' $file > $gp23_ref/$base_file_name".fasta"
sed -i "/^>/s/>/>$base_file_name"_"/g" $gp23_ref/$base_file_name".fasta"
done

cat $gp23_ref/*[0-9]_prot.fasta > ../results/gp23_environmental_refs_prot.fasta

sed -i '/^>/s/ /_/g' ../results/gp23_environmental_refs_prot.fasta

$usearch8 -derep_fulllength ../results/gp23_environmental_refs_prot.fasta -fastaout ../results/gp23_environmental_refs_prot_derep.fasta -sizeout

echo "=========================================="
echo " Rename isolates and references"
echo "=========================================="

python parsing_genbank_for_myoviridae_genomes_prot.py ../data/gp23_sequence_references/myoviridae_complete_genomes_downloaded.fasta  ../data/gp23_sequence_references/myoviridae_complete_genomes_prot.tsv

echo "Making more readable "
awk 'BEGIN {FS="\t"} NR <= 1 {next} {print ">"$1"_gi_"$2"\n"$4}' ../data/gp23_sequence_references/myoviridae_complete_genomes_prot.tsv > ../data/gp23_sequence_references/myoviridae_complete_genomes_prot.fasta
sed -i '/^>/s/ /_/g' ../data/gp23_sequence_references/myoviridae_complete_genomes_prot.fasta



$hmm_build --informat afa $Translated_seqs_dir/gp23_NCBI_ref.hmm $Andre_ref


echo "=========================================="
echo " Align protein references"
echo "=========================================="


## align environmental

$clustal_omega -i ../results/gp23_environmental_refs_prot_derep.fasta --hmm-in=$Translated_seqs_dir/gp23_NCBI_ref.hmm  -auto -o $Translated_seqs_dir/environmental_gp23_seqs_aligned_clustal.fasta

aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_clustal.fasta
#aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_clustal_edited.fasta 


## align all the gp23 I could find from genomes
$clustal_omega -i ../data/gp23_sequence_references/myoviridae_complete_genomes_prot.fasta --hmm-in=$Translated_seqs_dir/gp23_NCBI_ref.hmm  -auto -o $Translated_seqs_dir/gp23_myoviridae_genomes_aligned_clustal.fasta

aliview $Translated_seqs_dir/gp23_myoviridae_genomes_aligned_clustal.fasta
#aliview $Translated_seqs_dir/gp23_myoviridae_genomes_aligned_clustal_edited.fasta


## add the ncbi references to the environmental sequences
$clustal_omega --p1 $Translated_seqs_dir/environmental_gp23_seqs_aligned_clustal_edited.fasta --p2 $NCBI_gp23_ref_renamed  -auto -o $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_clustal.fasta

aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_clustal.fasta
#aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_clustal_edited.fasta

$clustal_omega --p1 $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_clustal_edited.fasta --p2 $Translated_seqs_dir/gp23_myoviridae_genomes_aligned_clustal_edited.fasta  -auto -o $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal.fasta

aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal.fasta
#aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited.fasta


$clustal_omega -i $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited.fasta --hmm-in=$Translated_seqs_dir/gp23_NCBI_ref.hmm  -auto -o $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned.fasta

aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned.fasta
#aliview $Translated_seqs_dir/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.fasta

## check using fasttree


## check alignment by using FastTree for estimation
sed -i 's/\;/_/g;s/\+//g' ../results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.fasta

/home/julia/FastTree ../results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.fasta  > ../results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.tree



java -jar /home/julia/FigTree/FigTree_v1.4.0/lib/figtree.jar ../results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.tree
