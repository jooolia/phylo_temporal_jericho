## migth be useful to work on sequence alignment with all of the data? 
## FOR PROTEINS
## set up some directories for keeping downloaded references

## 
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
NCBI_RdRP_ref=../data/MPL_sequence_references/NCBI_cd01699_RNA_dep_RNAP.fasta 
NCBI_RdRP_ref_renamed=../data/MPL_sequence_references/NCBI_cd01699_RNA_dep_RNAP_renamed.fasta

## need to download isolates
## ok downloaded RdRp of known isolates

## as isolate_name_gi_ginumber.fasta
## add them all together:
cat ../data/MPL_sequence_references/*_gi_*.fasta > ../data/MPL_sequence_references/marine_isolates_RdRp.fasta

## trying to rename these to something better

grep -c ">" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta

## search for the gi numbers

grep -E ">gi" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta

echo "Next step: changing the names."

grep -E ">gi\|[[:digit:]]+\|\w+\|(\w+\.[0-9])\| .* \[(.*)\]" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta
sed -i -E  "s/^>gi\|[[:digit:]]+\|\w+\|(\w+\.[0-9])\| .* \[(.*)\]/>\2_gi_\1/" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta


grep -E ">gi\|([[:digit:]]+)\:.*\[(.*)\]" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta
sed -i -E "s/^>gi\|([[:digit:]]+)\:.* \[(.*)\]/>\2_gi_\1/" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta 

grep -E ">gb\|(.+\.[0-9])\|\:.* \[(.*)\]" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta
sed -i -E "s/^>gb\|(.+\.[0-9])\|\:.* \[(.*)\]/>\2_gb_\1/" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta


sed -i -E "s/ /_/g" ../data/MPL_sequence_references/marine_isolates_RdRp.fasta

## but still need to rename all the reference seqquences. 

### Do alignments ###

## put the aligned and manipulated data in ../results/alignments/

align_dir=../results/alignments

mkdir -p $align_dir

MPL_ref=../data/MPL_sequence_references

echo "=========================================="
echo " Parse protein references"
echo "=========================================="


for file in $MPL_ref/*.gb
do
base_file_name=$(basename -s ".gb" $file)
python parsing_genbank_for_MPL_prot.py $file $MPL_ref/$base_file_name"_prot.tsv" 
done

echo "=========================================="
echo " Make sure names are human readable"
echo "=========================================="

for file in $MPL_ref/*_prot.tsv
do
echo $file
base_file_name=$(basename -s ".tsv" $file)
## rename header so that it is something more pleaseing 
awk 'BEGIN {FS="\t"} NR <= 1 {next} {print ">"$1"_gi_"$2"\n"$4}' $file > $MPL_ref/$base_file_name".fasta"
sed -i "/^>/s/>/>$base_file_name"_"/g" $MPL_ref/$base_file_name".fasta"
done

cat $MPL_ref/*_prot.fasta > ../results/MPL_environmental_refs_prot.fasta

## remove spaces
sed -i '/^>/s/ /_/g' ../results/MPL_environmental_refs_prot.fasta

$usearch8 -derep_fulllength ../results/MPL_environmental_refs_prot.fasta -fastaout ../results/MPL_environmental_refs_prot_derep.fasta -sizeout

## remove spaces
sed -i '/^>/s/ /_/g' ../results/MPL_environmental_refs_prot_derep.fasta


echo "=========================================="
echo " Rename isolates and references"
echo "=========================================="


$hmm_build --informat afa ../results/NCBI_RdRp_ref.hmm $NCBI_RdRP_ref


echo "=========================================="
echo " Align protein references"
echo "=========================================="


## align environmental

$clustal_omega -i ../results/MPL_environmental_refs_prot_derep.fasta --hmm-in=../results/NCBI_RdRp_ref.hmm  -auto -o $Translated_seqs_dir/environmental_RdRp_seqs_prot_aligned_clustal.fasta

aliview $Translated_seqs_dir/environmental_RdRp_seqs_prot_aligned_clustal.fasta
#aliview $Translated_seqs_dir/environmental_RdRp_seqs_prot_aligned_clustal_edited.fasta



$clustal_omega -i ../data/MPL_sequence_references/marine_isolates_RdRp.fasta --hmm-in=../results/NCBI_RdRp_ref.hmm -auto -o $Translated_seqs_dir/marine_isolates_RdRp_aligned_clustal.fasta

aliview $Translated_seqs_dir/marine_isolates_RdRp_aligned_clustal.fasta
#aliview $Translated_seqs_dir/marine_isolates_RdRp_aligned_clustal_edited.fasta




## add the ncbi references. 
$clustal_omega --p1 $Translated_seqs_dir/marine_isolates_RdRp_aligned_clustal_edited.fasta --p2 $NCBI_RdRP_ref_renamed  -auto -o $Translated_seqs_dir/marine_isolates_RdRp_seqs_aligned_with_references_clustal.fasta

aliview $Translated_seqs_dir/marine_isolates_RdRp_seqs_aligned_with_references_clustal.fasta
#aliview $Translated_seqs_dir/marine_isolates_RdRp_seqs_aligned_with_references_clustal_edited.fasta


$clustal_omega -i ../results/MPL_environmental_refs_prot_derep.fasta --p1 $Translated_seqs_dir/marine_isolates_RdRp_seqs_aligned_with_references_clustal_edited.fasta -auto -o $Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal.fasta

aliview $Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal.fasta
#aliview $Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.fasta


sed -i 's/\;/_/g;s/\+//g' $Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.fasta
## check alignment with a tree

## check alignment by using FastTree for estimation
/home/julia/FastTree $Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.fasta > $Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.tree

java -jar /home/julia/FigTree/FigTree_v1.4.0/lib/figtree.jar $Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.tree 

