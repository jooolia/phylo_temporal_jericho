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

#clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64
clustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64

#hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 
hmm_build=hmmbuild


## NCBI protein alignment
NCBI_RdRP_ref=../data/MPL_sequence_references/NCBI_cd01699_RNA_dep_RNAP.fasta 
NCBI_RdRP_ref_renamed=../data/MPL_sequence_references/NCBI_cd01699_RNA_dep_RNAP_renamed.fasta


RdRp_alignment=$Translated_seqs_dir/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.fasta
Sequences_MPL=$Translated_seqs_dir/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 


$hmm_build --informat afa RdRp.hmm $RdRp_alignment



echo "================================================"
echo "Split up file"
echo "================================================"			
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $Sequences_MPL > $Sequences_MPL"_1L.fasta"

split -l 20000 -d $Sequences_MPL"_1L.fasta" $Sequences_MPL"_1L.fasta"
# gives $Translated_seqs_dir/$basefilename"_clustal.fasta"[0-9][0-9]

for split_file in $Sequences_MPL"_1L.fasta"[0-9][0-9]
do

echo $split_file 


$clustal_omega -i $split_file  --hmm-in=RdRp.hmm --outfmt=fasta -v --force -o $split_file"_clustal_aligned.fasta" 

done




## need to do profile alignment for 00 and up until there are no more...

mv $Sequences_MPL"_1L.fasta00_clustal_aligned.fasta" $Sequences_MPL"initial_profile.fasta"
mv $Sequences_MPL"_1L.fasta01_clustal_aligned.fasta" $Sequences_MPL"initial_profile2.fasta"

$clustal_omega --p1=$Sequences_MPL"initial_profile.fasta"  --p2=$Sequences_MPL"initial_profile2.fasta" --outfmt=fasta -v --force -o $Sequences_MPL"_clustal_aligned_all_profile.fasta"


for files in  $Sequences_MPL"_1L.fasta"[0-9][0-9]"_clustal_aligned.fasta"
do 
echo $files

$clustal_omega --p1=$files --p2=$Sequences_MPL"_clustal_aligned_all_profile.fasta" --outfmt=fasta -v --force -o $Sequences_MPL"_clustal_aligned_all_profile.fasta"
done

