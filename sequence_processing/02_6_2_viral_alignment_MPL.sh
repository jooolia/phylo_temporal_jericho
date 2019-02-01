

source /home/labop/Data/Julia/Julia_python/bin/activate

### Nucleotide alignments for viruses
set -o nounset -o pipefail -o errexit

## Programs
Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated
trimal=/home/labop/trimal/source/trimal

clustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64
#clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64

#hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 
hmm_build=hmmbuild

usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32
#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32

mothur1_3=/home/labop/mothur/mothur
#mothur1_3=/home/julia/mothur/mothur

## Reference data lives here
MPL_ref=../data/MPL_sequence_references


## Pull down nucleotide reference material for MPL

## Use list of IDs to get genbank files of the Culley papers
##not sure how to fetch from search and then get sequences from NCBI. That would be cool. 
## this downloads the fastas to a directory
#python fetch_genomes.py $MPL_ref/Cully_environmental_seqeunce_data_gi_numbers.txt $MPL_ref/env_refs/

echo "=========================================="
echo " Fetch nucleotide references MPL"
echo "=========================================="

### MPL refs
for file in $MPL_ref/*.gb
do
base_file_name=$(basename -s ".gb" $file)
python parsing_genbank_for_MPL_nuc.py $file $MPL_ref/$base_file_name"_env_nuc.tsv" 
done

for file in $MPL_ref/*_env_nuc.tsv
do
echo $file
base_file_name=$(basename -s ".tsv" $file)
awk 'BEGIN {FS="\t"} NR <= 1 {next} {print ">"$1"_gi_"$2"\n"$4}' $file > $MPL_ref/$base_file_name".fasta"
sed -i "/^>/s/>/>$base_file_name"_"/g" $MPL_ref/$base_file_name".fasta"
done

## these are protein $MPL_ref/env_refs/*.fasta
#at $MPL_ref/env_refs/*.fasta > ../results/MPL_environmental_refs_nuc.fasta

cat  $MPL_ref/*_env_nuc.fasta > ../results/MPL_environmental_refs_nuc.fasta

## concatenate all of the environmental sequences together
sed -i '/^>/s/ /_/g' ../results/MPL_environmental_refs_nuc.fasta

$usearch8 -derep_fulllength ../results/MPL_environmental_refs_nuc.fasta -fastaout ../results/MPL_environmental_refs_nuc_derep.fasta -sizeout


echo "=========================================="
echo " Align nucleotide references"
echo "=========================================="

## Align the environmental seqs and use those...

mafft --localpair --maxiterate 1000 --thread 8 --op 3 --ep 1.5 --legacygappenalty  --anysymbol ../results/MPL_environmental_refs_nuc_derep.fasta > ../results/MPL_environmental_refs_nuc_derep_MAFFT.fasta

python send_email_python.py "../results/MPL_environmental_refs_nuc_derep_MAFFT.fasta"
read -p "Please edit ../results/MPL_environmental_refs_nuc_derep_MAFFT.fasta to become ../results/MPL_environmental_refs_nuc_derep_MAFFT_edited.fasta"
#aliview ../results/MPL_environmental_refs_nuc_derep_MAFFT.fasta

## downloads complete genomes from NCBI of picornaviridae
## gives file with today's date
python download_picornaviridae_genomes_from_NCBI.py 

## this date thing is maybe problematic
python parsing_genbank_for_picorna_genomes_nuc.py ../data/MPL_sequence_references/picornavirales_complete_genomes_downloaded.fasta  ../data/MPL_sequence_references/picornavirales_RdRp_from_genomes_nuc.tsv

awk 'BEGIN {FS="\t"} NR <= 1 {next} {print ">"$1"_gi_"$2"\n"$4}' ../data/MPL_sequence_references/picornavirales_RdRp_from_genomes_nuc.tsv > $MPL_ref/picornavirales_RdRp_from_genomes_nuc.fasta

mafft --localpair --maxiterate 1000 --thread 8 --op 3 --ep 1.5 --legacygappenalty $MPL_ref/picornavirales_RdRp_from_genomes_nuc.fasta > ../results/picornavirales_RdRp_from_genomes_nuc_MAFFT.fasta

python send_email_python.py "../results/picornavirales_RdRp_from_genomes_nuc_MAFFT.fasta"
read -p "Please edit ../results/picornavirales_RdRp_from_genomes_nuc_MAFFT.fasta to become ../results/picornavirales_RdRp_from_genomes_nuc_MAFFT_edited.fasta"
#aliview ../results/picornavirales_RdRp_from_genomes_nuc_MAFFT.fasta


set +u +e
deactivate