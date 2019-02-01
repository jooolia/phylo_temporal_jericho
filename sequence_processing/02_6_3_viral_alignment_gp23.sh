
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

gp23_ref=../data/gp23_sequence_references


## Pull down nucleotide reference material for gp23

echo "=========================================="
echo " Fetch nucleotide references"
echo "=========================================="

### gp23 refs

## searching in NCBI:
#txid336724[Organism:noexp] AND ("400"[SLEN] : "700"[SLEN])

## have searched for the references and then download all the genbank files associated with those links. 
for file in $gp23_ref/*.gb
do
base_file_name=$(basename -s ".gb" $file)
python parsing_genbank_for_gp23_nuc.py $file $gp23_ref/$base_file_name"_nuc.tsv" 
done

for file in $gp23_ref/*_nuc.tsv
do
echo $file
base_file_name=$(basename -s ".tsv" $file)
awk 'BEGIN {FS="\t"} NR <= 1 {next} {print ">"$1"_gi_"$2"\n"$4}' $file > $gp23_ref/$base_file_name".fasta"
sed -i "/^>/s/>/>$base_file_name"_"/g" $gp23_ref/$base_file_name".fasta"
done

cat $gp23_ref/*_nuc.fasta > ../results/gp23_environmental_refs_nuc.fasta
## concatenate all of the environmental sequences together

sed -i '/^>/s/ /_/g' ../results/gp23_environmental_refs_nuc.fasta

$usearch8 -derep_fulllength ../results/gp23_environmental_refs_nuc.fasta -fastaout ../results/gp23_environmental_refs_nuc_derep.fasta -sizeout

## Align the environmental seqs and use those...


mafft --localpair --maxiterate 1000 --thread 8 --op 3 --ep 1.5 --legacygappenalty ../results/gp23_environmental_refs_nuc_derep.fasta > ../results/gp23_environmental_refs_nuc_derep_MAFFT.fasta

python send_email_python.py "../results/gp23_environmental_refs_nuc_derep_MAFFT.fasta"
read -p "Please edit ../results/gp23_environmental_refs_nuc_derep_MAFFT.fasta to become ../results/gp23_environmental_refs_nuc_derep_MAFFT_edited.fasta"
#aliview ../results/gp23_environmental_refs_nuc_derep_MAFFT.fasta

echo "=========================================="
echo " Align nucleotide references gp23"
echo "=========================================="

## downloads complete genomes from NCBI of phycondaviridae
python download_myoviridae_genomes_from_NCBI.py 

python parsing_genbank_for_gp23_genomes_nuc.py ../data/gp23_sequence_references/myoviridae_complete_genomes_downloaded.fasta  ../data/gp23_sequence_references/myoviridae_complete_genomes_nuc.tsv

awk 'BEGIN {FS="\t"} NR <= 1 {next} {print ">gi_"$2"\n"$4}' ../data/gp23_sequence_references/myoviridae_complete_genomes_nuc.tsv > $gp23_ref/myoviridae_complete_genomes_nuc.fasta


mafft --localpair --maxiterate 1000 --thread 8 --op 3 --ep 1.5 --legacygappenalty $gp23_ref/myoviridae_complete_genomes_nuc.fasta > ../results/myoviridae_complete_genomes_nuc_MAFFT.fasta

python send_email_python.py "../results/myoviridae_complete_genomes_nuc_MAFFT.fasta"
read -p "Please edit ../results/myoviridae_complete_genomes_nuc_MAFFT.fasta to become ../results/myoviridae_complete_genomes_nuc_MAFFT_edited.fasta"
#aliview ../results/myoviridae_complete_genomes_nuc_MAFFT.fasta


set +u +e
deactivate


