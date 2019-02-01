#!/bin/bash
## Author: 
## Modified by: Julia Gustavsen
set -o nounset -o pipefail -o errexit


usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32

## Create blast database.

### GP23

cat ../data/gp23_sequence_references/myoviridae_complete_genomes_nuc.fasta  ../results/gp23_environmental_refs_nuc_derep.fasta  > ../data/gp23_sequence_references/seqs_for_gp23_blast_db.fasta

makeblastdb -in ../data/gp23_sequence_references/seqs_for_gp23_blast_db.fasta -out ../data/gp23_sequence_references/gp23_nuc -dbtype nucl -hash_index

#### MPL
cat ../data/MPL_sequence_references/picornavirales_RdRp_from_genomes_nuc.fasta  ../results/MPL_environmental_refs_nuc_derep.fasta > ../data/MPL_sequence_references/seqs_for_MPL_blast_db.fasta

makeblastdb -in ../data/MPL_sequence_references/seqs_for_MPL_blast_db.fasta -out ../data/MPL_sequence_references/MPL_nuc -dbtype nucl -hash_index