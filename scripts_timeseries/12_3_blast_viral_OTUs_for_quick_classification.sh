#!/bin/bash
## Author:  Author: Xi A Tian (alvinxtian@hotmail.com)
## Modified by: Julia Gustavsen


## directory containing NCBI databases
DB_DIR=/Data/blastdb
OUT_DIR=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined_OTUs_picked/Translated/Picked_OTUs/


gp23=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined_OTUs_picked/Translated/Picked_OTUs/Total_gp23_filtered_otus_95.00.fasta

AVS=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined_OTUs_picked/Translated/Picked_OTUs/AVS_concat_middle_missing_section_removed_otus_90.00.fasta

MPL=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined_OTUs_picked/Translated/Picked_OTUs/Total_MPL_filtered_otus_95.00.fasta



## blast the viral OTUs so that I have a best blast hit to give to quickly taxonomically id the OTUs


echo "blast no hit agaisnt nr"

tblastn -query $gp23 -db nt -outfmt "6 qseqid sseqid sscinames sdescr evalue bitscore" -evalue 1e-5 -out $OUT_DIR/gp23_blast_results.txt -num_threads 10 -max_target_seqs 1
tblastn -query $AVS -db nt -outfmt "6 qseqid sseqid sscinames sdescr evalue bitscore" -evalue 1e-5 -out $OUT_DIR/AVS_blast_results.txt -num_threads 10 -max_target_seqs 1
tblastn -query $MPL -db nt -outfmt "6 qseqid sseqid sscinames sdescr evalue bitscore" -evalue 1e-5 -out $OUT_DIR/MPL_blast_results.txt -num_threads 10 -max_target_seqs 1