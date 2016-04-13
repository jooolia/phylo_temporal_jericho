#!/bin/bash

## Author=Julia Gustavsen
## Last edited: 2 September 2014

## test out making trees for RdRP for Jericho

## Use a collection of the marine viruses and run them with the Mafft homolog 
## webserver (it is also a program too, but doesn't seem well supported)
## Then I make sure to identify the sequences so that they can be removed later
## I add the NCBI taxonomy and then trim the alignment to the desired length.
## I need to keep looking at it using aliview so that I don't make mistakes. 
## I cat together previous environmental sequences and then cluster them at 95% identity.
## These are then added to the previous alignments. Alignthe 454 sequences and
## see which ones are not alignable. Remove them, then remove gaps and add 454 
## sequences to the existing alignment using the MAFFT --add option. Trim the 
## alignment using trimal and then make a test tree using Fasttree. 
## Send the alignment to RaxML webserver.

source /home/labop/Data/Julia/Julia_python/bin/activate


set -o nounset -o pipefail -o errexit

## Programs
Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated

hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 

usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32

Database_dir=../../Taxonomic_databases

usearch7=/home/labop/usearch7/usearch7.0.1090_i86linux32
#usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32

## testing out the new vsearch for things requiring large amounts of memory
vsearch=/home/labop/vsearch/vsearch-1.0.6-linux-x86_64


#mothur1_3=/home/labop/mothur/mothur
mothur1_3=/home/julia/mothur/mothur




S18_miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.fasta
S18_miseq_amplicons_sed=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.ng_otus_97.00_sed.fasta
normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_18S.tsv
normalized_otu_table_top_100=../data/OTU_table_Jericho_time_series_18s_normalized_top_100.tsv
Times_series_18s=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.fasta
Times_series_18s_top_100=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series_top_100.fasta


sed '/^>/s/;/\./g;/^>/s/=/\./g ' $S18_miseq_amplicons > $S18_miseq_amplicons_sed 

 python filter_fasta_based_on_normalzed_data.py $S18_miseq_amplicons_sed $normalized_otu_table $Times_series_18s

 python filter_fasta_based_on_normalzed_data.py $S18_miseq_amplicons_sed $normalized_otu_table_top_100 $Times_series_18s_top_100

# ##re-running some parts of this
 $mothur1_3 "#count.seqs(name=$Times_series_18s); align.seqs(fasta=$Times_series_18s, reference=$Database_dir/silva.seed_v119_cut_for_18s, flip=t, processors=8); summary.seqs(fasta=current)"

 $mothur1_3 "#filter.seqs(fasta=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.align, trump=., vertical=T)"

 aliview ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter.fasta
# aliview ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.fasta

## I think this gets rid of a lot of the diversity. but could get rid of the gaps. 
 #$mothur1_3 "#filter.seqs(fasta=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.fasta, soft=50, vertical=T)"
$mothur1_3 "#filter.seqs(fasta=../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.fasta, trump=- vertical=T)"
# aliview ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.fasta

sed -i 's/\./_/g;s/\+//g' ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.fasta
## check alignment by using FastTree for estimation
/home/julia/FastTree -nt ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.fasta  > ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree

java -jar /home/julia/FigTree/FigTree_v1.4.0/lib/figtree.jar ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree
