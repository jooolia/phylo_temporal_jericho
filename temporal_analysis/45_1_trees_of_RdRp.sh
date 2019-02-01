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
clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64

hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 

usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32

## NCBI protein alignment
NCBI_RdRP_ref=../data/MPL_sequence_references/NCBI_cd01699_RNA_dep_RNAP.fasta 
NCBI_RdRP_ref_renamed=../data/MPL_sequence_references/NCBI_cd01699_RNA_dep_RNAP_renamed.fasta


usearch7=/home/julia/usearch/usearch7.0.1090_i86linux32
##mafft version

#lustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64
clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64

## trim it so that it is just the RdRp section


env_iso_and_ref_align=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.fasta

#$usearch7 -sortbylength ../results/usearch/95identity/RdRp_OTUs_95identity.fasta -output ../data/Alignments/RdRp_OTUs_95identity_over75aa.fasta -minseqlength 75

RdRp_miseq_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95.00.fasta

normalized_otu_table=../data/OTU_table_Jericho_time_series_normalized_MPL.tsv
normalized_otu_table_top_100=../data/OTU_table_Jericho_time_series_MPL_normalized_top_100.tsv
Times_series_RdRp=../results/Total_MPL_filtered_otus_95.00_for_Jericho_Time_series.fasta
Times_series_RdRp_top_100=../results/Total_MPL_filtered_otus_95.00_for_Jericho_Time_series_top_100.fasta

python filter_fasta_based_on_normalzed_data.py $RdRp_miseq_amplicons $normalized_otu_table $Times_series_RdRp

python filter_fasta_based_on_normalzed_data.py $RdRp_miseq_amplicons $normalized_otu_table_top_100 $Times_series_RdRp_top_100


echo "############################################################"
echo "### Make alignment with all OTUs from normalized data ######"
echo "############################################################"
## need to get only those that are in the normalized data. 

$clustal_omega -i $Times_series_RdRp --p1=$env_iso_and_ref_align --outfmt=fasta -v --force -o ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal.fasta

aliview ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal.fasta
#aliview ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta

sed -i 's/\;/_/g;s/\+//g' ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta

sed -i 's/\ //g' ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta


sed -i 's/=//g' ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta

## check alignment by using FastTree for estimation
/home/julia/FastTree ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta  > ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.tree

### model test beofre sending to RaxML
java -jar /home/julia/prottest-3.4-20140123/prottest-3.4.jar -i ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta

java -jar /home/julia/FigTree/FigTree_v1.4.0/lib/figtree.jar ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.tree 



## then do automated trimming
 /home/julia/trimal/source/trimal -in  ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta -out ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled.fasta -fasta -htmlout ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled.html -automated1

  aliview ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled.fasta

#### Cosmetic tidying before tree-building and visualizing ####
## remove semi-colons
sed 's/\;/_/g;s/\+//g' ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled.fasta  > ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled_no_colons.fasta  


## check alignment by using FastTree for estimation
/home/julia/FastTree ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled_no_colons.fasta > ../results/../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled_no_colons.tree

java -jar /home/julia/FigTree/FigTree_v1.4.0/lib/figtree.jar ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled_no_colons.tree 


### model test beofre sending to RaxML
#java -jar /home/julia/prottest-3.4-20140123/prottest-3.4.jar -i ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimalled_no_colons.tree

