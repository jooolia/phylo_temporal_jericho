#!/usr/bin/make -f

## matched target is the 'at' symbol
## $< is rule's first prereq
# $? is out of date prereqs
## $^ all prereqs

data_dir = ../data
concat_libs_dir = ../results/Divided_by_primers/All_libs_combined


Trimmed_paired_16S =../results/Divided_by_primers/All_libs_combined/Total_16s_R1_annotated_trimmed_paired.fastq
Alignment_16S_summary = ../results/Divided_by_primers/All_libs_combined/Total_16s_concat_annotated.good.summary
Trimmed_paired_16S_R2 = ../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired.fastq

R2_16s_fasta = ../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired.fasta
Alignment_16S_R2_summary = ../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired.good.summary
#No_gaps_16s_R2= ../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired.good.good.ng.fasta

R2_reverse_16s_aligned = ../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired_reversegood.good.align
R2_reverse_16s_pcr_seqs_no_gaps = ../results/Divided_by_primers/All_libs_combined/Total_16s_R2_annotated_trimmed_paired_reversegood.good.pcr.ng.fasta

concat_16s = ../results/Divided_by_primers/All_libs_combined/Total_16s_concat_annotated.fasta
No_gaps_16S = ../results/Divided_by_primers/All_libs_combined/Total_16s_concat_annotated.good.good.ng.fasta

Sequences_18s = ../results/Divided_by_primers/All_libs_combined/Total_18s_pear.assembled_lines_annotated_quality_trim.fastq
No_gaps_18s = ../results/Divided_by_primers/All_libs_combined/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng.fasta

Total_MPL_annotated_trimmed_paired = ../results/Divided_by_primers/All_libs_combined/Total_MPL_pear.assembled_lines_annotated_quality_trim.fastq
MPL_fasta = ../results/Divided_by_primers/All_libs_combined/Total_MPL_pear.assembled_lines_annotated_quality_trim.fasta

### aligned and filtered

MPL_enviro_refs_align = ../results/MPL_environmental_refs_nuc_derep_MAFFT_edited.fasta

gp23_enviro_refs_align = ../results/gp23_environmental_refs_nuc_derep_MAFFT_edited.fasta

MPL_aligned_and_no_gaps = ../results/Divided_by_primers/All_libs_combined/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng.fasta
gp23_aligned_and_no_gaps = ../results/Divided_by_primers/All_libs_combined/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng.fasta
Total_gp23_annotated_trimmed_paired = ../results/Divided_by_primers/All_libs_combined/Total_gp23_pear.assembled_lines_annotated_quality_trim.fastq

### Screen for contaminants
MPL_blast_database = ../data/MPL_sequence_references/MPL_nuc.nhd
gp23_no_contam = ../results/Divided_by_primers/All_libs_combined/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.fasta

## Translated
gp23_translated = ../results/Divided_by_primers/All_libs_combined/Translated/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa
## Data for making tables:

MPL_translated_filtered = ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa
MPL_translated_filtered_aligned = ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profile.fasta

MPL_translated_filtered_aligned_trimmed = ../results/Divided_by_primers/All_libs_combined/Translated/Total_MPL_pear.assembled_lines_annotated_quality_length_filt.fasta

gp23_translated_filtered = ../results/Divided_by_primers/All_libs_combined/Translated/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 

gp23_translated_filtered_aligned = ../results/Divided_by_primers/All_libs_combined/Translated/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa_clustal_aligned_all_profile.fasta

gp23_translated_filtered_aligned_trimmed = ../results/Divided_by_primers/All_libs_combined/Translated/Total_gp23_pear.assembled_annotated_quality_length_filt.fasta


### OTU tables 

S18_OTU_table = ../results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.ng_otus_97.00global_OTU_table.tsv

MPL_OTU_table = ../results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95global_OTU_table.tsv


### Protein references:
NCBI_gp23 = ../data/gp23_sequence_references/NCBI_ref_pfam07068_gp23_renamed.fasta
NCBI_RdRp = ../data/MPL_sequence_references/NCBI_cd01699_RNA_dep_RNAP_renamed.fasta

### Protein alignments
RdRp_env_and_ref_alignment = ../results/Divided_by_primers/All_libs_combined/Translated/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.fasta


gp23_env_and_ref_alignment = ../results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.fasta



all: $(Trimmed_paired_16S) $(R2_16s_fasta) $(Alignment_16S_R2_summary) $(R2_reverse_16s_aligned) $(R2_reverse_16s_pcr_seqs_no_gaps) $(concat_16s) $(Alignment_16S_summary) $(No_gaps_16S) $(No_gaps_18s)  $(MPL_aligned_and_no_gaps) $(gp23_aligned_and_no_gaps) $(MPL_blast_database) $(gp23_no_contam) $(gp23_translated) $(table_counts) $(RdRp_env_and_ref_alignment) $(MPL_translated_filtered) $(MPL_translated_filtered_aligned) $(MPL_translated_filtered_aligned_trimmed) $(gp23_env_and_ref_alignment) $(gp23_translated_filtered) $(gp23_translated_filtered_aligned) $(gp23_translated_filtered_aligned_trimmed) $(MPL_no_chimeras) $(gp23_no_chimeras) $(MPL_OTU_table) $(S18_OTU_table)  $(table_of_reads_by_pool) $(alluvial_flow_data_table) 


Rscript = /usr/bin/Rscript
### Initial processing of Miseq sequences for all libraries ####

$(Trimmed_paired_16S) : 01_Miseq_run_quality_trimming_and_splitting_by_primer.sh renaming_library_files_according_to_mapping_file.py Counting_primers_with_python.py count_sequence_length.py split_trimmed_libs_by_primer.py count_merged_primer_files.py annotate_lines_of_fastq_libs.py FastQC_plotting_all_libs_together.Rmd Miseq_data_report.rmd Quality_of_primers_by_libs.Rmd
	bash 01_Miseq_run_quality_trimming_and_splitting_by_primer.sh

## fake target
$(Sequences_18s) : $(Trimmed_paired_16S)
$(Total_gp23_annotated_trimmed_paired) : $(Trimmed_paired_16S)
$(Total_MPL_annotated_trimmed_paired) : $(Trimmed_paired_16S)


$(R2_16s_fasta) : 02_1_convert_16S.sh $(Trimmed_paired_16S) 
	bash 02_1_convert_16S.sh
### Align 16s sequences to reference ###

$(Alignment_16S_R2_summary) : 02_2_alignment_16S.sh $(R2_16s_fasta) 
	bash 02_2_alignment_16S.sh
### Screen 16s sequences ###

$(R2_reverse_16s_aligned) : 02_3_1_screen_for_direction_16s_R1_and_R2.sh
	bash 02_3_1_screen_for_direction_16s_R1_and_R2.sh

$(R2_reverse_16s_pcr_seqs_no_gaps) : 02_3_2_pcr-seqs_forward_and_reverse_16s.sh $(R2_reverse_16s_aligned)
	bash 02_3_2_pcr-seqs_forward_and_reverse_16s.sh


$(concat_16s) : 02_4_concat_16S.sh concatenate_R1_and_R2_for_non_merging_primers_already_aligned_forward.py concatenate_R1_and_R2_for_non_merging_primers_already_aligned_reverse.py $(R2_reverse_16s_pcr_seqs_no_gaps)
	bash 02_4_concat_16S.sh

 $(No_gaps_16S) : 02_5_screen_align_and_degap_concat_16S.sh $(concat_16s) 
	bash 02_5_screen_align_and_degap_concat_16S.sh

### Algin and screen 18S ###

$(No_gaps_18s) : 02_alignment_screening_and_degap_18S.sh $(Sequences_18s) $(Trimmed_paired_16S)
	bash 02_alignment_screening_and_degap_18S.sh

#### Aligning and screening viral sequences

$(MPL_enviro_refs_align) : 02_6_2_viral_alignment_MPL.sh $(wildcard ../data/MPL_sequence_references/*.gb)
	bash 02_6_2_viral_alignment_MPL.sh

$(gp23_enviro_refs_align) : 02_6_3_viral_alignment_gp23.sh $(wildcard ../data/gp23_sequence_references/*.gb)
	bash 02_6_3_viral_alignment_gp23.sh


$(MPL_aligned_and_no_gaps) : 02_7_2_align_and_filter_with_mothur_MPL.sh $(Total_MPL_annotated_trimmed_paired) $(Trimmed_paired_16S) $(MPL_enviro_refs_align)
	bash 02_7_2_align_and_filter_with_mothur_MPL.sh

$(gp23_aligned_and_no_gaps) : 02_7_3_align_and_filter_with_mothur_gp23.sh $(Total_gp23_annotated_trimmed_paired) $(Trimmed_paired_16S) $(gp23_enviro_refs_align)
	bash 02_7_3_align_and_filter_with_mothur_gp23.sh

$(MPL_blast_database) : 03_1_make_blast_database_for_viral_seqs.sh $(gp23_enviro_refs_align) $(MPL_enviro_refs_align) $(AVS_env_alignment)
	bash 03_1_make_blast_database_for_viral_seqs.sh

$(gp23_no_contam) : 03_2_blast_to_remove_contaminants.sh $(AVS1_aligned_and_no_gaps) $(MPL_aligned_and_no_gaps) $(gp23_aligned_and_no_gaps)
	bash 03_2_blast_to_remove_contaminants.sh

#### Translation of gp23, and MPL ######

$(gp23_translated) : 04_translation.sh $(gp23_no_contam)
	bash 04_translation.sh

$(AVS_R2_length_filt_clustal_cleaned) : 05_1_length_filter.sh 05_2_align_R1_and_R2_for_template.sh $(gp23_translated)
	bash 05_1_length_filter.sh

$(MPL_translated_filtered) : $(AVS_R2_length_filt_clustal_cleaned)
$(gp23_translated_filtered) : $(AVS_R2_length_filt_clustal_cleaned)


$(MPL_translated_filtered_aligned) : 07_1_alignment_MPL.sh $(RdRp_env_and_ref_alignment) $(MPL_translated_filtered) 
	bash 07_1_alignment_MPL.sh

$(MPL_translated_filtered_aligned_trimmed) : 07_1_2_edit_and_trim_MPL_alignment.sh $(MPL_translated_filtered_aligned)
	bash 07_1_2_edit_and_trim_MPL_alignment.sh

$(gp23_translated_filtered_aligned) : 07_2_gp23_alignment.sh $(gp23_env_and_ref_alignment) $(gp23_translated_filtered)
	bash 07_2_gp23_alignment.sh

$(gp23_translated_filtered_aligned_trimmed) : 07_2_1_edit_and_trim_gp23_alignment.sh $(gp23_translated_filtered_aligned)
	bash 07_2_1_edit_and_trim_gp23_alignment.sh



#### OTU picking of all amplicons. Also taxonomic classification of 16s and 18s OTUs ####
$(MPL_OTU_table) : 08_OTU_clustering_various_methods.sh $(gp23_translated_filtered_aligned_trimmed) $(MPL_translated_filtered_aligned_trimmed) get_intersection_of_non_chimeric_sequences.py 
	bash 08_OTU_clustering_various_methods.sh

$(S18_OTU_table) : 08_1_OTU_clustering_16s_and_18s.sh $(No_gaps_16S) $(No_gaps_18s) get_intersection_of_non_chimeric_sequences.py 
	bash 08_1_OTU_clustering_16s_and_18s.sh

# clean:
# 	rm testfile.txt
# 	echo 'testfile.txt has been deleted'
