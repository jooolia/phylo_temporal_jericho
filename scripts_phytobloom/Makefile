

Jericho_Pier_original_data = ../../JerichoTimeSeries/data/JerichoDataCollected.csv

Jericho_Pier_selected_data = ../results/Jericho_data_for_env_analysis.csv

gp23_OTU_table = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_annotated_quality_length_filt_otus_95.00global_OTU_table.tsv

16s_OTU_table = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ng_otus_97.00global_OTU_table.tsv

18s_OTU_table = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng_otus_97.00global_OTU_table.tsv


original_gp23_table = ../data/OTU_table_bloom_time_series_gp23.tsv

normalized_gp23_table = ../data/OTU_table_bloom_time_series_normalized_gp23.tsv

## could I make a collection of all of these normalized OTU tables?

normalized_OTU_table_16s = ../data/OTU_table_bloom_time_series_normalized_16S_R1.tsv

original_18s_taxonomy = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.taxonomy

original_16s_taxonomy = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.taxonomy


cleaned_16s_taxonomy = ../results/cleaned_up_16s_taxonomy_Jericho.csv

rarefaction_curves = ../figures/Rarefaction_curves_original_data.pdf

richness_over_time_16s =../figures/time_series_richness_all_amplicons_line_graphs001.pdf

taxonomic_over_time_16s = ../figures/16S_time_series_by_order_specific_OTUs001.pdf

richess_over_time_area = ../figures/time_series_top_10_each_OTUs_areaplots_all_amplicons_all_dates001.pdf

NMDS_16S = ../figures/16S_OTU_table_NMDS.pdf

richness_csv = ../results/amplicon_richness_by_date.csv

richness_csv_by_date = ../results/amplicon_richness_by_date.csv

Jericho_panal_chart = ../figures/Jericho_full_panel_chart.pdf

Jericho_pairs_plot = ../figures/pairs_of_Jericho_parameter.pdf

Frequency_OTUs_16S = ../figures/16S_frequency_of_otus_by_Class_no_unclassified.pdf

community_similarity_over_time = ../figures/time_series_community_similarity_all_amplicons_line_graphs001.pdf

PDF = $(wildcard ../figures/*.pdf)
PNG = $(patsubst ../figures/%.pdf,../figures/%.png,$(PDF))

community_similarity_regressed = ../figures/Community_similarity_matrix_regressed_over_time_one_pager.png

taxonomic_over_time = ../figures/Dinoflagellates_over_bloom.pdf

mpl_pulled_out_for_oligotyping = ../results/derep_MPL_OTU_1588like_reads_aligned.fasta

oligotypes_all_times_figures = ../figures/count_oligo_MPL288_log_y_all_times.pdf

#JAG_theme = ../../JAG_black_presentation.R
JAG_theme = 
## paper

## for manuscript make --dry-run JAG_theme=../../JAG_manuscript_figure.R
## for presentation make --dry-run JAG_theme=../../JAG_black_presentation.R

ifeq ($(JAG_theme),../../JAG_manuscript_figure.R)
        cp ../figures/*.pdf ../figures_presentation
        cp ../figures/*.png ../figures_presentation
endif

CONVERTFLAGS := -density 300 -resize 25% -trim

PICLIST := $(wildcard ../figures/*.pdf)
PDFLIST := $(PICLIST:../figures/%.pdf=../figures/%.png)

#multi_pngs :=  %.o %-1.o %-2.o %-3.o %-4.o %-5.o %-6.o %-7.o %-8.o %-9.o

#png_to_generate = %.o../figures/$(multi_pngs:.o=.png)

## make paper and presentation instead of "all"
all: $(Jericho_Pier_selected_data) $(normalized_OTU_table_16s) $(normalized_gp23_table) $(cleaned_16s_taxonomy)  $(richness_over_time_16s)  $(richness_csv) $(Jericho_panal_chart) $(Dissimilarity_over_time_plots) $(community_similarity_over_time) $(community_similarity_regressed) $(mpl_pulled_out_for_oligotyping) $(taxonomic_over_time_16s) $(taxonomic_over_time) $(oligotypes_all_times_figures) $(richess_over_time_area) $(PDFLIST)
#$(rarefaction_curves) $(Jericho_pairs_plot) $(Frequency_OTUs_16S)  print_plots.html  $(NMDS_16S)
## make sure all pdfs are available for Rmarkdown
%.png : %.pdf
	convert $(CONVERTFLAGS) $< $@

## Process environmental data for Jericho libs

$(Jericho_Pier_selected_data): 001_environmental_data_for_Jericho.R $(Jericho_Pier_original_data)
	Rscript 001_environmental_data_for_Jericho.R

### Initial processing of Miseq sequences for all libraries ####
$(normalized_gp23_table): 01_pull_in_OTUs_specific_to_project_and_normalize.R 002_pull_small_table_of_VC_and_pool_number.R $(gp23_OTU_table) $(16s_OTU_table) $(18s_OTU_table) $(Jericho_Pier_selected_data) 
	Rscript 01_pull_in_OTUs_specific_to_project_and_normalize.R
	Rscript 01_2_pull_in_OTUs_from_total_time_series_and_normalize.R
	Rscript 002_pull_small_table_of_VC_and_pool_number.R

$(cleaned_16s_taxonomy): 01_5_taxonomic_info_reformat.R $(original_18s_taxonomy) $(original_16s_taxonomy)
	Rscript 01_5_taxonomic_info_reformat.R

#$(rarefaction_curves) : 02_rarefaction_curves.R $(original_gp23_table) $(normalized_gp23_table) 
#	Rscript 02_rarefaction_curves.R

$(Jericho_panal_chart): 03_panel_chart_environmental_data.R $(Jericho_Pier_selected_data)
	Rscript 03_panel_chart_environmental_data.R $(JAG_theme)

## Making plots 
$(richness_over_time_16s): 04_richness_of_OTUs_over_time.R $(normalized_gp23_table)
	Rscript 04_richness_of_OTUs_over_time.R $(JAG_theme)


$(richess_over_time_area): 05_OTUs_over_time_and_heatmaps.R $(normalized_gp23_table)
	Rscript 05_OTUs_over_time_and_heatmaps.R  $(JAG_theme)


$(taxonomic_over_time_16s): 05_5_taxonomic_plots.R $(normalized_gp23_table) $(cleaned_16s_taxonomy)
	Rscript --verbose 05_5_taxonomic_plots.R $(JAG_theme)

# $(Jericho_pairs_plot) : 06_pairs_analysis_of_environmental_data.R $(Jericho_Pier_selected_data) $(richness_csv_by_date)
# 	Rscript 06_pairs_analysis_of_environmental_data.R

# $(NMDS_16S) : 07_NMDS_of_samples.R $(normalized_gp23_table)
# 	Rscript 07_NMDS_of_samples.R

# $(Frequency_OTUs_16S) : 11_rank_abundance_by_rare_common_or_abundant.R $(normalized_gp23_table)
# 	Rscript 11_rank_abundance_by_rare_common_or_abundant.R

# # prints out the corragrams. Maybe I will have it just do that and rename the other one? the temporal one gives the same 

$(community_similarity_regressed): 14_community_similarity_matrices_regressed_over_time.R  $(normalized_gp23_table) $(Jericho_Pier_selected_data)
	Rscript 14_community_similarity_matrices_regressed_over_time.R $(JAG_theme)

$(community_similarity_over_time): 15_diversity_over_time_and_shared_between_times.R  $(normalized_gp23_table) $(Jericho_Pier_selected_data)
	Rscript 15_diversity_over_time_and_shared_between_times.R $(JAG_theme)

$(taxonomic_over_time): 53_viral_clades_subset_with_18s.R $(normalized_gp23_table) $(Jericho_Pier_selected_data) $(taxonomic_over_time_16s)
	Rscript 53_viral_clades_subset_with_18s.R 

$(mpl_pulled_out_for_oligotyping): 69_oligotyping_of_specific_MPL_otus.sh $(gp23_OTU_table) $(normalized_gp23_table)
#	bash 69_oligotyping_of_specific_MPL_otus.sh

: 70_align_chosen_MPL_OTUs.sh
: 71_entropy_analysis_of_MPL.sh
: 72_oligotyping_of_MPL_select_OTUs.sh
: 73_pull_out_group_A_nucleotide_and_oligotype.sh
: 74_pull_out_rqaphidophte_nucleotide_and_oligotype.sh
: 75_align_MPL_compare_persistent_vs_ephemeral.sh
: 75_pull_out_cyanos_and_do_oligotyping.sh
: 76_pull_out_gp23_and_oligotype.sh	

$(oligotypes_all_times_figures): 66_plot_oligotypes_over_time.R $(mpl_pulled_out_for_oligotyping)
	Rscript 66_plot_oligotypes_over_time.R 

clean_figs:
	rm ../figures/*.pdf ../figures/*.png