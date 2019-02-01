#!/usr/bin/make

Jericho_Pier_original_data = ../data/JerichoDataCollected.csv

Jericho_Pier_selected_data = ../results/Jericho_data_for_env_analysis.csv

gp23_OTU_table = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_annotated_quality_length_filt_otus_95.00global_OTU_table.tsv

16s_OTU_table = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_concat_annotated.good.good.ng_otus_97.00global_OTU_table.tsv

18s_OTU_table = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ng_otus_97.00global_OTU_table.tsv

original_gp23_table = ../data/OTU_table_Jericho_time_series_gp23.tsv

normalized_gp23_table = ../data/OTU_table_Jericho_time_series_normalized_gp23.tsv

## could I make a collection of all of these normalized OTU tables?
normalized_OTU_table_16s = ../data/OTU_table_Jericho_time_series_normalized_16S_R1.tsv

RdRp_miseq_amplicons = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95.00.fasta

gp23_miseq_amplicons = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_annotated_quality_length_filt_otus_95.00.fasta

## alignments :

RdRp_env_iso_and_ref_align = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/environmental_RdRp_marine_isolates_and_ref_NCBI_aligned_clustal_edited.fasta

gp23_env_iso_and_ref_align = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.fasta

original_18s_taxonomy = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_18s_pear.assembled_lines_annotated_quality_trim.good.filter.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.taxonomy

original_16s_taxonomy = ../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_16s_forwardgood.good.pcr.ngnon_chimeras_ref97.00.nr_v119.wang.taxonomy

cleaned_16s_taxonomy = ../results/cleaned_up_16s_taxonomy_Jericho.csv

rarefaction_curves = ../figures/Rarefaction_montage.png

top_20_otus = ../figures/time_series_top_20_OTUs_barplots_all_amplicons_no_high_res001.pdf

richness_over_time_16s = ../figures/16S_proportional_area_time_series001.pdf

richness_cor = ../results/amplicon_richness_by_date.csv

Jericho_panel_chart = ../figures/Jericho_full_panel_chart.pdf

JAG_theme = ../JAG_manuscript_figure.R

RdRp_ref_95_tree = ../results/RAxML_bipartitions.RdRptree

gp23_tree_w_miseq = ../results/gp23_95_miseq_data_with_env_iso_and_ref_align_clustal_edited_trimmed_Filee_1L_RAxML_bipartitions.result

S18_tree_w_miseq = ../results/Total_18s_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree

S16_tree_w_miseq = ../results/Total_16s_R1_filtered_otus_97.00_for_Jericho_Time_series.filter_edited.filter.tree

S16_heatmap_by_tree = ../figures/S16_97_miseq_data_heatmap_no_high_res.pdf

phylogenetic_diversity = ../figures/phylogenetic_diversity_over_time_all_amplicons_line_graphs001.pdf

amplicon_env_cor = ../results/amplicon_richness_by_date.csv

group_A_with_HAKA = ../figures/group_A_with_Raphidophytes.pdf

overall_bar_and_phylo_virus = ../figures/viral_groups_over_time_no_high_res_prop.pdf

gp23_groups_with_16s_orders = ../results/gp23_groups_with_16S_orders_cor.csv

pd_com_shared_over_time_fig = ../figures/PD_Shared_OTUS_comm_sim_over_time_all_amplicons_line_graphs_labelled001.pdf

pers_and_eph_figs = ../figures/S16_miseq_data_heatmap_persistant_and_ephemeral_no_high_res.pdf

light_data = ../data/ibis_shortwave_irradiance_2010_2011.csv

light_data_fixed = ../data/ibis_shortwave_irradiance_2010_2011_date_fixed.csv

mantel_test_csv = ../results/mantel_tests_otus_and_env.csv

zoomed_in_tree = ../figures/gp23_zoomed_in_tree1.pdf

bootstrapped_tree = ../figures/gp23_tree_with_bootstrap_and_points.pdf

RDA_mpl = ../figures/RDA_normalized_MPL_group_sums_rda_001.pdf

bubbleplot_16s = ../figures/S16_97_miseq_data_bubbleplot_no_high_res.pdf

Rscript = /usr/bin/Rscript

all: Makefile.gml $(Jericho_Pier_selected_data) $(light_data_fixed) $(Jericho_panel_chart) $(normalized_gp23_table) $(cleaned_16s_taxonomy) $(rarefaction_curves) $(top_20_otus) $(richness_over_time_16s) $(Frequency_OTUs_16S) $(RDA_mpl) $(mantel_test_csv) $(RdRp_tree_w_miseq) $(gp23_tree_w_miseq) $(S18_tree_w_miseq) $(S16_tree_w_miseq) $(RdRp_ref_95_tree) $(S16_heatmap_by_tree) $(phylogenetic_diversity) $(overall_bar_and_phylo_virus) $(pers_and_eph_figs) $(amplicon_env_cor) $(zoomed_in_tree) $(bootstrapped_tree) $(gp23_groups_with_16s_orders) $(group_A_with_HAKA) $(pd_com_shared_over_time_fig)  

$(Jericho_Pier_selected_data): 01_environmental_data_for_Jericho.R $(Jericho_Pier_original_data)
	$(Rscript) 01_environmental_data_for_Jericho.R

$(light_data_fixed): 02_light_data_for_Jericho.R $(light_data)
	$(Rscript) 02_light_data_for_Jericho.R

# ### Initial processing of Miseq sequences for all libraries ####
$(normalized_gp23_table): 03_pull_in_OTU_specific_to_project_and_normalize.R $(gp23_OTU_table) $(16s_OTU_table) $(18s_OTU_table) $(Jericho_Pier_selected_data)
	$(Rscript) 03_pull_in_OTU_specific_to_project_and_normalize.R

$(16s_OTU_table_top_100): 04_get_top_100_OTUs.R $(normalized_gp23_table)
	$(Rscript) 04_get_top_100_OTUs.R

$(cleaned_16s_taxonomy): 05_taxonomic_info_reformat.R ../../making_colour_palette_for_all_studies.R $(original_18s_taxonomy) $(original_16s_taxonomy) $(normalized_gp23_table)
	$(Rscript) 05_taxonomic_info_reformat.R
	$(Rscript) ../../making_colour_palette_for_all_studies.R

# ## actually depends on the original OTU table, but that is also pulled in when the normalization happens. 
$(rarefaction_curves): 02_rarefaction_curves.R $(normalized_gp23_table) $(JAG_theme)
	$(Rscript) 02_rarefaction_curves.R $(JAG_theme)
	bash 02_2_montage_of_rarefaction_curves.sh

$(Jericho_panel_chart): 07_panel_chart_environmental_data.R $(Jericho_Pier_selected_data) $(JAG_theme)
	$(Rscript) 07_panel_chart_environmental_data.R $(JAG_theme)

$(top_20_otus): 08_OTUs_over_time_and_heatmaps.R $(normalized_gp23_table) $(JAG_theme) $(cleaned_16s_taxonomy) 
	$(Rscript) 08_OTUs_over_time_and_heatmaps.R $(JAG_theme)


$(richness_over_time_16s): 09_taxonomic_plots.R  $(normalized_gp23_table) $(JAG_theme) $(cleaned_16s_taxonomy) 
	$(Rscript) 09_taxonomic_plots.R $(JAG_theme)

$(RDA_mpl): 21_RDAs_of_euks_bac_and_viruses.R $(Jericho_Pier_selected_data) $(normalized_gp23_table) $(cleaned_16s_taxonomy)
	$(Rscript) 21_RDAs_of_euks_bac_and_viruses.R


$(mantel_test_csv): 25_Mantel_test.R
	$(Rscript) 25_Mantel_test.R

$(RdRp_tree_w_miseq): 45_1_trees_of_RdRp.sh $(normalized_gp23_table) $(16s_OTU_table_top_100) $(RdRp_miseq_amplicons) filter_fasta_based_on_normalzed_data.py $(RdRp_env_iso_and_ref_align)
	 bash 45_1_trees_of_RdRp.sh

$(gp23_tree_w_miseq): 45_2_trees_of_gp23.sh $(normalized_gp23_table) $(16s_OTU_table_top_100) $(gp23_miseq_amplicons) filter_fasta_based_on_normalzed_data.py $(gp23_env_iso_and_ref_align)
	 bash 45_2_trees_of_gp23.sh

$(S18_tree_w_miseq): 45_4_trees_of_18s.sh $(normalized_gp23_table) $(16s_OTU_table_top_100) filter_fasta_based_on_normalzed_data.py 
	 bash 45_4_trees_of_18s.sh

$(S16_tree_w_miseq): 45_5_trees_16s.sh $(normalized_gp23_table) $(16s_OTU_table_top_100) filter_fasta_based_on_normalzed_data.py
	 bash 45_5_trees_16s.sh		

$(S16_heatmap_by_tree): 46_heatmaps_ordered_by_trees.R $(RdRp_ref_95_tree) $(gp23_tree_w_miseq)
	$(Rscript) 46_heatmaps_ordered_by_trees.R $(JAG_theme)

$(bubbleplot_16s): 46_2_bubble_plots.R $(normalized_gp23_table) $(gp23_tree_w_miseq)
 	$(Rscript) 46_2_bubble_plots.R

$(phylogenetic_diversity): 46_5_phylogenetic_diversity_over_time.R $(RdRp_ref_95_tree) $(gp23_tree_w_miseq)
	$(Rscript) 46_5_phylogenetic_diversity_over_time.R $(JAG_theme)	
	
$(pers_and_eph_figs): 46_6_heatmaps_persistant_or_ephemeral_ordered_by_trees.R $(RdRp_ref_95_tree) $(gp23_tree_w_miseq)
	$(Rscript) 46_6_heatmaps_persistant_or_ephemeral_ordered_by_trees.R $(JAG_theme)

$(overall_bar_and_phylo_virus): 47_top_20_OTUs_over_time_by_phylogenetic_grouping.R $(gp23_tree_w_miseq)
	$(Rscript) 47_top_20_OTUs_over_time_by_phylogenetic_grouping.R $(JAG_theme)

$(zoomed_in_tree): 47_2_divided_up_trees.R $(gp23_tree_w_miseq)
	$(Rscript) 47_2_divided_up_trees.R  $(JAG_theme)

$(bootstrapped_tree): 47_3_bootstrapped_tree.R  $(gp23_tree_w_miseq)
	$(Rscript) 47_3_bootstrapped_tree.R  $(JAG_theme)

$(amplicon_env_cor): 50_correlations_with_PD_and_env.R $(normalized_gp23_table) $(gp23_tree_w_miseq)
 	$(Rscript) 50_correlations_with_PD_and_env.R

$(gp23_groups_with_16s_orders): 52_viral_clades_correlate_to_env_params.R $(normalized_gp23_table) $(gp23_tree_w_miseq) $(Jericho_Pier_selected_data)
	$(Rscript) 52_viral_clades_correlate_to_env_params.R $(JAG_theme)

$(group_A_with_HAKA): 53_viral_clades_subset_with_18s.R $(normalized_gp23_table) $(cleaned_16s_taxonomy) $(gp23_tree_w_miseq) $(RDA_mpl)
	$(Rscript) 53_viral_clades_subset_with_18s.R $(JAG_theme)

$(pd_com_shared_over_time_fig): 60_combining_all_plots_for_PD_shared_and_comm_sim.R $(normalized_gp23_table) $(gp23_tree_w_miseq)
	$(Rscript) 60_combining_all_plots_for_PD_shared_and_comm_sim.R $(JAG_theme)

# Render a figure of this pipeline
Makefile.gml: Makefile
	makefile2graph > Makefile.gv
	gv2gml -o Makefile.gml Makefile.gv

# Render a GraphViz file to PNG
#%.png: %.gv
#	dot -Tpng -o $*.png $*.gv


# .PHONY: clean_figs

# clean_figs:
# 	rm ../figures/*.pdf ../figures/*.png
