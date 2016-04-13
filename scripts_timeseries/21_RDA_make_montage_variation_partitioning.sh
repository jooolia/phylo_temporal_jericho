## montage of rarefaction curves



  convert ../figures/Variation_partitioning_18s_ord_sum_chem_bio_time_richness_forward_selection.png -fill black \
-font Helvetica-bold -pointsize 14 \
label:' A) Eukaryotes' \
-gravity northwest -geometry +4+4 \
-composite ../figures/Variation_partitioning_18s_ord_sum_chem_bio_time_richness_forward_selection_annotated.png

  convert ../figures/Variation_partitioning_16s_ord_sum_chem_bio_time_richness_forward_selection.png -fill black \
-font Helvetica-bold -pointsize 14 \
label:' B) Bacteria ' \
-gravity northwest -geometry +4+4 \
-composite ../figures/Variation_partitioning_16s_ord_sum_chem_bio_time_richness_forward_selection_annotated.png

  convert ../figures/Variation_partitioning_MPL_group_sum_chem_bio_time_forward_selection.png -fill black \
-font Helvetica-bold -pointsize 14 \
label:' C) Marine picorna-like viruses ' \
-gravity northwest -geometry +4+4 \
-composite ../figures/Variation_partitioning_MPL_group_sum_chem_bio_time_forward_selection_annotated.png

  convert ../figures/Variation_partitioning_gp23_group_sum_chem_bio_time_forward_selection.png -fill black \
-font Helvetica-bold -pointsize 14 \
label:' D) T4-like myoviruses ' \
-gravity northwest -geometry +4+4 \
-composite ../figures/Variation_partitioning_gp23_group_sum_chem_bio_time_forward_selection_annotated.png



montage -mode concatenate -density 300 ../figures/Variation_partitioning_18s_ord_sum_chem_bio_time_richness_forward_selection_annotated.png ../figures/Variation_partitioning_16s_ord_sum_chem_bio_time_richness_forward_selection_annotated.png  ../figures/Variation_partitioning_MPL_group_sum_chem_bio_time_forward_selection_annotated.png ../figures/Variation_partitioning_gp23_group_sum_chem_bio_time_forward_selection_annotated.png ../figures/variation_partitioning_montage.png


