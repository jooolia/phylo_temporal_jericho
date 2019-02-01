library(phyloseq)
library(ape)

normalized_MPL_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_MPL.tsv", row.names="VC_number")

OTU <- otu_table(normalized_MPL_OTUs , taxa_are_rows = FALSE)

tree <- read.tree("../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.tree")

#tree <- root(tree, "gi_29345244_Enterobacteria_phage_T4_", resolve.root=TRUE)

phy_test <- phyloseq(OTU, tree)

dist_unifrac <- UniFrac(phy_test, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
