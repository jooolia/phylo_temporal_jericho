
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
library(phyloseq)
library(ggplot2)
library(plyr)
library(reshape2)

## read in normalized OTU table
## uses last row as VC
OTU_table_gp23 <- read.delim("../data/OTU_table_Jericho_time_series_normalized_gp23.tsv", row.names="VC_number")

## might have to ignore the VC_number last column?
test_otu <- otu_table(OTU_table_gp23, taxa_are_rows=FALSE)

### get relevant VC numbers
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)
row.names(Jericho_data) <- Jericho_data$VC_number

Library_metadata <- read.csv("../../JerichoAndSOGsequencing/Library_list_with_barcode_and_PCR_amplicons.csv", nrows=61)

env_data <- sample_data(subset(Jericho_data, select=-c(VC_number)))
# the NAs are therefore not part of the time series. Will deal with them later. 



sample_names(test_otu)
sample_names(env_data)
physeq <- phyloseq(test_otu, env_data)

#plot_bar(physeq)
plot_bar(physeq, x="season")
plot_bar(physeq, x="Date")
#plot_bar(physeq, facet_grid=~season)
plot_heatmap(physeq)
plot_heatmap(physeq, sample.order="Date", na.value="white")
#plot_richness(physeq)
plot_richness(physeq, x="season")
plot_richness(physeq, x="Date")
plot_richness(physeq, measures="shannon")
plot_net(physeq)

estimate_richness(physeq)
ig <- make_network(physeq, max.dist=0.9)
plot_network(ig, physeq, color="season",)

dist_methods <- unlist(distance("list"))
print(dist_methods)

# These require tree
dist_methods[(1:3)]
# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]
# This is the user-defined method:
dist_methods["designdist"]
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]


plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for (i in dist_methods) {
 # Calculate distance matrix
 iDist <- distance(physeq, method = i)
 # Calculate ordination
 iMDS <- ordinate(physeq, "MDS", distance = iDist)
 ## Make plot Don't carry over previous plot (if error, p will be blank)
 p <- NULL
 # Create plot, store as temp variable, p
 p <- plot_ordination(physeq, iMDS, color="season"
                   #   , color = "PO4", shape = "Chl_a"
                      )
 # Add title to each plot
 p <- p + ggtitle(paste("MDS using distance method ", i, sep = ""))
 # Save the graphic to file.
 plist[[i]] = p
}

## shade according to season...
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2 , color = season
                   #, shape = Enterotype
              )
           )
p = p + geom_point(size = 10, alpha = 0.5)
p = p + facet_wrap(~distance, scales = "free")
p = p + ggtitle("MDS on various distance metrics for MPL")
p

print(plist[["jsd"]])
print(plist[["jaccard"]])
print(plist[["bray"]])
print(plist[["gower"]])
print(plist[["chao"]])

chao_dist <- ordinate(physeq, method="PCoA", distance = "chao")

p <- plot_ordination(physeq,chao_dist, color="season"
                     #   , color = "PO4", shape = "Chl_a"
)+ geom_point(size = 10, alpha = 0.5)
# Add title to each plot
p <- p + ggtitle("MDS using distance method Chao with PCoA")
p


### try with 18s and the taxonomic stuff

normalized_18s_OTUs <- read.delim("../data/OTU_table_Jericho_time_series_normalized_18S.tsv", row.names="VC_number")

taxonomy_18s <- read.csv( "../results/cleaned_up_18s_taxonomy_Jericho.csv", row.names=1)

taxonomy_16s <- read.csv( "../results/cleaned_up_16s_taxonomy_Jericho.csv", row.names=1)

test_otu <- otu_table(normalized_18s_OTUs, taxa_are_rows=FALSE)

rownames(taxonomy_18s) <- taxonomy_18s$otu_number
taxonomy_18s <-taxonomy_18s[,2:7]
taxonomy_18s <- as.matrix(taxonomy_18s)
tax_mat <- tax_table(taxonomy_18s)

physeq <- phyloseq(test_otu, env_data, tax_mat)

plot_bar(physeq)
plot_bar(physeq, x="season")
plot_bar(physeq, x="Date")
plot_bar(physeq, fill="Genus")
plot_bar(physeq, fill="Domain")
plot_bar(physeq, fill="Phylum")
plot_bar(physeq, x="Date", fill="Phylum")+ geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(physeq,x="Date", facet_grid=~Phylum)
plot_heatmap(physeq)
plot_heatmap(physeq, sample.order="season", na.value="white")
plot_heatmap(physeq, taxa.order="Phylum", na.value="white", sample.order="season",low = "#66CCFF", high = "#000033")
plot_richness(physeq)
plot_richness(physeq, x="season")
plot_richness(physeq, x="Phylum")
plot_richness(physeq, x="Date")

estimate_richness(physeq)

plot_net(physeq)

