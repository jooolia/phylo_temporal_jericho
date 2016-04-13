## histogram of degree distribution

library(igraph)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(cowplot)


graph_file <- read.graph("../results/LSA_tables/network_LSA_overall.graphml",format="graphml")

## degree distribution

degree_df <- as.data.frame(V(graph_file)$degree)
names(degree_df) <- "degree"
prop.table(table(degree_df))

## with degree under 5
sum(prop.table(table(degree_df))[1:5])


pdf("../figures/network_degree_histogram.pdf")
degree_overall <- ggplot(degree_df, aes(x=degree))+geom_histogram(binwidth=2)
degree_overall
dev.off()

## host -virus

graph_file <- read.graph("../results/LSA_tables/network_LSA_euk_to_MPL.graphml",format="graphml")

## degree distribution

degree_df <- as.data.frame(V(graph_file)$degree)
names(degree_df) <- "degree"
prop.table(table(degree_df))

## with degree under 5
sum(prop.table(table(degree_df))[1:5])

pdf("../figures/network_degree_histogram_euk_mpl.pdf")
degree_euk_mpl <- ggplot(degree_df, aes(x=degree))+geom_histogram(binwidth=2)
degree_euk_mpl
dev.off()

### gp23 - bac

graph_file <- read.graph("../results/LSA_tables/network_LSA_bac_to_gp23.graphml",format="graphml")

## degree distribution

degree_df <- as.data.frame(V(graph_file)$degree)
names(degree_df) <- "degree"
prop.table(table(degree_df))

## with degree under 5
sum(prop.table(table(degree_df))[1:5])

pdf("../figures/network_degree_histogram_bac_gp23.pdf")
degree_bac_gp23 <- ggplot(degree_df, aes(x=degree))+geom_histogram(binwidth=2)
degree_bac_gp23
dev.off()

### euk bac

graph_file <- read.graph("../results/LSA_tables/network_LSA_euk_to_bac.graphml",format="graphml")

## degree distribution

degree_df <- as.data.frame(V(graph_file)$degree)
names(degree_df) <- "degree"
prop.table(table(degree_df))

## with degree under 5
sum(prop.table(table(degree_df))[1:5])

pdf("../figures/network_degree_histogram_euk_bac.pdf")
degree_euk_bac <- ggplot(degree_df, aes(x=degree))+geom_histogram(binwidth=2)
degree_euk_bac
dev.off()


pdf("../figures/network_degree_histogram_all.pdf")
plot_grid(degree_overall,
          degree_euk_mpl,
          degree_bac_gp23,
          degree_euk_bac,
          ncol=2, 
          align="hv",
          labels = c("A) Overall                                  ",
                     "B) Eukaryotes to marine picorna-like viruses",
                     "C) Bacteria to T4-like myoviruses           ",
                     "D) Eukaryotes to bacteria                   "), hjust = -0.1, 
          label_size = 10
)
dev.off()


labels = c("A) Overall",
           "B) Eukaryotes to marine picorna-like viruses",
           "C) Bacteria to T4-like myoviruses",
           "D) Eukaryotes to bacteria")
