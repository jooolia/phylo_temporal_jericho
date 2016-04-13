## starting here https://github.com/idekerlab/cy-rest-R/blob/develop/r_markdown/basic1.Rmd

library(RJSONIO)
library(igraph)
library(httr)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

load(file="../../phyla_palette_euks.txt")
load(file="../../class_palette_euks.txt")
load(file="../../order_palette_euks.txt")
load(file="../../family_palette_euks.txt")
load(file="../../phyla_palette_bac.txt")
load(file="../../class_palette_bac.txt")
load(file="../../order_palette_bac.txt")
load(file="../../family_palette_bac.txt")


time_otu_table <- read.delim("../results/LSA_tables/normalized_all_with_env_for_LSA.txt", row.names=1)

times_in_graph_file <- colnames(time_otu_table)

source("./cyrest_cytoscape_functions.R")
# Basic settings
## this was customized in my cytoscape
port.number = 1234
resetCytoscapeSession(port.number)
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")

print(base.url)

version.url = paste(base.url, "version", sep="/")
cytoscape.version = GET(version.url)
cy.version = fromJSON(rawToChar(cytoscape.version$content))
print(cy.version)

## look at default

default.style.url = paste(base.url, "styles/default", sep="/")
GET(url=default.style.url)
print(default.style.url)


top_10_MPL <- read.delim("../data/OTU_table_Jericho_time_series_MPL_normalized_top_10.tsv", row.names = 1)
top5_MPL_OTUs <- names(sort(colSums(top_10_MPL),decreasing = TRUE))[1:5]
top5_MPL_OTUs <- gsub("OTU", "MPLOTU", top5_MPL_OTUs)

top_10_gp23 <- read.delim("../data/OTU_table_Jericho_time_series_gp23_normalized_top_10.tsv", row.names = 1)
top5_gp23_OTUs <- names(sort(colSums(top_10_gp23),decreasing = TRUE))[1:5]
top5_gp23_OTUs <- gsub("OTU", "gp23OTU", top5_gp23_OTUs)


top_10_16s <- read.delim("../data/OTU_table_Jericho_time_series_16s_R1_normalized_top_10.tsv", row.names = 1)
top5_16s_OTUs <- names(sort(colSums(top_10_16s),decreasing = TRUE))[1:5]
top5_16s_OTUs <- gsub("OTU", "S16", top5_16s_OTUs )
top5_16s_OTUs <- gsub(".size.*.", "", top5_16s_OTUs)

top_10_18s <- read.delim("../data/OTU_table_Jericho_time_series_18s_normalized_top_10.tsv", row.names = 1)
top5_18s_OTUs <- names(sort(colSums(top_10_18s),decreasing = TRUE))[1:5]
top5_18s_OTUs <- gsub("OTU", "S18", top5_18s_OTUs )
top5_18s_OTUs <- gsub(".size.*.", "", top5_18s_OTUs)

top5_OTUs <- c(top5_MPL_OTUs,
               top5_gp23_OTUs,
               top5_16s_OTUs,
               top5_18s_OTUs)


## make graphs for all of these
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
# graphml_for_cytoscape <- graphml_for_cytoscape[!graphml_for_cytoscape=="network_LSA_overall.graphml"] 

## just want to test out small portion

## no need to keep rereading it
graph_file <- read.graph("../results/LSA_tables/network_LSA_overall_between_like.graphml",format="graphml")
  V(graph_file)$name
## cytoscape needs to be open
for (i in seq_along(top5_OTUs)){
 if (top5_OTUs[i] %in%   V(graph_file)$name){
 print(i)
 print(top5_OTUs[i])
 OTU_and_first_neighbours <- induced.subgraph(graph=graph_file,vids=unlist(neighborhood(graph=graph_file,order=1,nodes=top5_OTUs[i]))) ## node and first neighbours
 V(OTU_and_first_neighbours )$degree <- degree(OTU_and_first_neighbours)
#  test2 <- induced.subgraph(graph=graph_file,vids=unlist(neighborhood(graph=graph_file,order=2,nodes="S18_1"))) ## node and first and second neighbours
 # plot(OTU_and_first_neighbours ) 
# plot(test2)
 make_cytoscape_network(base.url, OTU_and_first_neighbours, paste0(top5_OTUs[i], "_between_like"))
 make_cytoscape_network_with_edges_sized_time(base.url, 
                                              OTU_and_first_neighbours,
                                              times_in_graph_file,
                                              paste0(top5_OTUs[i], "_between_like"))
 }
 }

graph_file <- read.graph("../results/LSA_tables/network_LSA_overall_between_unlike.graphml",format="graphml")
## cytoscape needs to be open
for (i in seq_along(top5_OTUs)){
 if (top5_OTUs[i] %in%  V(graph_file)$name){
 print(i)
 print(top5_OTUs[i])
 OTU_and_first_neighbours <- induced.subgraph(graph=graph_file,vids=unlist(neighborhood(graph=graph_file,order=1,nodes=top5_OTUs[i]))) ## node and first neighbours
 V(OTU_and_first_neighbours )$degree <- degree(OTU_and_first_neighbours)
 #  test2 <- induced.subgraph(graph=graph_file,vids=unlist(neighborhood(graph=graph_file,order=2,nodes="S18_1"))) ## node and first and second neighbours
 # plot(OTU_and_first_neighbours ) 
 # plot(test2)
 make_cytoscape_network(base.url, OTU_and_first_neighbours, paste0(top5_OTUs[i], "_between_unlike"))
 make_cytoscape_network_with_edges_sized_time(base.url, 
                                              OTU_and_first_neighbours,
                                              times_in_graph_file,
                                              paste0(top5_OTUs[i], "_between_unlike"))
 }
}

graph_file <- read.graph("../results/LSA_tables/network_LSA_overall.graphml",format="graphml")
## cytoscape needs to be open
for (i in seq_along(top5_OTUs)){
 if (top5_OTUs[i] %in%  V(graph_file)$name){
 print(i)
 print(top5_OTUs[i])
 OTU_and_first_neighbours <- induced.subgraph(graph=graph_file,vids=unlist(neighborhood(graph=graph_file,order=1,nodes=top5_OTUs[i]))) ## node and first neighbours
 V(OTU_and_first_neighbours )$degree <- degree(OTU_and_first_neighbours)
 #  test2 <- induced.subgraph(graph=graph_file,vids=unlist(neighborhood(graph=graph_file,order=2,nodes="S18_1"))) ## node and first and second neighbours
 # plot(OTU_and_first_neighbours ) 
 # plot(test2)
 make_cytoscape_network(base.url, OTU_and_first_neighbours, paste0(top5_OTUs[i], "_overall"))
 make_cytoscape_network_with_edges_sized_time(base.url, 
                                              OTU_and_first_neighbours,
                                              times_in_graph_file,
                                              paste0(top5_OTUs[i], "_overall"))
 }
}


saveCytoscapeSession(filepath="../figures/networks_top_5_otus.cys")
