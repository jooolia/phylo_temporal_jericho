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


## make graphs for all of these
graphml_for_cytoscape <- list.files(path = "../results/LSA_tables", pattern = "network_LSA.*.graphml", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
graphml_for_cytoscape <- graphml_for_cytoscape[!graphml_for_cytoscape=="network_LSA_overall.graphml"] 
graphml_for_cytoscape <- graphml_for_cytoscape[!graphml_for_cytoscape=="network_LSA_overall_between_unlike.graphml"] 
graphml_for_cytoscape <- graphml_for_cytoscape[!graphml_for_cytoscape=="network_LSA_overall_positive_edges.graphml"] 
## just want to test out small portion
 # graphml_for_cytoscape <-graphml_for_cytoscape[3]

## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../results/LSA_tables/",graphml_for_cytoscape[i]),format="graphml")
 print(graphml_for_cytoscape[i])
 make_cytoscape_network(base.url, graph_file, graphml_for_cytoscape[i])
 make_cytoscape_network_with_modules(base.url, graph_file, graphml_for_cytoscape[i])
}


saveCytoscapeSession(filepath="../figures/networks_overall_and_with_modules.cys")