## starting here https://github.com/idekerlab/cy-rest-R/blob/develop/r_markdown/basic1.Rmd

library(RJSONIO)
library(igraph)
library(httr)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

# Basic settings
## this was customized in my cytoscape
port.number = 8080
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

make_cytoscape_network <- function (base.url, graph1, network_name) {
 # Name of this new style
 style.name = "EssentialityAndDegree"
 
 # Delete the existing style for fresh start...
 style.url = paste(base.url, "styles", sep="/")
 style.delete.url = paste(style.url, style.name, sep="/")
 DELETE(url=style.delete.url)
 
 # Define default values
 def.node.color <- list(
  visualProperty = "NODE_FILL_COLOR",
  value = "#999999"
 )
 
 def.node.size <- list(
  visualProperty = "NODE_SIZE",
  value = 10
 )
 
 def.node_label_colour <- list(
  visualProperty = "NODE_LABEL_COLOR",
  value = "#000000"
 )
 
 def.node.border.width <- list(
  visualProperty = "NODE_BORDER_WIDTH",
  value = 0
 )
 
 def.edge.width <- list(
  visualProperty = "EDGE_WIDTH",
  value = 8
 )
 
 def.edge.color <- list(
  visualProperty = "EDGE_STROKE_UNSELECTED_PAINT",
  value = "#aaaaaa"
 )
 
 def.edge.transparency = list(
  visualProperty="EDGE_TRANSPARENCY",
  value = 50
 )
 
 def.node.transparency = list(
  visualProperty="NODE_TRANSPARENCY",
  value = 200
 )
 
 line.style <- list(
  visualProperty = "EDGE_LINE_TYPE",
  value = "EQUAL_DASH"
 )
 
 background_colour <- list(
  visualProperty="NETWORK_BACKGROUND_PAINT",
  value = "#FFFFFF"
 )
 
 node_label_size <- list(
  visualProperty="NODE_LABEL_FONT_SIZE",
  value = "20"
 )
 
 
 # works: "DOT", ZIGZAG "LONG_DASH" "EQUAL_DASH"
 # nope :"SINE", "WAVE, "DASHED""dash"
 
 defaults <- list(def.node.color, def.node.size, 
                  def.edge.color, def.node.border.width, 
                  def.edge.width, def.node.transparency, 
                  def.node_label_colour,
                  def.edge.transparency,
                  background_colour,
                  node_label_size 
                  #, line.style
 )
 
 mappings = list()
 
 # Node Size Mapping
 min.degree = min(V(graph1)$degree)
 max.degree = max(V(graph1)$degree)
 
 point1 = list(
  value=min.degree,
  lesser= "10.0",
  equal="10.0",
  greater="10.0"
 )
 
 point2 = list(
  value=max.degree,
  lesser= "100.0",
  equal="100.0",
  greater="100.0"
 )
 
 node.size.continuous.points = list(point1, point2)
 
 node.size = list(
  mappingType="continuous",
  mappingColumn="degree",
  mappingColumnType="Double",
  visualProperty="NODE_SIZE",
  points = node.size.continuous.points
 )
 
 node.label = list(
  mappingType="passthrough",
  mappingColumn="name",
  mappingColumnType="String",
  visualProperty="NODE_LABEL"
 )
 

 
 ## would be nice to use my generated palette for this. 
 phyla.mappings = list()
 phyla.mappings[[1]] <- list(key = "Fusobacteria", value = "#D6491F")
 phyla.mappings[[2]] <- list(key = "Chloroflexi", value = "#D34641")
 phyla.mappings[[3]] <- list(key = "Tenericutes", value = "#864B30")
 phyla.mappings[[4]] <- list(key = "Cryptophyceae", value = "ZIGZAG")
 phyla.mappings[[5]] <- list(key = "Archaeplastida", value = "#E99449")
 phyla.mappings[[6]] <- list(key = "Firmicutes", value = "#8B2761")
 phyla.mappings[[7]] <- list(key = "Deferribacteres", value = "#939533")
 phyla.mappings[[8]] <- list(key = "Proteobacteria", value = "#C59B33")
 phyla.mappings[[9]] <- list(key = "Planctomycetes", value = "#C57FB8")
 phyla.mappings[[10]] <- list(key = "Cyanobacteria", value = "#589857")
 phyla.mappings[[11]] <- list(key = "unclassified 18S", value = "#3B3B3B")
 phyla.mappings[[12]] <- list(key = "unclassified 16S", value = "#A8A8A8")
 phyla.mappings[[13]] <- list(key = "Opisthokonta", value = "#6CC5CE")
 phyla.mappings[[14]] <- list(key = "SAR", value = "#DE2BA2")
 phyla.mappings[[15]] <- list(key = "Actinobacteria", value = "#193D32")
 phyla.mappings[[16]] <- list(key = "Bacteroidetes", value = "#7488C0")
 phyla.mappings[[17]] <- list(key = "Amoebozoa", value = "#585671")
 phyla.mappings[[18]] <- list(key = "Haptophyta", value = "#0AC757") 
 phyla.mappings[[19]] <- list(key = "Excavata", value = "#D16477") 
 
 node_colour_style = list(
  mappingType="discrete",
  mappingColumn="Phylum",
  mappingColumnType="String",
  visualProperty="NODE_FILL_COLOR",
  map = phyla.mappings
 )
 
 mappings = list(node.size, 
                 node.label,
                 node_colour_style
 )
 
 style <- list(title=style.name, defaults = defaults, mappings = mappings)
 style.JSON <- toJSON(style)
 
 style.url = paste(base.url, "styles", sep="/")
 POST(url=style.url, body=style.JSON, encode = "json")
 
 #### utility functions #####
 #
 # Returns edge attributes for member edges.
 #
 toCytoscape <- function (igraphobj) {
  # Extract graph attributes
  graph_attr = graph.attributes(igraphobj)
  
  # Extract nodes
  node_count = length(V(igraphobj))
  if('name' %in% list.vertex.attributes(igraphobj)) {
   V(igraphobj)$id <- V(igraphobj)$name
  } else {
   V(igraphobj)$id <- as.character(c(1:node_count))
  }
  
  nodes <- V(igraphobj)
  v_attr = vertex.attributes(igraphobj)
  v_names = list.vertex.attributes(igraphobj)
  
  nds <- array(0, dim=c(node_count))
  for(i in 1:node_count) {
   if(i %% 1000 == 0) {
    print(i)
   }
   nds[[i]] = list(data = mapAttributes(v_names, v_attr, i))
  }
  
  edges <- get.edgelist(igraphobj)
  edge_count = ecount(igraphobj)
  e_attr <- edge.attributes(igraphobj)
  e_names = list.edge.attributes(igraphobj)
  
  attr_exists = FALSE
  e_names_len = 0
  if(identical(e_names, character(0)) == FALSE) {
   attr_exists = TRUE
   e_names_len = length(e_names)
  }
  e_names_len <- length(e_names)
  
  eds <- array(0, dim=c(edge_count))
  for(i in 1:edge_count) {
   st = list(source=toString(edges[i,1]), target=toString(edges[i,2]))
   
   # Extract attributes
   if(attr_exists) {
    eds[[i]] = list(data=c(st, mapAttributes(e_names, e_attr, i)))
   } else {
    eds[[i]] = list(data=st)
   }
   
   if(i %% 1000 == 0) {
    print(i)
   }
  }
  
  el = list(nodes=nds, edges=eds)
  
  x <- list(data = graph_attr, elements = el)
  print("Done.  To json Start...")
  return (toJSON(x))
 }
 
 mapAttributes <- function(attr.names, all.attr, i) {
  attr = list()
  cur.attr.names = attr.names
  attr.names.length = length(attr.names)
  
  for(j in 1:attr.names.length) {
   if(is.na(all.attr[[j]][i]) == FALSE) {
    #       attr[j] = all.attr[[j]][i]
    attr <- c(attr, all.attr[[j]][i])
   } else {
    cur.attr.names <- cur.attr.names[cur.attr.names != attr.names[j]]
   }
  }
  names(attr) = cur.attr.names
  return (attr)
 }
 
 ## create network from graph
 
 cygraph <- toCytoscape(graph1)
 network.url = paste(base.url, "networks", sep="/")
 res <- POST(url=network.url, body=cygraph, encode="json")
 
 # Extract SUID of the new network
 network.suid = unname(fromJSON(rawToChar(res$content)))
 network.suid
 
 
 # Apply style
 apply.style.url = paste(base.url, "apply/styles", style.name , toString(network.suid), sep="/")
 GET(apply.style.url)
 
 # Apply force-directed layout --need this first
 # Tweak Layout parameters
 layout.params = list(
  name="unweighted",
  value=TRUE
 )
 layout.params.url = paste(base.url, "apply/layouts/kamada-kawai/parameters", sep="/")
 PUT(layout.params.url, body=toJSON(list(layout.params)), encode = "json")
 
 apply.layout.url = paste(base.url, "apply/layouts/kamada-kawai", toString(network.suid), sep="/")
 
 GET(apply.layout.url)
 
 layout.params = list(
  name="allegro-spring-electric",
  gravity=1600,
  gravityTypeSelection="Circular"
  #value=TRUE
 )
 
 ## want to save as a file...
 layout.params.url = paste(base.url, "apply/layouts/allegro-spring-electric/parameters", sep="/")
 PUT(layout.params.url, body=toJSON(list(layout.params)), encode = "json")
 
 apply.layout.url = paste(base.url, "apply/layouts/allegro-spring-electric", toString(network.suid), sep="/")
 GET(apply.layout.url)
 
 network.image.url = paste(
  base.url,
  "networks",
  toString(network.suid),
  "views/first.png",
  sep="/"
 )
 print(network.image.url)
 
 # Toggle graphics details
 lod.url = paste(base.url, "ui/lod", sep="/")
 PUT(lod.url)
 
 network.image.url_pdf = paste(
  base.url,
  "networks",
  toString(network.suid),
  "views/first.pdf",
  sep="/"
 )
 print(network.image.url_pdf)
 
 ## from python notebook https://github.com/idekerlab/cy-rest-python/blob/develop/basic/CytoscapeREST_Basic2.ipynb
 ## ha! that works!
 download.file(network.image.url, paste0("../figures/", network_name, ".png"))
 download.file(network.image.url_pdf, paste0("../figures/", network_name, ".pdf"))
 
}

make_cytoscape_network_with_modules <- function (base.url, graph1, network_name) {
 # Name of this new style
 style.name = "EssentialityAndMods"
 
 # Delete the existing style for fresh start...
 style.url = paste(base.url, "styles", sep="/")
 style.delete.url = paste(style.url, style.name, sep="/")
 DELETE(url=style.delete.url)
 
 # Define default values
 def.node.color <- list(
  visualProperty = "NODE_FILL_COLOR",
  value = "#999999"
 )
 
 def.node.size <- list(
  visualProperty = "NODE_SIZE",
  value = 10
 )
 
 def.node_label_colour <- list(
  visualProperty = "NODE_LABEL_COLOR",
  value = "#000000"
 )
 
 def.node.border.width <- list(
  visualProperty = "NODE_BORDER_WIDTH",
  value = 0
 )
 
 def.edge.width <- list(
  visualProperty = "EDGE_WIDTH",
  value = 8
 )
 
 def.edge.color <- list(
  visualProperty = "EDGE_STROKE_UNSELECTED_PAINT",
  value = "#aaaaaa"
 )
 
 def.edge.transparency = list(
  visualProperty="EDGE_TRANSPARENCY",
  value = 50
 )
 
 def.node.transparency = list(
  visualProperty="NODE_TRANSPARENCY",
  value = 200
 )
 
 line.style <- list(
  visualProperty = "EDGE_LINE_TYPE",
  value = "EQUAL_DASH"
 )
 
 background_colour <- list(
  visualProperty="NETWORK_BACKGROUND_PAINT",
  value = "#FFFFFF"
 )
 
 node_label_size <- list(
  visualProperty="NODE_LABEL_FONT_SIZE",
  value = "20"
 )
 
 
 # works: "DOT", ZIGZAG "LONG_DASH" "EQUAL_DASH"
 # nope :"SINE", "WAVE, "DASHED""dash"
 
 defaults <- list(def.node.color, def.node.size, 
                  def.edge.color, def.node.border.width, 
                  def.edge.width, def.node.transparency, 
                  def.node_label_colour,
                  def.edge.transparency,
                  background_colour,
                  node_label_size 
                  #, line.style
 )
 
 mappings = list()
 
 # Node Size Mapping
 min.degree = min(V(graph1)$degree)
 max.degree = max(V(graph1)$degree)
 
 point1 = list(
  value=min.degree,
  lesser= "10.0",
  equal="10.0",
  greater="10.0"
 )
 
 point2 = list(
  value=max.degree,
  lesser= "100.0",
  equal="100.0",
  greater="100.0"
 )
 
 node.size.continuous.points = list(point1, point2)
 
 node.size = list(
  mappingType="continuous",
  mappingColumn="degree",
  mappingColumnType="Double",
  visualProperty="NODE_SIZE",
  points = node.size.continuous.points
 )
 
 node.label = list(
  mappingType="passthrough",
  mappingColumn="name",
  mappingColumnType="String",
  visualProperty="NODE_LABEL"
 )
 

 ### now want to also plot the modules as a separate file
 number_of_modules <- length(unique(V(graph1)$membership))
 module_colours <- colorRampPalette(brewer.pal(n = 12, name = "Set3"))(number_of_modules) 
 #module_colours <- colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(number_of_modules)
 #print(module_colours)
 
 str(graph1)
 print(min(V(graph1)$degree))
 print(V(graph1)$membership)
 ## generate a colour for each module
 module_mappings = list()
 for (i in 1:length(module_colours)){
  module_mappings[[i]] <- list(key = as.character(paste0(unique(V(graph1)$membership)[i], ".0")), value = module_colours[i])
 }
 print(module_mappings)
 
 node_colour_style = list(
  mappingType="discrete",
  mappingColumn="membership",
  mappingColumnType="Double",
  visualProperty="NODE_FILL_COLOR",
  map = module_mappings
 )
 
 mappings = list(node.size, 
                 node.label,
                 node_colour_style
 )
 
 style <- list(title=style.name, defaults = defaults, mappings = mappings)
 style.JSON <- toJSON(style)
 
 style.url = paste(base.url, "styles", sep="/")
 POST(url=style.url, body=style.JSON, encode = "json")
 
 #### utility functions #####
 #
 # Returns edge attributes for member edges.
 #
 toCytoscape <- function (igraphobj) {
  # Extract graph attributes
  graph_attr = graph.attributes(igraphobj)
  
  # Extract nodes
  node_count = length(V(igraphobj))
  if('name' %in% list.vertex.attributes(igraphobj)) {
   V(igraphobj)$id <- V(igraphobj)$name
  } else {
   V(igraphobj)$id <- as.character(c(1:node_count))
  }
  
  nodes <- V(igraphobj)
  v_attr = vertex.attributes(igraphobj)
  v_names = list.vertex.attributes(igraphobj)
  
  nds <- array(0, dim=c(node_count))
  for(i in 1:node_count) {
   if(i %% 1000 == 0) {
    print(i)
   }
   nds[[i]] = list(data = mapAttributes(v_names, v_attr, i))
  }
  
  edges <- get.edgelist(igraphobj)
  edge_count = ecount(igraphobj)
  e_attr <- edge.attributes(igraphobj)
  e_names = list.edge.attributes(igraphobj)
  
  attr_exists = FALSE
  e_names_len = 0
  if(identical(e_names, character(0)) == FALSE) {
   attr_exists = TRUE
   e_names_len = length(e_names)
  }
  e_names_len <- length(e_names)
  
  eds <- array(0, dim=c(edge_count))
  for(i in 1:edge_count) {
   st = list(source=toString(edges[i,1]), target=toString(edges[i,2]))
   
   # Extract attributes
   if(attr_exists) {
    eds[[i]] = list(data=c(st, mapAttributes(e_names, e_attr, i)))
   } else {
    eds[[i]] = list(data=st)
   }
   
   if(i %% 1000 == 0) {
    print(i)
   }
  }
  
  el = list(nodes=nds, edges=eds)
  
  x <- list(data = graph_attr, elements = el)
  print("Done.  To json Start...")
  return (toJSON(x))
 }
 
 mapAttributes <- function(attr.names, all.attr, i) {
  attr = list()
  cur.attr.names = attr.names
  attr.names.length = length(attr.names)
  
  for(j in 1:attr.names.length) {
   if(is.na(all.attr[[j]][i]) == FALSE) {
    #       attr[j] = all.attr[[j]][i]
    attr <- c(attr, all.attr[[j]][i])
   } else {
    cur.attr.names <- cur.attr.names[cur.attr.names != attr.names[j]]
   }
  }
  names(attr) = cur.attr.names
  return (attr)
 }
 
 ## create network from graph
 
 cygraph <- toCytoscape(graph1)
 network.url = paste(base.url, "networks", sep="/")
 res <- POST(url=network.url, body=cygraph, encode="json")
 
 # Extract SUID of the new network
 network.suid = unname(fromJSON(rawToChar(res$content)))
 network.suid
 
 
 # Apply style
 apply.style.url = paste(base.url, "apply/styles", style.name , toString(network.suid), sep="/")
 GET(apply.style.url)
 
 # Apply force-directed layout --need this first
 # Tweak Layout parameters
 layout.params = list(
  name="unweighted",
  value=TRUE
 )
 layout.params.url = paste(base.url, "apply/layouts/kamada-kawai/parameters", sep="/")
 PUT(layout.params.url, body=toJSON(list(layout.params)), encode = "json")
 
 apply.layout.url = paste(base.url, "apply/layouts/kamada-kawai", toString(network.suid), sep="/")
 
 GET(apply.layout.url)
 
 layout.params = list(
  name="allegro-spring-electric",
  gravity=1600,
  gravityTypeSelection="Circular"
  #value=TRUE
 )
 
 ## want to save as a file...
 layout.params.url = paste(base.url, "apply/layouts/allegro-spring-electric/parameters", sep="/")
 PUT(layout.params.url, body=toJSON(list(layout.params)), encode = "json")
 
 apply.layout.url = paste(base.url, "apply/layouts/allegro-spring-electric", toString(network.suid), sep="/")
 GET(apply.layout.url)
 
 network.image.url = paste(
  base.url,
  "networks",
  toString(network.suid),
  "views/first.png",
  sep="/"
 )
 print(network.image.url)
 
 # Toggle graphics details
 lod.url = paste(base.url, "ui/lod", sep="/")
 PUT(lod.url)
 
 network.image.url_pdf = paste(
  base.url,
  "networks",
  toString(network.suid),
  "views/first.pdf",
  sep="/"
 )
 print(network.image.url_pdf)
 
 ## from python notebook https://github.com/idekerlab/cy-rest-python/blob/develop/basic/CytoscapeREST_Basic2.ipynb
 ## ha! that works!
 download.file(network.image.url, paste0("../figures/", network_name, "_with_modules.png"))
 download.file(network.image.url_pdf, paste0("../figures/", network_name, "_with_modules.pdf"))
}



## make graphs for all of these

graphml_for_cytoscape <- list.files(path = "../figures/", pattern = "*.graphml$", all.files = FALSE,
                                    full.names = FALSE, recursive = FALSE,
                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#graphml_for_cytoscape <- "network_LSA_normalized_18sfor_LSA.txtfirst20.txt.lsa.graphml"


## cytoscape needs to be open
for (i in 1:length(graphml_for_cytoscape)){
 graph_file <- read.graph(paste0("../figures/",graphml_for_cytoscape[i]),format="graphml")
 make_cytoscape_network(base.url, graph_file, graphml_for_cytoscape[i])
 make_cytoscape_network_with_modules(base.url, graph_file, graphml_for_cytoscape[i])
}


## make graphs like fig 1 in tara oceans paper
#graph_file <- read.graph("../results/LSA_tables/network_LSA_normalized_all_for_LSA.txt.lsa.graphml",format="graphml")

## test with one that would have all the types...
graph_file <- read.graph("../figures/network_cooccurence_spearman_merged_all_all_types_proportion_above0_01.graphml", format="graphml")

## Make plot of Number of edges vs. the tupe of relationship
library(stringr)
only_different <- delete.edges(graph_file, E(graph_file) [str_split(ends(graph_file, 1), "_")[[1]][1]== str_split(ends(graph_file, 1), "_", n=2)[[2]][1]])
## so want to find edges -euk to euk S18 to S18, Bac to euk, euk ot vir, 

edges <- as.data.frame(ends(graph_file, E(graph_file)))

weight <- E(graph_file)$LS
#weight <-  as.numeric(E(graph_file)$weight)

edges <- cbind(edges, weight)

edges_between_like <- edges[str_split_fixed(edges[,1], "_", n=2)[,1]== str_split_fixed(edges[,2], "_", n=2)[,1],]

euk_to_euk <- edges_between_like [str_split_fixed(edges_between_like [,1], "_", n=2)[,1]=="S18",]
bac_to_bac <- edges_between_like [str_split_fixed(edges_between_like [,1], "_", n=2)[,1]=="S16",]

edges_between_unlike <- edges[str_split_fixed(edges[,1], "_", n=2)[,1]!= str_split_fixed(edges[,2], "_", n=2)[,1],]

euk_to_bac <- edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="S18" |str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="S16" ,]
euk_to_bac <- euk_to_bac[str_split_fixed(euk_to_bac [,2], "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_bac [,2], "_", n=2)[,1]=="S16" ,]

euk_to_AVS <- edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="S18" |str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="AVS" ,]
euk_to_AVS <- euk_to_AVS[str_split_fixed(euk_to_AVS[,2], "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_AVS[,2], "_", n=2)[,1]=="AVS" ,]

euk_to_MPL <- edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="S18" |str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="MPL" ,]
euk_to_MPL <- euk_to_MPL[str_split_fixed(euk_to_MPL[,2], "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_MPL[,2], "_", n=2)[,1]=="MPL" ,]

bac_to_gp23 <- edges_between_unlike[str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="S16" |str_split_fixed(edges_between_unlike [,1], "_", n=2)[,1]=="gp23" ,]
bac_to_gp23 <- bac_to_gp23[str_split_fixed(bac_to_gp23[,2], "_", n=2)[,1]=="S16" |str_split_fixed(bac_to_gp23[,2], "_", n=2)[,1]=="gp23" ,]

## count negative vs. positive and do a histogram
count_edges_between_like <- c(pos=dim(filter(edges_between_like, weight > 0))[1],
                              neg=dim(filter(edges_between_like, weight < 0))[1])
count_edges_between_unlike <- c(pos=dim(filter(edges_between_unlike, weight > 0))[1],
                                neg=dim(filter(edges_between_unlike, weight < 0))[1])

count_euk_to_euk <- c(pos=dim(filter(euk_to_euk, weight > 0))[1],
                      neg=dim(filter(euk_to_euk, weight < 0))[1])

count_bac_to_bac <- c(pos=dim(filter(bac_to_bac, weight > 0))[1],
                      neg=dim(filter(bac_to_bac, weight < 0))[1])

count_euk_to_bac<- c(pos=dim(filter(euk_to_bac, weight > 0))[1],
                     neg=dim(filter(euk_to_bac, weight < 0))[1])

count_euk_to_AVS <- c(pos=dim(filter(euk_to_AVS, weight > 0))[1],
                      neg=dim(filter(euk_to_AVS, weight < 0))[1])

count_euk_to_MPL <- c(pos=dim(filter(euk_to_MPL, weight > 0))[1],
                      neg=dim(filter(euk_to_MPL, weight < 0))[1])

count_bac_to_gp23 <- c(pos=dim(filter(bac_to_gp23, weight > 0))[1],
                       neg=dim(filter(bac_to_gp23, weight < 0))[1])

counts_for_barplot <- rbind(count_edges_between_like,
                            count_edges_between_unlike,
                            count_euk_to_euk,
                            count_bac_to_bac,
                            count_euk_to_bac,
                            count_euk_to_AVS,
                            count_euk_to_MPL,
                            count_bac_to_gp23
) 
## test some negs
#counts_for_barplot[,2] <- c(30,50,5,10,20,30,50,60)

melted_counts <- melt(counts_for_barplot, value.name="Number_of_Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"

pdf("../figures/overall_network_barplot_of_edges_spearman.pdf")
ggplot(melted_counts, aes(x=Between_which_communities, y=Number_of_Edges, fill=Sign_of_relationship))+
 geom_bar(stat="identity")+
 theme_classic()+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

