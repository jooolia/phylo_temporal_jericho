library(igraph)

## want to make lists of which taxa the viral OTUs are alway seen with. 

graph_file <- read.graph("../results/LSA_tables/network_LSA_overall.graphml",format="graphml")

## find all the connections iwth mpl
MPL_subgraph <- induced.subgraph(graph=graph_file,
                                 vids=unlist(neighborhood(graph=graph_file,                                                                                       order=1,
                                                                           nodes= V(graph_file)$type == "MPL")
                                 )
)

## Ok then from the subgraph just determine all the 18s and their affiliations.
MPL_connected_18s_Phylum <- unique(V(MPL_subgraph)$Phylum[V(MPL_subgraph)$type == "S18"])
MPL_connected_18s_Class <- unique(V(MPL_subgraph)$Class[V(MPL_subgraph)$type == "S18"])
MPL_connected_18s_Family <- unique(V(MPL_subgraph)$Family[V(MPL_subgraph)$type == "S18"])
MPL_connected_18s_Genus <- unique(V(MPL_subgraph)$Genus[V(MPL_subgraph)$type == "S18"])

MPL_connected_16s_Phylum <- unique(V(MPL_subgraph)$Phylum[V(MPL_subgraph)$type == "S16"])
MPL_connected_16s_Class <- unique(V(MPL_subgraph)$Class[V(MPL_subgraph)$type == "S16"])
MPL_connected_16s_Family <- unique(V(MPL_subgraph)$Family[V(MPL_subgraph)$type == "S16"])
MPL_connected_16s_Genus <- unique(V(MPL_subgraph)$Genus[V(MPL_subgraph)$type == "S16"])



gp23_subgraph <- induced.subgraph(graph=graph_file,
                                 vids=unlist(neighborhood(graph=graph_file,                                                                                       order=1,
                                                          nodes= V(graph_file)$type == "gp23")
                                 )
)

## Ok then from the subgraph just determine all the 18s and their affiliations.
gp23_connected_18s_Phylum <- unique(V(gp23_subgraph)$Phylum[V(gp23_subgraph)$type == "S18"])
gp23_connected_18s_Class <- unique(V(gp23_subgraph)$Class[V(gp23_subgraph)$type == "S18"])
gp23_connected_18s_Family <- unique(V(gp23_subgraph)$Family[V(gp23_subgraph)$type == "S18"])
gp23_connected_18s_Genus <- unique(V(gp23_subgraph)$Genus[V(gp23_subgraph)$type == "S18"])

gp23_connected_16s_Phylum <- unique(V(gp23_subgraph)$Phylum[V(gp23_subgraph)$type == "S16"])
gp23_connected_16s_Class <- unique(V(gp23_subgraph)$Class[V(gp23_subgraph)$type == "S16"])
gp23_connected_16s_Family <- unique(V(gp23_subgraph)$Family[V(gp23_subgraph)$type == "S16"])
gp23_connected_16s_Genus <- unique(V(gp23_subgraph)$Genus[V(gp23_subgraph)$type == "S16"])


