## import the lsa and get it ready for cytoscape
library(igraph)
library(plyr)
library(dplyr)
library(stringr)
library(qvalue)



taxonomic_16_and_18s <- read.csv("../results/cleaned_up_taxonomies_with_renamed_OTUs.csv")

## 
normalized_otu_all <- read.delim("../results/LSA_tables/normalized_all_with_env_for_LSA.txt", row.names=1, na.strings = "NA", stringsAsFactors = FALSE)

lsa_file <- read.delim("../results/LSA_tables/normalized_all_with_env_for_LSA.txt.lsa")

## need to do Q locally because not running with the LSA

calculated_q_value <- qvalue(lsa_file$P)
lsa_file$Q <- calculated_q_value$qvalues

## rerun with Q!!

good_lsa <- lsa_file %>%
 filter(Q < 0.02 ) %>% 
 filter(P < 0.01)  %>% 
 filter(abs(LS) > 0.5)


good_lsa$X <- gsub(".size.*.", "",good_lsa$X)
good_lsa$Y <- gsub(".size.*.", "",good_lsa$Y)

good_lsa$X <- gsub("18SOTU", "S18",good_lsa$X)
good_lsa$X <- gsub("16SOTU", "S16",good_lsa$X)

good_lsa$Y <- gsub("18SOTU", "S18",good_lsa$Y)
good_lsa$Y <- gsub("16SOTU", "S16",good_lsa$Y)

## for now remove AVS
good_lsa <-  good_lsa %>% 
 filter(str_split_fixed(X, "_", n=2)[,1]!= "AVSOTU" ) %>%
 filter(str_split_fixed(Y, "_", n=2)[,1]!= "AVSOTU" )


rownames(normalized_otu_all) <- gsub(".size.*.", "",rownames(normalized_otu_all))
rownames(normalized_otu_all) <- gsub("18SOTU", "S18",rownames(normalized_otu_all))
rownames(normalized_otu_all) <- gsub("16SOTU", "S16",rownames(normalized_otu_all))
## So from this want to pull out the interactions that I'm interested in
## virus only
## virus spec

good_lsa_delay_0 <- filter(good_lsa, Delay == 0) 
good_lsa_delay_pos <- filter(good_lsa, Delay > 0)
good_lsa_delay_neg <- filter(good_lsa, Delay < 0)

good_lsa_pos_cor <- filter(good_lsa, LS > 0)
good_lsa_neg_cor <- filter(good_lsa, LS < 0)

### between like ### 
lsa_between_like <- good_lsa[str_split_fixed(good_lsa$X, "_", n=2)[,1]==str_split_fixed(good_lsa$Y, "_", n=2)[,1],]

lsa_euk_to_euk <- lsa_between_like[str_split_fixed(lsa_between_like$X, "_", n=2)[,1]=="S18",]
lsa_bac_to_bac <- lsa_between_like[str_split_fixed(lsa_between_like$X, "_", n=2)[,1]=="S16",]
lsa_MPL_to_MPL <- lsa_between_like[str_split_fixed(lsa_between_like$X, "_", n=2)[,1]=="MPLOTU",]
#lsa_AVS_to_AVS <- lsa_between_like[str_split_fixed(lsa_between_like$X, "_", n=2)[,1]=="AVSOTU",]
lsa_gp23_to_gp23<- lsa_between_like[str_split_fixed(lsa_between_like$X, "_", n=2)[,1]=="gp23OTU",]


### between unlike ####
lsa_between_unlike <- good_lsa[str_split_fixed(good_lsa$X, "_", n=2)[,1]!= str_split_fixed(good_lsa$Y, "_", n=2)[,1],]

euk_to_bac <- lsa_between_unlike[str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="S16" ,]
euk_to_bac <- euk_to_bac[str_split_fixed(euk_to_bac$Y, "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_bac$Y, "_", n=2)[,1]=="S16" ,]

#  euk_to_AVS <- lsa_between_unlike[str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="AVSOTU" ,]
#  euk_to_AVS <- euk_to_AVS[str_split_fixed(euk_to_AVS$Y, "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_AVS$Y, "_", n=2)[,1]=="AVSOTU" ,]

euk_to_MPL <- lsa_between_unlike[str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="MPLOTU" ,]
euk_to_MPL <- euk_to_MPL[str_split_fixed(euk_to_MPL$Y, "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_MPL$Y, "_", n=2)[,1]=="MPLOTU" ,]

# euk_to_MPL_to_AVS <- lsa_between_unlike[str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="MPLOTU"|str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="AVSOTU" ,] 
# euk_to_MPL_to_AVS <- euk_to_MPL_to_AVS[str_split_fixed(euk_to_MPL_to_AVS$Y, "_", n=2)[,1]=="S18" |str_split_fixed(euk_to_MPL_to_AVS$Y, "_", n=2)[,1]=="MPLOTU"|str_split_fixed(euk_to_MPL_to_AVS$Y, "_", n=2)[,1]=="AVSOTU" ,] 
#   
bac_to_gp23 <- lsa_between_unlike[str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="S16" | str_split_fixed(lsa_between_unlike$X, "_", n=2)[,1]=="gp23OTU" ,]
bac_to_gp23 <- bac_to_gp23[str_split_fixed(bac_to_gp23$Y, "_", n=2)[,1]=="S16" | str_split_fixed(bac_to_gp23$Y, "_", n=2)[,1]=="gp23OTU" ,]


## for AVS try with Chlorophyta, Dinophyta, Haptophyta, Heterokonta as in Zhong and Jacquet 2014.
OTUs_18s_compare_AVS <- subset(taxonomic_16_and_18s, Phylum == "SAR" | Phylum == "Haptophyta" | Order == "Chlorophyta")
OTUs_18s_compare_AVS <- subset(OTUs_18s_compare_AVS, !(Order == "Dinoflagellata" | Order == "Ciliophora" | Family == "Diatomea"))


#  lsa_euks_subset_for_AVS <- subset(lsa_between_unlike, (X %in% OTUs_18s_compare_AVS$otu_number|Y %in% OTUs_18s_compare_AVS$otu_number))
#  lsa_euks_subset_for_AVS <- lsa_euks_subset_for_AVS[str_split_fixed(lsa_euks_subset_for_AVS$X, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_euks_subset_for_AVS$X, "_", n=2)[,1]=="AVSOTU" ,]
#  lsa_euks_subset_for_AVS <- lsa_euks_subset_for_AVS[str_split_fixed(lsa_euks_subset_for_AVS$Y, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_euks_subset_for_AVS$Y, "_", n=2)[,1]=="AVSOTU" ,]
#  
## for MPL need diatoms, raphidophytes, Stramenopiles; Labyrinthulomycetes; Thraustochytriaceae

OTUs_18s_compare_MPL <- subset(taxonomic_16_and_18s, Phylum == "SAR" | Phylum == "Haptophyta" )
lsa_euks_subset_for_MPL <- subset(lsa_between_unlike, (X %in% OTUs_18s_compare_MPL$otu_number|Y %in% OTUs_18s_compare_MPL$otu_number))
lsa_euks_subset_for_MPL <- lsa_euks_subset_for_MPL[str_split_fixed(lsa_euks_subset_for_MPL$X, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_euks_subset_for_MPL$X, "_", n=2)[,1]=="MPLOTU" ,]
lsa_euks_subset_for_MPL <- lsa_euks_subset_for_MPL[str_split_fixed(lsa_euks_subset_for_MPL$Y, "_", n=2)[,1]=="S18" |str_split_fixed(lsa_euks_subset_for_MPL$Y, "_", n=2)[,1]=="MPLOTU" ,]


time_otu_table <- read.delim("../results/LSA_tables/normalized_all_with_env_for_LSA.txt", row.names=1)

times_in_graph_file <- colnames(time_otu_table)


make_igraph_for_cytoscape <- function (lsa_network, names_amplicons) {
 
 
  # lsa_network <-  head(good_lsa, n=200)
 
 ## need to add dates and stuff on to lsa netwrok
 ## just swithc the vlaues to 1 and then filter by the lsa columsn left...
 ### 4 aspect add in time
 ## want to add in the dates when the OTUs are found
 ## also want to add in dates when OTUs are found. 
 
 ## so for each date want to show if it is present or absent
 ## so could use an overall OTU table by date. Probably something I used to make the LSA in the first place
 normalized_otu_all_subset <- normalized_otu_all[(rownames(normalized_otu_all) %in% c(lsa_network$X,lsa_network$Y)),]
 # normalized_otu_all_subset$OTU_names <- rownames(normalized_otu_all_subset)
 
 #for each node you want to see if it is present at that date
 
 for (col in seq(ncol(normalized_otu_all_subset))){ 
  
  is_there <- c()
  merged_original_otus_1_col <- dplyr::select(normalized_otu_all_subset, col)
  merged_original_otus_1_col[,(names(normalized_otu_all_subset)[col])] <- as.numeric(merged_original_otus_1_col[,(names(normalized_otu_all_subset)[col])])
  merged_original_otus_1_col$OTU_names <- rownames(normalized_otu_all_subset) 
  col_name <- names(normalized_otu_all_subset)[col]
  ## because of non-standard evaluation...??
  call <- substitute(filter(merged_original_otus_1_col,
                            target > 1), list(target = as.name(col_name))) 
  filtered_original_otus_1_row <- eval(call) 
  
  
  ### want to add in this added up, so see how many times seen. :) Could do with matrix I guess. 
  
  for (row in seq(nrow(lsa_network))){
   row_of_interest <- lsa_network[row,]
   ## if these significantly correlated otus are found at one time point give add the date found to a data frame
   if ((row_of_interest$X  %in%  filtered_original_otus_1_row$OTU_names) & (row_of_interest$Y %in%  filtered_original_otus_1_row$OTU_names)){
    #print("yes")
    is_there <- c(is_there,1)
   }
   else{
    #print("no")
    is_there <- c(is_there,0)
   }
  }
  ## column bind to LSA network
  lsa_network <- cbind(lsa_network, is_there)
  names(lsa_network)[ncol(lsa_network)] <- names(filtered_original_otus_1_row)[1]
 }
 
 #str(lsa_network[, colnames(normalized_otu_all_subset)])
 lsa_network$edges_presence_width <- rowSums(lsa_network[, colnames(normalized_otu_all_subset)])
 ## add up 0s from the time col
 
 ### so for the edges for x date that are > 0 keep those rows that and then make network from that.
 
 for (date in seq_along(times_in_graph_file)){
  print(date)
  print(times_in_graph_file[date])
  
  ## need to use the standard evaluation of dplyr in this loop. 
  string_date <- paste(times_in_graph_file[date], ">",  "0")
  only_date_rows <- filter_(lsa_network, string_date)
if (nrow(only_date_rows) > 0){
 graph1 <- graph.data.frame(only_date_rows, directed=TRUE)
 
 V( graph1)$Domain <- as.character(taxonomic_16_and_18s$Domain[match(V( graph1)$name,
                                                                     taxonomic_16_and_18s$otu_number)])
 V( graph1)$Phylum <- as.character(taxonomic_16_and_18s$Phylum[match(V( graph1)$name,
                                                                     taxonomic_16_and_18s$otu_number)])
 V( graph1)$Class <- as.character(taxonomic_16_and_18s$Class[match(V( graph1)$name,
                                                                   taxonomic_16_and_18s$otu_number)])
 V( graph1)$Order <- as.character(taxonomic_16_and_18s$Order[match(V( graph1)$name,
                                                                   taxonomic_16_and_18s$otu_number)])
 V( graph1)$Family <- as.character(taxonomic_16_and_18s$Family[match(V( graph1)$name,
                                                                     taxonomic_16_and_18s$otu_number)])
 V( graph1)$Genus <- as.character(taxonomic_16_and_18s$Genus[match(V( graph1)$name,
                                                                   taxonomic_16_and_18s$otu_number)])
 
 
 V(graph1)$type <- ifelse(str_split_fixed(V(graph1)$name, "_", n=2)[,1]=="S18","S18",
                          ifelse(str_split_fixed(V(graph1)$name, "_", n=2)[,1]=="S16", 
                                 "S16",
                                 ifelse(str_split_fixed(V(graph1)$name, "_", n=2)[,1]=="MPLOTU", "MPL",
                                        ifelse(str_split_fixed(V(graph1)$name, "_", n=2)[,1]=="gp23OTU", "gp23", "other"))))
 
 
 E(graph1)$weight_LS <- only_date_rows$LS
 E(graph1)$delay <- only_date_rows$Delay
 
 # 2. Calculate some statistics and assign then to the graph
 
 graph1$name = paste0("Network (BA Model) LSA of ",names_amplicons)
 graph1$density = graph.density(graph1)
 V(graph1)$degree <- degree(graph1)
 V(graph1)$closeness <- closeness(graph1)
 V(graph1)$betweenness <- betweenness(graph1)
 V(graph1)$page_rank <- page.rank(graph1)$vector
 #  V(graph1)$community <- label.propagation.community(graph1)$membership
 E(graph1)$betweenness <- edge.betweenness(graph1)
 
 # 3 calculate communiites
 
 graph1$diameter <- diameter(graph1)
 graph1$farthest_nodes <- farthest.nodes(graph1)
 graph1$largest_clique <- largest.cliques(graph1)
 
 ## Many networks consist of modules which are densely connected themselves 
 ## but sparsely connected to other modules.
 wc <- cluster_edge_betweenness(graph1)
 wc <- cluster_walktrap(graph1)
 V(graph1)$membership <- wc$membership
 graph1$modularity <- modularity(wc)
 
 write.graph(graph1,
             file=paste("../results/LSA_tables/time_LSA_",
                        names_amplicons,
                        times_in_graph_file[date],
                        ".graphml",
                        sep=""),
             format="graphml")
}
 }
 
}

make_igraph_for_cytoscape(good_lsa_pos_cor, "overall_positive_edges")
make_igraph_for_cytoscape(good_lsa_neg_cor, "overall_negative_edges") 
make_igraph_for_cytoscape(lsa_between_like, "overall_between_like")
make_igraph_for_cytoscape(lsa_euk_to_euk, "euk_to_euk")
make_igraph_for_cytoscape(lsa_bac_to_bac, "bac_to_bac")
make_igraph_for_cytoscape(lsa_gp23_to_gp23, "gp23_to_gp23")
make_igraph_for_cytoscape(lsa_MPL_to_MPL, "MPL_to_MPL")
#make_igraph_for_cytoscape(lsa_AVS_to_AVS, "AVS_to_AVS")

make_igraph_for_cytoscape(lsa_between_unlike, "overall_between_unlike")
#make_igraph_for_cytoscape(euk_to_AVS, "euk_to_AVS")
make_igraph_for_cytoscape(euk_to_bac, "euk_to_bac")
make_igraph_for_cytoscape(euk_to_MPL, "euk_to_MPL")
make_igraph_for_cytoscape(bac_to_gp23, "bac_to_gp23")
#make_igraph_for_cytoscape(euk_to_MPL_to_AVS, "euk_to_MPL_to_AVS")

## subset taxonomically
#make_igraph_for_cytoscape(lsa_euks_subset_for_AVS, "euks_subset_for_AVS")
make_igraph_for_cytoscape(lsa_euks_subset_for_MPL, "euks_subset_for_MPL") #maybe too big?

make_igraph_for_cytoscape(good_lsa, "overall") 
make_igraph_for_cytoscape(good_lsa_delay_0, "delay_0")
make_igraph_for_cytoscape(good_lsa_delay_pos, "delay_pos")
make_igraph_for_cytoscape(good_lsa_delay_neg, "delay_neg")


