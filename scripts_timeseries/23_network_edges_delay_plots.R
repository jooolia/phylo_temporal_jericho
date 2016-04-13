## find number of edges between seasons
## find number of edges between sites X days apart


library(igraph)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)


graph_file <- read.graph("../results/LSA_tables/network_LSA_overall.graphml",format="graphml")


edges <- as.data.frame(ends(graph_file, E(graph_file)))

weight <- E(graph_file)$weight_LS

delay <- E(graph_file)$delay


edges <- cbind(edges, weight,delay)

pos_edges <- filter(edges, weight > 0)
neg_edges <- filter(edges, weight < 0)


## count negative vs. positive and do a histogram
count_pos_edges_lag <- c(pos=dim(filter(pos_edges, delay > 0))[1],
                         no_delay=dim(filter(pos_edges, delay == 0))[1],
                              neg=dim(filter(pos_edges, delay < 0))[1])
count_neg_edges_lag <- c(pos=dim(filter(neg_edges, delay > 0))[1],
                         no_delay=dim(filter(neg_edges, delay == 0))[1],
                         neg=dim(filter(neg_edges, delay < 0))[1])


counts_for_barplot <- rbind(count_pos_edges_lag,
                            count_neg_edges_lag
) 

melted_counts <- melt(counts_for_barplot, value.name="Number_of_Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"


## need to standardize this effort...by potential sample

pdf("../figures/overall_network_barplot_of_edges_by_lags.pdf", width = 15, height = 11)

ggplot(melted_counts, aes(x=Between_which_communities, y=Number_of_Edges, fill=Sign_of_relationship))+
 geom_bar(stat="identity")+
  ggtitle("Edges by  lags")+
 theme_classic()+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


## percent

percent_pos_edges_lag <- c(pos=(dim(filter(pos_edges, delay > 0))[1])/ecount(graph_file)*100,
                         no_delay=(dim(filter(pos_edges, delay == 0))[1])/ecount(graph_file)*100,
                         neg=(dim(filter(pos_edges, delay < 0))[1])/ecount(graph_file)*100)
percent_neg_edges_lag <- c(pos=(dim(filter(neg_edges, delay > 0))[1])/ecount(graph_file)*100,
                         no_delay=(dim(filter(neg_edges, delay == 0))[1])/ecount(graph_file)*100,
                         neg=(dim(filter(neg_edges, delay < 0))[1])/ecount(graph_file)*100)


counts_for_barplot <- rbind(percent_pos_edges_lag,
                            percent_neg_edges_lag
) 

melted_counts <- melt(counts_for_barplot, value.name="Percent_Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"


## need to standardize this effort...by potential sample

pdf("../figures/overall_network_barplot_of_edges_by_lags_percent.pdf", width = 15, height = 11)

ggplot(melted_counts, aes(x=Between_which_communities, y=Percent_Edges, fill=Sign_of_relationship))+
 geom_bar(stat="identity")+
 ggtitle("Edges by  lags overall")+
 theme_classic()+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


## percent by edge type (not overall number of edges)

percent_pos_edges_lag <- c(pos=(dim(filter(pos_edges, delay > 0))[1])/nrow(pos_edges)*100,
                           no_delay=(dim(filter(pos_edges, delay == 0))[1])/nrow(pos_edges)*100,
                           neg=(dim(filter(pos_edges, delay < 0))[1])/nrow(pos_edges)*100)
percent_neg_edges_lag <- c(pos=(dim(filter(neg_edges, delay > 0))[1])/nrow(neg_edges)*100,
                           no_delay=(dim(filter(neg_edges, delay == 0))[1])/nrow(neg_edges)*100,
                           neg=(dim(filter(neg_edges, delay < 0))[1])/nrow(neg_edges)*100)


counts_for_barplot <- rbind(percent_pos_edges_lag,
                            percent_neg_edges_lag
) 

melted_counts <- melt(counts_for_barplot, value.name="Percent_Edges")
names(melted_counts)[1] <- "Between_which_communities"
names(melted_counts)[2] <- "Sign_of_relationship"


## need to standardize this effort...by potential sample

pdf("../figures/overall_network_barplot_of_edges_by_lags_percent_each_type.pdf", width = 15, height = 11)

ggplot(melted_counts, aes(x=Between_which_communities, y=Percent_Edges, fill=Sign_of_relationship))+
 geom_bar(stat="identity")+
 ggtitle("Edges by  lags by type")+
 theme_classic()+
 theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
