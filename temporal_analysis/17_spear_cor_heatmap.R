
library(reshape2)
library(ggplot2)

spearman_correlations <- read.csv("../results/correlation_tests_otus_and_env.csv")[11:18,]
spearman_correlations <- spearman_correlations[!grepl("phyto", spearman_correlations$X),]
spearman_correlations <- spearman_correlations[!grepl("hetero", spearman_correlations$X),]
spearman_correlations <- spearman_correlations[,!grepl("phyto", colnames(spearman_correlations))]
spearman_correlations <- spearman_correlations[,!grepl("hetero", colnames(spearman_correlations))]
## had gotten rid of the diagonal, that is why it is 0 for the columns. 
spearman_correlations <- spearman_correlations[,!grepl("lagged", colnames(spearman_correlations))]


## need to copy 4 correlations  for the column MPL OTUs - the 16s and the gp23 comparisons

#spearman_correlations[1,"normalized_MPL_OTUs"] <- spearman_correlations[3,"normalized_16s_OTUs"]
#spearman_correlations[5,"normalized_MPL_OTUs"] <- spearman_correlations[3,"normalized_16s_OTUs_lagged"]

## new df_with_stars replaced
spearman_correlations_sig <- spearman_correlations
spearman_correlations_sig$Average_SiO2 <- ""
spearman_correlations_sig$Average_SiO2[grepl("\\*",  spearman_correlations$Average_SiO2)] <- 0.05
spearman_correlations_sig$Average_SiO2[grepl("\\*\\*", spearman_correlations$Average_SiO2)] <- 0.01
spearman_correlations_sig$Average_SiO2[grepl("\\*\\*\\*", spearman_correlations$Average_SiO2)] <- 0.001

spearman_correlations_sig$Average_chl_a <- ""
spearman_correlations_sig$Average_chl_a[grepl("\\*",  spearman_correlations$Average_chl_a)] <- 0.05
spearman_correlations_sig$Average_chl_a[grepl("\\*\\*", spearman_correlations$Average_chl_a)] <- 0.01
spearman_correlations_sig$Average_chl_a[grepl("\\*\\*\\*", spearman_correlations$Average_chl_a)] <- 0.001

spearman_correlations_sig$normalized_MPL_OTUs <- ""
spearman_correlations_sig$normalized_MPL_OTUs[grepl("\\*",  spearman_correlations$normalized_MPL_OTUs)] <- 0.05
spearman_correlations_sig$normalized_MPL_OTUs[grepl("\\*\\*", spearman_correlations$normalized_MPL_OTUs)] <- 0.01
spearman_correlations_sig$normalized_MPL_OTUs[grepl("\\*\\*\\*", spearman_correlations$normalized_MPL_OTUs)] <- 0.001

spearman_correlations_sig$normalized_gp23_OTUs <- ""
spearman_correlations_sig$normalized_gp23_OTUs[grepl("\\*",  spearman_correlations$normalized_gp23_OTUs)] <- 0.05
spearman_correlations_sig$normalized_gp23_OTUs[grepl("\\*\\*", spearman_correlations$normalized_gp23_OTUs)] <- 0.01
spearman_correlations_sig$normalized_gp23_OTUs[grepl("\\*\\*\\*", spearman_correlations$normalized_gp23_OTUs)] <- 0.001


spearman_correlations_sig <- subset(spearman_correlations_sig, select = c(X, Average_SiO2,Average_chl_a, normalized_MPL_OTUs,normalized_gp23_OTUs))


## get rid of the stars

spearman_correlations$normalized_16s_OTUs <- as.numeric(gsub(" ?\\*?", "", spearman_correlations$normalized_16s_OTUs))
spearman_correlations$normalized_18s_OTUs <- as.numeric(gsub(" ?\\*?", "", spearman_correlations$normalized_18s_OTUs))
spearman_correlations$normalized_MPL_OTUs <- as.numeric(gsub(" ?\\*?", "", spearman_correlations$normalized_MPL_OTUs))
spearman_correlations$normalized_gp23_OTUs <- as.numeric(gsub(" ?\\*?", "", spearman_correlations$normalized_gp23_OTUs))
spearman_correlations$Average_chl_a <- as.numeric(gsub(" ?\\*?", "", spearman_correlations$Average_chl_a))
spearman_correlations$Average_SiO2 <- as.numeric(gsub(" ?\\*?", "", spearman_correlations$Average_SiO2))

melted_cormat <- melt(spearman_correlations)
#melted_cormat$value[melted_cormat$value == 0] <- NA

#spearman_correlations_sig <- spearman_correlations_sig[,c(1,6:9)]
#names(spearman_correlations_sig) <- c("X",  "normalized_16s_OTUs","normalized_18s_OTUs","normalized_MPL_OTUs", "normalized_gp23_OTUs")
melted_cormat_sig <- melt(spearman_correlations_sig, id.vars = "X")
melted_cormat_sig$value <- as.numeric(melted_cormat_sig$value)


## make names more readable
melted_cormat$variable <- gsub("_", " ", melted_cormat$variable)
melted_cormat$X <- gsub("_", " ", melted_cormat$X)
melted_cormat_sig$variable <- gsub("_", " ", melted_cormat_sig$variable)
melted_cormat_sig$X <- gsub("_", " ", melted_cormat_sig$X)

## remove normalized but make sure to mention in figure legend
melted_cormat$variable <- gsub("normalized ", " ", melted_cormat$variable)
melted_cormat$X <- gsub("normalized ", " ", melted_cormat$X)
melted_cormat_sig$variable <- gsub("normalized ", " ", melted_cormat_sig$variable)
melted_cormat_sig$X <- gsub("normalized ", " ", melted_cormat_sig$X)


melted_cormat_sig <- melted_cormat_sig[!is.na(melted_cormat_sig$value),]
sig_pal <- c(
  #"#000000",
  "#999999"
)

pdf("../figures/spearman_correlations_heatmap.pdf",
    width=6)

ggplot(data = melted_cormat,
       aes(y = X,
           x = variable,
           fill = value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue",
                       high = "red",
                       mid = "white", 
                       midpoint = 0,
                       limit = c(-1,1),
                       space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal(base_size = 8)+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()+
  geom_text(aes(variable,
                X,
                label = value),
            color = "black",
            size = 2) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())+
  geom_tile(
    data = melted_cormat_sig,
    aes(y = X,
        x = variable,
        colour=as.factor(value)),
    width = 0.95,
    height = 0.95,
    alpha = 0,
    size = 0.5)+
  scale_colour_manual(values=sig_pal,
                      name = "P-value",
                      labels = c("< 0.05"))
dev.off()

