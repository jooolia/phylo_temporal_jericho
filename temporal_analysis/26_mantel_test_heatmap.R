
library(reshape2)
library(ggplot2)

mantel_tests <- read.csv("../results/mantel_tests_otus_and_env.csv")[,1:5]

## new df_with_stars replaced
mantel_tests_sig <- mantel_tests
mantel_tests_sig$normalized_16s_OTU_sig <- ""
mantel_tests_sig$normalized_16s_OTU_sig[grepl("\\*",  mantel_tests$normalized_16s_OTUs)] <- 0.05
mantel_tests_sig$normalized_16s_OTU_sig[grepl("\\*\\*", mantel_tests$normalized_16s_OTUs)] <- 0.01
mantel_tests_sig$normalized_16s_OTU_sig[grepl("\\*\\*\\*", mantel_tests$normalized_16s_OTUs)] <- 0.001

mantel_tests_sig$normalized_18s_OTUs_sig <- ""
mantel_tests_sig$normalized_18s_OTUs_sig[grepl("\\*",  mantel_tests$normalized_18s_OTUs)] <- 0.05
mantel_tests_sig$normalized_18s_OTUs_sig[grepl("\\*\\*", mantel_tests$normalized_18s_OTUs)] <- 0.01
mantel_tests_sig$normalized_18s_OTUs_sig[grepl("\\*\\*\\*", mantel_tests$normalized_18s_OTUs)] <- 0.001

mantel_tests_sig$normalized_MPL_OTUs_sig <- ""
mantel_tests_sig$normalized_MPL_OTUs_sig[grepl("\\*",  mantel_tests$normalized_MPL_OTUs)] <- 0.05
mantel_tests_sig$normalized_MPL_OTUs_sig[grepl("\\*\\*", mantel_tests$normalized_MPL_OTUs)] <- 0.01
mantel_tests_sig$normalized_MPL_OTUs_sig[grepl("\\*\\*\\*", mantel_tests$normalized_MPL_OTUs)] <- 0.001

mantel_tests_sig$normalized_gp23_OTUs_sig <- ""
mantel_tests_sig$normalized_gp23_OTUs_sig[grepl("\\*",  mantel_tests$normalized_gp23_OTUs)] <- 0.05
mantel_tests_sig$normalized_gp23_OTUs_sig[grepl("\\*\\*", mantel_tests$normalized_gp23_OTUs)] <- 0.01
mantel_tests_sig$normalized_gp23_OTUs_sig[grepl("\\*\\*\\*", mantel_tests$normalized_gp23_OTUs)] <- 0.001

## get rid of the stars

mantel_tests$normalized_16s_OTUs <- as.numeric(gsub(" ?\\*?", "", mantel_tests$normalized_16s_OTUs))
mantel_tests$normalized_18s_OTUs <- as.numeric(gsub(" ?\\*?", "", mantel_tests$normalized_18s_OTUs))
mantel_tests$normalized_MPL_OTUs <- as.numeric(gsub(" ?\\*?", "", mantel_tests$normalized_MPL_OTUs))
mantel_tests$normalized_gp23_OTUs <- as.numeric(gsub(" ?\\*?", "", mantel_tests$normalized_gp23_OTUs))


melted_cormat <- melt(mantel_tests)


mantel_tests_sig <- mantel_tests_sig[,c(1,6:9)]
names(mantel_tests_sig) <- c("X",  "normalized_16s_OTUs","normalized_18s_OTUs","normalized_MPL_OTUs", "normalized_gp23_OTUs")
melted_cormat_sig <- melt(mantel_tests_sig, id.vars = "X")
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
sig_pal <- c("#000000","#999999")

pdf("../figures/mantel_tests_heatmap.pdf",
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
                       name="Mantel tests") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(variable,
                X,
                label = value),
            color = "black",
            size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )+
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
                      labels = c("< 0.01",
                                 "< 0.05"))
dev.off()

