
library(vegan)
library(ggplot2)
library(reshape2)

RdRp_PBD <- read.csv("../results/RdRp_only_miseq_phylogenetic_diversity_tablebeta_div.csv", row.names=1)
gp23_PBD <- read.csv("../results/gp23_only_miseq_phylogenetic_diversity_tablebeta_div.csv", row.names=1)
AVS_PBD <- read.csv("../results/AVS_only_miseq_phylogenetic_diversity_tablebeta_div.csv", row.names=1)
S18_PBD <- read.csv("../results/S18_phylogenetic_diversity_tablebeta_div.csv", row.names=1)
S16_PBD <- read.csv("../results/S16_phylogenetic_diversity_tablebeta_div.csv", row.names=1)

if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
 point_colour <- "black"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
  point_colour <- "white"
 }
}


## try with 18s and MPL
union_vcs_MPL_and_18s <- intersect(rownames(RdRp_PBD), rownames(S18_PBD))

Jericho_MPL <- droplevels(subset(RdRp_PBD, rownames(RdRp_PBD) %in% union_vcs_MPL_and_18s))
colnames(Jericho_MPL) <- gsub("X","", colnames(Jericho_MPL) )
Jericho_MPL_ready <- droplevels(subset(Jericho_MPL, select = colnames(Jericho_MPL) %in% union_vcs_MPL_and_18s))

Jericho_18s <- droplevels(subset(S18_PBD, rownames(S18_PBD) %in% union_vcs_MPL_and_18s))
colnames(Jericho_18s) <- gsub("X","", colnames(Jericho_18s) )
Jericho_18s_ready <- droplevels(subset(Jericho_18s, select = colnames(Jericho_18s) %in% union_vcs_MPL_and_18s))

mantel(Jericho_18s_ready, Jericho_MPL_ready)
mantel(Jericho_18s_ready, Jericho_MPL_ready,method="spear")


S18_diss <- melt(as.matrix(Jericho_18s_ready))

MPL_diss <- melt(as.matrix(Jericho_MPL_ready))

test_beta <- merge(S18_diss, MPL_diss, by=c("Var1", "Var2"))
#test_beta <- filter(test_beta, value.x >0 & value.y >0)

lm_eqn = function(m) {
 
 l <- list(a = format(coef(m)[1], digits = 2),
           b = format(abs(coef(m)[2]), digits = 2),
           r2 = format(summary(m)$r.squared, digits = 3));
 
 if (coef(m)[2] >= 0)  {
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
 } else {
  eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
 }
 
 as.character(as.expression(eq));                 
}

ggplot(test_beta, aes(x=value.x, y=value.y))+geom_point(colour=point_colour)+theme_JAG_presentation()
ggplot(test_beta, aes(x=value.x, y=value.y))+geom_point(colour=point_colour)+stat_smooth(method="lm", se=FALSE) +geom_text(aes(x = min(value.x), y = max(value.y), colour=point_colour, label = lm_eqn(lm(value.y ~ value.x, test_beta))), parse = TRUE)+theme_JAG_presentation()


## S16 and gp23

## 
union_vcs_gp23_and_16s <- intersect(rownames(RdRp_PBD), rownames(S18_PBD))

Jericho_gp23 <- droplevels(subset(RdRp_PBD, rownames(RdRp_PBD) %in% union_vcs_gp23_and_16s))
colnames(Jericho_gp23) <- gsub("X","", colnames(Jericho_gp23) )
Jericho_gp23_ready <- droplevels(subset(Jericho_gp23, select = colnames(Jericho_gp23) %in% union_vcs_gp23_and_16s))

Jericho_16s <- droplevels(subset(S18_PBD, rownames(S18_PBD) %in% union_vcs_gp23_and_16s))
colnames(Jericho_16s) <- gsub("X","", colnames(Jericho_16s) )
Jericho_16s_ready <- droplevels(subset(Jericho_16s, select = colnames(Jericho_16s) %in% union_vcs_gp23_and_16s))

mantel(Jericho_16s_ready, Jericho_gp23_ready)
mantel(Jericho_16s_ready, Jericho_gp23_ready,method="spear")


S18_diss <- melt(as.matrix(Jericho_16s_ready))

gp23_diss <- melt(as.matrix(Jericho_gp23_ready))

test_beta <- merge(S18_diss, gp23_diss, by=c("Var1", "Var2"))
test_beta <- filter(test_beta, value.x >0 & value.y >0)


ggplot(test_beta, aes(x=value.x, y=value.y))+geom_point(colour=point_colour)+theme_JAG_presentation()
ggplot(test_beta, aes(x=value.x, y=value.y))+geom_point(colour=point_colour)+stat_smooth(method="lm", se=FALSE) +geom_text(aes(x = min(value.x), y = max(value.y),colour=point_colour, label = lm_eqn(lm(value.y ~ value.x, test_beta))), parse = TRUE)+theme_JAG_presentation()
