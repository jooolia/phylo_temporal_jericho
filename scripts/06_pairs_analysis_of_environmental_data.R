library(ggplot2)
library(GGally)

## Read in data needed:
## calculated parameters:
Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)
Jericho_data$Date <- as.Date(Jericho_data$Date)
## amplicon richness calculated in earlier script.
amplicon_richness <- read.csv("../results/amplicon_richness_by_date.csv")
amplicon_richness$Date <- as.Date(amplicon_richness$Date)

## subset data that is interesting to me:
keep <- c("Date",
 "Average_viral_abundance",
          "Average_bacterial_abundance",
          "Average_chl_a",
          "Average_PO4",
          "Average_SiO2",
          "Average_NO3_NO2",
          "Temperature_YSI",
          "Salinity_ppt_YSI",
          "Dissolved_oxygen_percent",  
          "pH",
 "day_length")
## might be useful to do some more cleaning before
## keep only certain rows. 
params_Jericho <- droplevels(Jericho_data[,keep])


pairs(params_Jericho)
pairs(params_Jericho, log = "xy")

pdf("../figures/pairs_of_Jericho_parameter%03d.pdf", onefile = FALSE)
pairs(params_Jericho, main="Pairs of Jericho environmental data")
pairs(params_Jericho, log = "xy", main="Log-scaled Pairs of Jericho environmental data")
dev.off()

pdf("../figures/pairs_of_Jericho_parameter_with_regression_lines.pdf")
pairs( params_Jericho, panel=function(x,y){
 points(x,y)
 abline(lm(y~x), col='red')
}, 
main="Pairs of Jericho environmental data with linear regressions")

dev.off()

## Compare richness to the environmental parameters =================

## merge in to the environmental data excluding VC_number and season
environmental_and_amplicon_richness <- merge(params_Jericho, amplicon_richness, by="Date", all=TRUE)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
 usr <- par("usr"); on.exit(par(usr))
 par(usr = c(0, 1, 0, 1))
 r <- cor(x, y, use="complete.obs")
 txt <- format(c(r, 0.123456789), digits = digits)[1]
 txt <- paste0(prefix, txt)
 if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
 text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}
pairs(environmental_and_amplicon_richness, lower.panel=panel.smooth,upper.panel=panel.cor, main="Pairs of Jericho environmental data and Amplicon richness with correlations")

pdf("../figures/pairs_of_Jericho_parameter_and_amplicon_richness_with_regression_lines%03d.pdf", width=30, height=30, onefile = FALSE)
pairs( environmental_and_amplicon_richness, panel=function(x,y){
 points(x,y)
 abline(lm(y~x), col='red')
}, 
main="Pairs of Jericho environmental data and Amplicon richness with linear regressions")
pairs(environmental_and_amplicon_richness, lower.panel=panel.smooth,upper.panel=panel.cor, main="Pairs of Jericho environmental data and Amplicon richness with correlations")

dev.off()

## try deleting that high value of Chlorophyll
environmental_and_amplicon_richness$Average_chl_a[27] <- NA

pdf("../figures/pairs_of_Jericho_parameter_and_amplicon_richness_with_regression_lines_with_high_chlorophyll_deleted%03d.pdf", width=30, height=30, onefile = FALSE)
pairs( environmental_and_amplicon_richness, panel=function(x,y){
 points(x,y)
 abline(lm(y~x), col='red')
}, 
main="Pairs of Jericho environmental data and Amplicon richness with linear regressions")
pairs(environmental_and_amplicon_richness, lower.panel=panel.smooth,upper.panel=panel.cor, main="Pairs of Jericho environmental data and Amplicon richness with correlations")
dev.off()

pdf("../figures/pairs_of_Jericho_parameter_and_amplicon_richness_with_smoothing_with_high_chlorophyll_deleted.pdf", width=30, height=30)
pairs(environmental_and_amplicon_richness, panel=panel.smooth, main="Pairs of Jericho environmental data and Amplicon richness with smoothing")
dev.off()