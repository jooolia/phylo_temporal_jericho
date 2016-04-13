#Purpose: Create a panel chart from Jericho Pier data that I collected over 1 year.
#Author: Julia Gustavsen
#Date Created: 2 September 2013
#Date modified: 20 October 2014
#To-do: Also make graphs by season, and with the blitzes so that I can look for patterns. Fix pH data so that it can be added back. Fix chlorophyl so that the scale is better for most of them. 

### want to add in error bars at some point. 

library(ggplot2)
library(scales)
library(grid)
library(reshape2)
library(gridExtra)
library(dplyr)
library(cowplot)



args <- commandArgs(TRUE)
inputFile <- args[1]

## test to see if input file is given, so I can decide whether to use this argument or the orginal one
if (!file_test("-f", inputFile)) {
 print("input theme not defined, using orginal one for manuscript.")
 source("../../JAG_manuscript_figure.R")
 path_colour <- "black"
 line_colour <- "black"
} else {
 print("Cool you passed a nice theme file to this script")
 source(inputFile)
 if (inputFile == "../../JAG_black_presentation.R"){
  path_colour <- "white"
  line_colour <- "white"
 }
}

season_line <- geom_vline(xintercept = as.numeric(c(as.Date("2010-03-22"), 
                                                    as.Date("2010-06-22"), 
                                                    as.Date("2010-09-22"),
                                                    as.Date("2010-12-22"),
                                                    as.Date("2011-03-22"),
                                                    as.Date("2011-06-22"))),
                          colour="grey",
                          size=1.5)
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)
date_scaling <-   scale_x_date(limits = c(as.Date("2011-06-20"), as.Date("2011-07-06")))


Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis.csv", row.names=1)

Jericho_data$Date <- as.Date(Jericho_data$Date)

irradiation_data <- read.csv("../../JerichoTimeSeries/data/ibis_shortwave_irradiance_2010_2011_date_fixed.csv",row.names =
                              1)
irradiation_data$Date <- as.Date(irradiation_data$Date)

Jericho_data <- merge(Jericho_data, irradiation_data[, c("Mean","Date")], by = c("Date", "Date"))
colnames(Jericho_data)[dim(Jericho_data)[2]]<- "Mean_irradiation"

###Create panel chart of all the parameters #####

#Jericho_enviro_and_bio_data <- Jericho_data[,c("Date",
#                                                 "Temperature_YSI",
#                                                 "Salinity_ppt_YSI",
#                                                 "Dissolved_oxygen_percent",
#                                                 "pH",
#                                                "Average_viral_abundance",
#                                                "Standard_error_viral_abundance", 
#                                                "Average_bacterial_abundance",
#                                                "Standard_error_bacterial_abundance",
#                                                "Average_chl_a",
#                                                "Standard_error_chl_a",
#                                                "Average_PO4",
#                                                "Standard_error_PO4",
#                                                "Average_SiO2",
#                                                "Standard_error_SiO2",
#                                                "Average_NO3_NO2",
#                                                "Standard_error_NO3_NO2",
#                                                "Secchi_disk_disappears",
#                                                "Secchi_disk_reappears"
#                                                #,
#                            #"Mean_irradiation"
#                            )]


#To relabel the facet names I need to make a separate list that will be used by the function (Jericho_labeller) below to rename the facets. 
Jericho_names <- list( "Temperature_YSI" = "Temperature (C)",
                       "Salinity_ppt_YSI" = "Salinity (psu)",
                       "Dissolved_oxygen_percent" = "Dissolved O2 (%)",
                       "pH",
                       "Average_viral_abundance" = "Viral abundance (viruses/mL)",
                       "Average_bacterial_abundance" = "Bacterial abundance (cells/mL)",
                       "Average_chl_a"="Chlorophyll a (ug/L)",
                       "Average_PO4" = "Phosphate (um/L)", 
                       "Average_SiO2"= "Silicate (um/L)",
                       "Average_NO3_NO2"= "Nitrate+Nitrite (um/L)"
)

Jericho_labeller <- function(variable,value){
 return(Jericho_names[value])
}

#Reshape the data into a long format so that I can use the facet wrapping on it.
Jericho_enviro_and_bio_long <-melt(Jericho_data,id=("Date"), na.rm= TRUE)

## For presentations use like this:
# source("JAG_black_presentation.R")
# Panel_chart_defaults <- ggplot(Jericho_enviro_and_bio_long, 
#                                aes(Date,value)) +
#   geom_line(colour=line_colour)+
#   geom_point(colour="grey50")+
#   facet_grid(variable~., 
#              scales = "free_y", 
#              labeller=Jericho_labeller) +
#   theme_JAG_presentation()+
#   theme(strip.text.y = element_text(angle=0),
#         axis.line = element_line(line_colour))



## option so that viral and bacterial plots have better scales
fancy_scientific <- function(l) {
 # turn in to character string in scientific notation
 l <- format(l, scientific = TRUE)
 # quote the part before the exponent to keep all the digits
 l <- gsub("^(.*)e+", "'\\1'e", l)
 # turn the 'e+' into plotmath format
 l <- gsub("e", "%*%10^", l)
 ## plus needs to be escaped with \\
 l <- gsub("\\+", "", l)
 # return this as an expression
 parse(text=l)
} 




### Plots of all variables for full data set. ### 
environmental_plot <- function (long_data_set, parameter, ylimits, ylabel) {
 Temperature_plot <- ggplot(droplevels(subset(long_data_set,
                                              variable == parameter)), 
                            aes(Date,as.numeric(value)))+
  spring_bloom_line+
  season_line+
  geom_line(colour=line_colour, size=1.5)+
  geom_point(colour="grey50", size=3)+
  date_scaling +
  ylab(ylabel)+
  scale_y_continuous(limits=ylimits, expand=c(0,0))+
  ggtitle(NULL)+
  xlab("Date")+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        axis.title.y = element_text(face="bold", size=13, colour = 'black', angle = 90, vjust = 0.5))
}

par(mar=c(0,0,0,0))
par(omd=c(0,0,0,0))
temperature_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Temperature_YSI", c(12,20), ylabel =  expression(Temperature~(degree ~ C)))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

salinity_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Salinity_ppt_YSI", c(5,15), "Salinity \n(psu)\n")+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)




viral_abundance_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long,
                                                   "Average_viral_abundance",
                                                   c(1,60000000),
                                                   "Viral\n abundance \n(viruses/mL)")+
 
 scale_y_continuous(labels=fancy_scientific,
                    limits=c(1000000,80000000),
                    expand=c(0,0))+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

bacterial_abundance_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_bacterial_abundance", c(1,60000000), "Bacterial \n abundance\n (cells/mL)")+scale_y_continuous(labels=fancy_scientific,limits=c(1,7000000), expand=c(0,0))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

chlorophyll_a_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_chl_a", c(0,200), expression("Chlorophyll a (" * mu * "g/L)")) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

chlorophyll_a_full_dataset_under_10 <- environmental_plot(Jericho_enviro_and_bio_long, "Average_chl_a", c(0,10), expression("Chlorophyll a (" * mu * "g/L)")) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)




phosphate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_PO4", c(0,1.5), expression("Phosphate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

silicate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_SiO2", c(40,80), expression("Silicate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

nitrate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_NO3_NO2", c(0,80), expression("Nitrate+Nitrite (" * mu * "M)"))+scale_y_continuous(limits=c(0,6))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

do_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Dissolved_oxygen_percent", c(50,150), "Dissolved \noxygen (%)")+scale_y_continuous(limits=c(50,150))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

ph_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "pH", c(7,10), "pH")+scale_y_continuous(limits=c(7,10))

# irradiation_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Mean_irradiation", c(0,350), "Mean_irradiation units?")+scale_y_continuous(expand=c(0,0), limits=c(0,350))


## widths are not the same, how to change?
## amazing!!
# Get the widths
gA <- ggplot_gtable(ggplot_build(temperature_full_dataset))
gB <- ggplot_gtable(ggplot_build(salinity_full_dataset))
gC <- ggplot_gtable(ggplot_build(viral_abundance_full_dataset))
gD <- ggplot_gtable(ggplot_build(bacterial_abundance_full_dataset ))
gE <- ggplot_gtable(ggplot_build(chlorophyll_a_full_dataset))
gF <- ggplot_gtable(ggplot_build(silicate_full_dataset))
gG <- ggplot_gtable(ggplot_build(phosphate_full_dataset))
gH <- ggplot_gtable(ggplot_build(nitrate_full_dataset ))
gI <- ggplot_gtable(ggplot_build(do_full_dataset))
gJ <- ggplot_gtable(ggplot_build(ph_full_dataset))
#gI <- ggplot_gtable(ggplot_build(irradiation_full_dataset))
maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3], 
                     gC$widths[2:3], gD$widths[2:3],
                     gE$widths[2:3], gF$widths[2:3], 
                     gG$widths[2:3], gH$widths[2:3],
                     gI$widths[2:3],gJ$widths[2:3]
)

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth
gG$widths[2:3] <- maxWidth
gH$widths[2:3] <- maxWidth
gI$widths[2:3] <- maxWidth
gJ$widths[2:3] <- maxWidth
# Arrange the charts
grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
             gI,gJ, 
             ncol=1)

### Print graph to file ###

pdf("../figures/Jericho_full_panel_chart.pdf", width = 15, height = 25, onefile = FALSE )
# grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
#              gI,gJ,
#              ncol=1)
plot_grid(
 gA, gB, gC, gD,gE, gF, gG, gH,
          gI,gJ,
          ncol=1,
          labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
          align = "v")

dev.off()


### Plots With standard error of the mean ####

environmental_plot_with_se <- function(data_frame, data_frame_mean, data_frame_se, ylimits,ylabel) { 
 
 param_with_se <- data.frame(Sample_name=data_frame[,1], 
                             Means=as.numeric(data_frame_mean),
                             param_se=as.numeric(data_frame_se))
 
 param_plot_line <- ggplot(param_with_se, 
                           aes(Sample_name, Means))+
  spring_bloom_line+
  season_line+
  geom_line(colour=line_colour, size = 1.5)+
  geom_point(colour="grey50", size = 3) +
  date_scaling +
  ylab(ylabel)+
  scale_y_continuous(limits=ylimits, expand=c(0,0))+
  ggtitle(NULL)+
  xlab("Date")+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        axis.title.y = element_text(face="bold", size=13, colour = 'black'))
 
 param_plot_line
 param_for_graph_with_se <- aes(ymax = Means + param_se,
                                ymin= Means - param_se)
 param_plot_line_error <- param_plot_line +
  geom_errorbar(param_for_graph_with_se, 
                position="dodge", 
                width=0.25,
                colour=line_colour) 
 return(param_plot_line_error)
 
 
 
}

#missing dates in time series. 
missing_VA_date <- Jericho_data$Date[is.na(Jericho_data$Average_viral_abundance)]

viral_abundance_full_dataset_with_se <- environmental_plot_with_se(Jericho_data,
                                                                   Jericho_data[,"Average_viral_abundance"],
                                                                   Jericho_data[,"Standard_error_viral_abundance"], c(1,60000000), "Viral \n abundance \n(viruses/mL)")+
 scale_y_continuous(labels=fancy_scientific,limits=c(1000000,80000000), expand=c(0,0))+
 annotate("text", x = missing_VA_date , y = 1000000, label="*", size=6)+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)


#missing dates in time series. 
missing_BA_date <- Jericho_data$Date[is.na(Jericho_data$Average_bacterial_abundance)]

bacterial_abundance_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_bacterial_abundance"],Jericho_data[,"Standard_error_bacterial_abundance"], c(1,70000000), "Bacterial \n abundance\n (cells/mL)")+
 scale_y_continuous(labels=fancy_scientific,limits=c(1,7000000),
                    expand=c(0,0))+
 annotate("text", x = missing_BA_date , y = 1, label="*", size=6)+
 theme(axis.ticks.x = element_blank(),
       axis.text.x = element_blank(),
       axis.title.x=element_blank())+
 xlab(NULL)

chlorophyll_a_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_chl_a"],Jericho_data[,"Standard_error_chl_a"], c(0,200), expression("Chlorophyll a (" * mu * "g/L)")) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

chlorophyll_a_full_dataset_under_10_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_chl_a"],Jericho_data[,"Standard_error_chl_a"], c(0,200), expression("Chlorophyll a (" * mu * "g/L)")) +coord_cartesian(ylim=c(0, 10))


 #ylab(expression(mu * "g/L"))+
#  theme(panel.grid.major = element_blank(),
#        axis.text.x = element_text(size = 8, colour = "black"),
#        plot.background = element_rect(size=2,linetype="solid",color="grey50"),
#        axis.text.y = element_text(size = 8, colour = "black"))


vp <- grid::viewport(width = 0.5, 
                     height = 0.5,
                     x = 0.95,
                     y = 0.9, 
                     just = c("right",
                              "top"))
pushViewport(vp)

print(chlorophyll_a_full_dataset_with_se) 
print(chlorophyll_a_full_dataset_under_10_with_se, vp = vp)


# Create a transparent theme object
transparent_theme <- theme(
 axis.title.x = element_blank(),
 axis.title.y = element_text(size = 10, colour = "black"),
 #axis.text.x = element_blank(), 
 #axis.text.y = element_blank(),
 #axis.ticks = element_blank(),
 axis.text.x = element_text(size = 9, colour = "black"),
 axis.text.y = element_text(size = 9, colour = "black"),
 panel.grid = element_blank(),
 #axis.line = element_blank(),
 panel.background = element_rect(fill = "transparent",colour = NA),
 plot.background = element_rect(size=2,linetype="solid",color="grey50")
 #plot.background = element_rect(fill = "transparent",colour = NA)
 )

p1 <- chlorophyll_a_full_dataset_with_se# see previous sections for the scatterPlot

# Box plot of the x variable
p2 <- chlorophyll_a_full_dataset_under_10_with_se+
 transparent_theme+ ylab(expression("(" * mu * "g/L)"))



p2_grob = ggplotGrob(p2)

ymax1 <- 190
ymin1 <- 30
xmin1 <- min(Jericho_data$Date)+6
xmax1 <- max(Jericho_data$Date)+1.5

chlorophyll_with_small_one <- p1 + annotation_custom(grob = p2_grob,xmin = as.numeric(xmin1), xmax = as.numeric(xmax1),
                       ymin = as.numeric(ymin1), ymax = as.numeric(ymax1))




phosphate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_PO4"],Jericho_data[,"Standard_error_PO4"], c(0,1.5), expression("Phosphate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

silicate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_SiO2"],Jericho_data[,"Standard_error_SiO2"],  c(45,75), expression("Silicate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

nitrate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_NO3_NO2"],Jericho_data[,"Standard_error_NO3_NO2"],  c(0,20), expression("Nitrate+Nitrite (" * mu * "M)"))+scale_y_continuous(limits=c(0,6))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

# ggplotGrob(print(chlorophyll_a_full_dataset_under_10_with_se, vp = vp))

# pushViewport(viewport(layout = grid.layout(2, 2)))
# vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
# print(gA, vp = vplayout(1, 1:2))
# print(p2, vp = vplayout(2, 1))
# print(p3, vp = vplayout(2, 2))

gA <- ggplot_gtable(ggplot_build(temperature_full_dataset))
gB <- ggplot_gtable(ggplot_build(salinity_full_dataset))
gC <- ggplot_gtable(ggplot_build(viral_abundance_full_dataset_with_se))
gD <- ggplot_gtable(ggplot_build(bacterial_abundance_full_dataset_with_se ))
gE <- ggplot_gtable(ggplot_build(chlorophyll_a_full_dataset_with_se))
gE <- ggplot_gtable(ggplot_build(chlorophyll_with_small_one))
gF <- ggplot_gtable(ggplot_build(silicate_full_dataset_with_se))
gG <- ggplot_gtable(ggplot_build(phosphate_full_dataset_with_se))
gH <- ggplot_gtable(ggplot_build(nitrate_full_dataset_with_se ))
gI <- ggplot_gtable(ggplot_build(do_full_dataset))
gJ <- ggplot_gtable(ggplot_build(ph_full_dataset))
#gI <- ggplot_gtable(ggplot_build(irradiation_full_dataset))


maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3], 
                     gC$widths[2:3], gD$widths[2:3],
                     gE$widths[2:3], gF$widths[2:3], 
                     gG$widths[2:3], gH$widths[2:3],
                     gI$widths[2:3],gJ$widths[2:3])

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth
gG$widths[2:3] <- maxWidth
gH$widths[2:3] <- maxWidth
gI$widths[2:3] <- maxWidth
gJ$widths[2:3] <- maxWidth






# Arrange the charts
grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
             gI, gJ,
             ncol=1)

### Print graph to file ###

pdf("../figures/Jericho_full_panel_chart_with_se.pdf", width = 10, height = 20, onefile = FALSE )
# grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
#              gI,gJ,
#              ncol=1)

plot_grid(
 gA, gB, gC, gD,gE, gF, gG, gH,
 gI,gJ,
 ncol=1,
 labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
 align = "v")

dev.off()



pdf(
 "../figures/Jericho_env_with_se%03d.pdf", width = 20, height = 15, onefile = FALSE
)
temp_plot <- environmental_plot(Jericho_enviro_and_bio_long, "Temperature_YSI", c(12,20), ylabel =  expression(Temperature~(degree ~ C)))
temp_plot
sal_plot <- environmental_plot(Jericho_enviro_and_bio_long, "Salinity_ppt_YSI", c(5,15), "Salinity \n(psu)\n")
sal_plot
environmental_plot_with_se(Jericho_data,Jericho_data[,"Average_viral_abundance"], Jericho_data[,"Standard_error_viral_abundance"], c(1,60000000), "Viral \n abundance \n(viruses/mL)")+scale_y_continuous(labels=fancy_scientific,limits=c(1000000,80000000), expand=c(0,0))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_bacterial_abundance"],Jericho_data[,"Standard_error_bacterial_abundance"], c(1,70000000), "Bacterial \n abundance\n (cells/mL)")+scale_y_continuous(labels=fancy_scientific,limits=c(1,7000000), expand=c(0,0))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_chl_a"],Jericho_data[,"Standard_error_chl_a"], c(0,200), expression("Chlorophyll a (" * mu * "g/L)")) 

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_SiO2"],Jericho_data[,"Standard_error_SiO2"],  c(45,75), expression("Silicate (" * mu * "M)"))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_PO4"],Jericho_data[,"Standard_error_PO4"], c(0,1.5), expression("Phosphate (" * mu * "M)"))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_NO3_NO2"],Jericho_data[,"Standard_error_NO3_NO2"],  c(0,20), expression("Nitrate+Nitrite (" * mu * "M)"))+scale_y_continuous(expand=c(0,0), limits=c(0,6))

do_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Dissolved_oxygen_percent", c(50,150), "Dissolved \noxygen (%)")+scale_y_continuous(limits=c(50,150))
do_full_dataset

ph_full_dataset 

dev.off()

## Plot nutrients on one graph together ####

target <- c("Average_PO4",
            "Average_SiO2",
            "Average_NO3_NO2")

target_se <- c( "Standard_error_PO4",  "Standard_error_SiO2", "Standard_error_NO3_NO2")
nuts_data <- filter(Jericho_enviro_and_bio_long, variable %in% target)


nutrients_plot_line <- ggplot( nuts_data, 
                               aes(Date, as.numeric(value), group=variable, colour=variable))+ 
 spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(limits=c(0,80),expand=c(0,0))+ggtitle(NULL)+
 ylab("uM")+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())



pdf("../figures/Jericho_nutrients.pdf", width = 15, height = 11, onefile = FALSE )
nutrients_plot_line
dev.off()

## Plot chlorophyll as a small subsection ####



chla_data <- filter(Jericho_enviro_and_bio_long, variable == "Average_chl_a")


chla_plot_full_line <- ggplot( chla_data, 
                               aes(Date, as.numeric(value), group=variable, colour=variable))+ 
 spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(limits=c(0,200),expand=c(0,0))+ggtitle(NULL)+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())


chla_plot_mini_line <- ggplot( chla_data, 
                               aes(Date, as.numeric(value), group=variable, colour=variable))+ 
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(limits=c(0,10),expand=c(0,0))+ggtitle(NULL)+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())


pdf("../figures/Jericho_chl_a.pdf", width = 15, height = 11, onefile = FALSE )
grid.arrange(chla_plot_full_line,
             chla_plot_mini_line,
             ncol=1)
dev.off()


pdf("../figures/Jericho_pH.pdf", width = 15, height = 11, onefile = FALSE )
ph_plot_full_line <- ggplot( droplevels(subset(Jericho_enviro_and_bio_long, variable=="pH")), 
                             aes(Date, as.numeric(value), group=variable))+ 
 spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous()+ggtitle(NULL)+
 ylab("ph")+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())
ph_plot_full_line 
dev.off()

pdf("../figures/Jericho_DO.pdf", width = 15, height = 11, onefile = FALSE )
DO_plot_full_line <- ggplot( droplevels(subset(Jericho_enviro_and_bio_long, variable=="Dissolved_oxygen_percent")), 
                             aes(Date, as.numeric(value), group=variable))+ 
 spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(expand=c(0,0))+ggtitle(NULL)+
 ylab("Dissolved oxygen(%)")+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())
DO_plot_full_line
dev.off()


pdf("../figures/Jericho_secchi.pdf", width = 15, height = 11, onefile = FALSE )

sechi_plot_full_line <- ggplot( droplevels(subset(Jericho_enviro_and_bio_long,
                                                  (variable=="Secchi_disk_disappears"))),
                                aes(x=Date,
                                    y=as.numeric(value),
                                    group=variable))+
 #spring_bloom_line+
 #season_line+
 geom_jitter(position = position_jitter(width = .1),
             size=6,colour=line_colour) +
 scale_y_reverse(limits=c(5,0))+
 ggtitle(NULL)+
 theme_JAG_presentation()
sechi_plot_full_line

dev.off()


## Plot chlorophyll as a small subsection ####
chla_data <- filter(Jericho_enviro_and_bio_long,
                    variable == "Average_chl_a")

chla_plot_full_line <- ggplot(chla_data,
                              aes(Date, as.numeric(value),
                                  group = variable)) +
 season_line+
 spring_bloom_line+
 geom_line(size=1.5) +
 geom_point(colour = "grey50") +
 date_scaling+
 scale_y_continuous(limits = c(0,200),
                    expand = c(0,0))+ 
 ggtitle(NULL)+
 xlab(NULL)+
 ylab(expression("Chlorophyll a (" * mu * "g/L)"))+
 theme_JAG_presentation()+
 theme(panel.grid.major = element_blank()
       #axis.text.x = element_text(size = 8, colour = "black"),
       #plot.background = element_rect(size=2,linetype="solid",color="grey50"),
       #axis.text.y = element_text(size = 8, colour = "black")
 )


chla_plot_mini_line <- ggplot(chla_data,
                              aes(Date, as.numeric(value),
                                  group = variable)) +
 season_line+
 spring_bloom_line+
 geom_line(size=1.5) +
 geom_point(colour = "grey50") +
 date_scaling+
 #scale_y_continuous(limits = c(0,6),
 #                   expand = c(0,0))+
 coord_cartesian(ylim=c(0,10))+
 ggtitle(NULL)+
 ylab(expression("Chlorophyll a (" * mu * "g/L)"))+
 theme_JAG_presentation()+
 theme(panel.grid.major = element_blank(),
       plot.background = element_rect(size=2,linetype="solid",color="grey50"),
       axis.text.x = element_text(size = 8, colour = "black"),
       #plot.background = element_rect(size=2,linetype="solid",color="grey50"),
       axis.text.y = element_text(size = 8, colour = "black"))


vp <- grid::viewport(width = 0.5, 
                     height = 0.5,
                     x = 0.95,
                     y = 0.85, 
                     just = c("right",
                              "top"))

# print(chla_plot_mini_line) 
# print(chla_plot_full_line, vp = vp)

print(chla_plot_full_line) 
print(chla_plot_mini_line, vp = vp)


pdf(
 "../figures/Jericho_chl_a%03d.pdf", width = 15, height = 11, onefile = FALSE
)
#grid.arrange(chla_plot_full_line,
#             chla_plot_mini_line,
#             ncol = 1)
print(chla_plot_full_line) 
print(chla_plot_mini_line, vp = vp)
dev.off()


#############################################
### Do the whole time series ########
############################################

season_line <- geom_vline(xintercept = as.numeric(c(as.Date("2010-03-22"), 
                                                    as.Date("2010-06-22"), 
                                                    as.Date("2010-09-22"),
                                                    as.Date("2010-12-22"),
                                                    as.Date("2011-03-22"),
                                                    as.Date("2011-06-22"))),
                          colour="grey",
                          size=1.5)
spring_bloom_line <- geom_vline(xintercept = as.numeric(as.Date("2011-04-08")),
                                colour="green",
                                size=1)


# HAKA_bloom <- geom_rect(aes(xmin = as.Date("2011-06-20"),
#                             xmax = as.Date("2011-07-06"), 
#                             ymin = -Inf,
#                             ymax = Inf),
#                         alpha = 0.2,
#                         fill = "burlywood3",
#                         colour=NA) 

HAKA_bloom <- annotate("rect",xmin = as.Date("2011-06-20"),
                       xmax = as.Date("2011-07-06"), 
                       ymin = -Inf,
                       ymax = Inf,
                       alpha = 0.4,
                       fill = "burlywood3",
                       colour=NA) 



date_scaling <-  scale_x_date(breaks = date_breaks("month"), 
                              labels = date_format("%b'%y"),
                              limits = c(as.Date("2010-06-22"),
                                         as.Date("2011-07-19")))



Jericho_data <- read.csv("../results/Jericho_data_for_env_analysis_all_time_series.csv", row.names=1)

Jericho_data$Date <- as.Date(Jericho_data$Date)




irradiation_data <- read.csv("../../JerichoTimeSeries/data/ibis_shortwave_irradiance_2010_2011_date_fixed.csv",row.names =
                              1)
irradiation_data$Date <- as.Date(irradiation_data$Date)

Jericho_data <- merge(Jericho_data, irradiation_data[, c("Mean","Date")], by = c("Date", "Date"))
colnames(Jericho_data)[dim(Jericho_data)[2]]<- "Mean_irradiation"


#Reshape the data into a long format so that I can use the facet wrapping on it.
Jericho_enviro_and_bio_long <-melt(Jericho_data,id=("Date"), na.rm= TRUE)

## option so that viral and bacterial plots have better scales
fancy_scientific <- function(l) {
 # turn in to character string in scientific notation
 l <- format(l, scientific = TRUE)
 # quote the part before the exponent to keep all the digits
 l <- gsub("^(.*)e+", "'\\1'e", l)
 # turn the 'e+' into plotmath format
 l <- gsub("e", "%*%10^", l)
 ## plus needs to be escaped with \\
 l <- gsub("\\+", "", l)
 # return this as an expression
 parse(text=l)
} 

### Plots of all variables for full data set. ### 
environmental_plot <- function (long_data_set, parameter, ylimits, ylabel) {
 Temperature_plot <- ggplot(droplevels(subset(long_data_set,
                                              variable == parameter)), 
                            aes(Date,as.numeric(value)))+
  HAKA_bloom+
    spring_bloom_line+
  season_line+
  geom_line(colour=line_colour, size =1.5)+
  geom_point(colour="grey50", size = 3)+
  date_scaling +
  ylab(ylabel)+scale_y_continuous(limits=ylimits, expand=c(0,0))+ggtitle(NULL)+   xlab("Date")+
  theme_JAG_presentation()+
  theme(axis.title.y = element_text(face="bold", size=13, colour = 'black', angle = 90, vjust = 0.5))
}

par(mar=c(0,0,0,0))
par(omd=c(0,0,0,0))
temperature_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Temperature_YSI", c(0,25), ylabel =  expression(Temperature~(degree ~ C)))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

salinity_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Salinity_ppt_YSI", c(0,30), "Salinity \n(psu)\n")+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

viral_abundance_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_viral_abundance", c(1,60000000), "Viral \n abundance \n(viruses/mL)")+scale_y_continuous(labels=fancy_scientific,limits=c(1000000,80000000), expand=c(0,0))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

bacterial_abundance_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_bacterial_abundance", c(1,60000000), "Bacterial \n abundance\n (cells/mL)")+scale_y_continuous(labels=fancy_scientific,limits=c(1,7000000), expand=c(0,0))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

chlorophyll_a_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_chl_a", c(0,200), expression("Chlorophyll a (" * mu * "g/L)")) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

chlorophyll_a_full_dataset_under_10 <- environmental_plot(Jericho_enviro_and_bio_long, "Average_chl_a", c(0,10), expression("Chlorophyll a (" * mu * "g/L)")) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

phosphate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_PO4", c(0,3), expression("Phosphate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

silicate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_SiO2", c(0,80), expression("Silicate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

nitrate_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Average_NO3_NO2", c(0,80), expression("Nitrate+Nitrite (" * mu * "M)"))+scale_y_continuous(expand=c(0,0), limits=c(0,30))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

do_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Dissolved_oxygen_percent", c(50,150), "Dissolved \noxygen (%)")+scale_y_continuous(limits=c(50,150))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

missing_pH_date <- Jericho_data$Date[is.na(Jericho_data$pH)]

ph_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "pH", c(7,10), "pH")+scale_y_continuous(limits=c(7,10))+
 annotate("text", x = missing_pH_date , y = 7, label="*", size=6)


# irradiation_full_dataset <- environmental_plot(Jericho_enviro_and_bio_long, "Mean_irradiation", c(0,350), "Mean_irradiation units?")+scale_y_continuous(expand=c(0,0), limits=c(0,350))


## widths are not the same, how to change?
## amazing!!
# Get the widths
gA <- ggplot_gtable(ggplot_build(temperature_full_dataset))
gB <- ggplot_gtable(ggplot_build(salinity_full_dataset))
gC <- ggplot_gtable(ggplot_build(viral_abundance_full_dataset))
gD <- ggplot_gtable(ggplot_build(bacterial_abundance_full_dataset ))
gE <- ggplot_gtable(ggplot_build(chlorophyll_a_full_dataset))
gF <- ggplot_gtable(ggplot_build(silicate_full_dataset))
gG <- ggplot_gtable(ggplot_build(phosphate_full_dataset))
gH <- ggplot_gtable(ggplot_build(nitrate_full_dataset ))
gI <- ggplot_gtable(ggplot_build(do_full_dataset  ))
gJ <- ggplot_gtable(ggplot_build(ph_full_dataset ))
#gI <- ggplot_gtable(ggplot_build(irradiation_full_dataset))
maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3], 
                     gC$widths[2:3], gD$widths[2:3],
                     gE$widths[2:3], gF$widths[2:3], 
                     gG$widths[2:3], gH$widths[2:3],
                     gI$widths[2:3],gJ$widths[2:3]
)

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth
gG$widths[2:3] <- maxWidth
gH$widths[2:3] <- maxWidth
gI$widths[2:3] <- maxWidth
gJ$widths[2:3] <- maxWidth
# Arrange the charts
grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
             gI, gJ,
             ncol=1)

### Print graph to file ###

pdf("../figures/Jericho_all_time_series_full_panel_chart.pdf", width = 15, height = 20, onefile = FALSE )
# grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
#              gI,gJ,
#              ncol=1)
plot_grid(
 gA, gB, gC, gD,gE, gF, gG, gH,
 gI,gJ,
 ncol=1,
 labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
 align = "v")

dev.off()


### Plots With standard error of the mean ####

environmental_plot_with_se <- function(data_frame, data_frame_mean, data_frame_se, ylimits,ylabel) { 
 
 param_with_se <- data.frame(Sample_name=data_frame[,1], 
                             Means=as.numeric(data_frame_mean),
                             param_se=as.numeric(data_frame_se))
 
 param_plot_line <- ggplot(param_with_se, 
                           aes(Sample_name, Means))+ 
  HAKA_bloom+
    spring_bloom_line+
  season_line+
  geom_line(colour=line_colour, size = 1.5)+
  geom_point(colour="grey50", size = 3) +
  date_scaling +
  ylab(ylabel)+
  scale_y_continuous(limits=ylimits, expand=c(0,0))+
  ggtitle(NULL)+
  xlab("Date")+
  theme_JAG_presentation()+
  theme(panel.grid.major = element_blank(),
        axis.title.y = element_text(face="bold", size=13, colour = 'black', angle = 90, vjust = 0.5))
 
 param_plot_line
 param_for_graph_with_se <- aes(ymax = Means + param_se,
                                ymin= Means - param_se)
 param_plot_line_error <- param_plot_line +
  geom_errorbar(param_for_graph_with_se, 
                position="dodge", 
                width=0.25,
                colour=line_colour)
 return(param_plot_line_error)
}


#missing dates in time series. 
missing_VA_date <- Jericho_data$Date[is.na(Jericho_data$Average_viral_abundance)]


viral_abundance_full_dataset_with_se <- environmental_plot_with_se(Jericho_data,
                                                                   Jericho_data[,"Average_viral_abundance"],
                                                                   Jericho_data[,"Standard_error_viral_abundance"], c(1,60000000), "Viral \n abundance \n(viruses/mL)")+
 scale_y_continuous(labels=fancy_scientific,limits=c(1000000,80000000), expand=c(0,0))+
 annotate("text", x = missing_VA_date , y = 1000000, label="*", size=6)+
 theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

#missing dates in time series. 

missing_BA_date <- Jericho_data$Date[is.na(Jericho_data$Average_bacterial_abundance)]

bacterial_abundance_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_bacterial_abundance"],Jericho_data[,"Standard_error_bacterial_abundance"], c(1,70000000), "Bacterial \n abundance\n (cells/mL)")+ 
 annotate("text", x = missing_BA_date , y = 1, label="*", size=6)+
 scale_y_continuous(labels=fancy_scientific,limits=c(1,7000000), expand=c(0,0))+ 
 theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+
 xlab(NULL)

chlorophyll_a_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_chl_a"],Jericho_data[,"Standard_error_chl_a"], c(0,200), expression("Chlorophyll a (" * mu * "g/L)")) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

chlorophyll_a_full_dataset_under_10_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_chl_a"],Jericho_data[,"Standard_error_chl_a"], c(0,200), expression("Chlorophyll a (" * mu * "g/L)")) + coord_cartesian(ylim=c(0, 10))+ xlab(NULL)

# Create a transparent theme object
transparent_theme <- theme(
 axis.title.x = element_blank(),
 axis.title.y = element_text(size = 10, colour = "black"),
 #axis.text.x = element_blank(), 
 #axis.text.y = element_blank(),
 #axis.ticks = element_blank(),
 axis.text.x = element_text(size = 9, colour = "black"),
 axis.text.y = element_text(size = 9, colour = "black"),
 panel.grid = element_blank(),
 #axis.line = element_blank(),
 panel.background = element_rect(fill = "transparent",colour = NA),
 plot.background = element_rect(size=2,linetype="solid",color="grey50")
 #plot.background = element_rect(fill = "transparent",colour = NA)
)

p1 <- chlorophyll_a_full_dataset_with_se# see previous sections for the scatterPlot

# Box plot of the x variable
p2 <- chlorophyll_a_full_dataset_under_10_with_se+
 transparent_theme+ ylab("(ug/L)")

p2_grob = ggplotGrob(p2)

ymax1 <- 190
ymin1 <- 30
xmin1 <- min(Jericho_data$Date)+6
xmax1 <- max(Jericho_data$Date)-60

chlorophyll_with_small_one <- p1 + annotation_custom(grob = p2_grob,xmin = as.numeric(xmin1), xmax = as.numeric(xmax1),
                                                     ymin = as.numeric(ymin1), ymax = as.numeric(ymax1))
chlorophyll_with_small_one 


phosphate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_PO4"],Jericho_data[,"Standard_error_PO4"], c(0,4), expression("Phosphate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

silicate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_SiO2"],Jericho_data[,"Standard_error_SiO2"],  c(0,80), expression("Silicate (" * mu * "M)"))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)

nitrate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_NO3_NO2"],Jericho_data[,"Standard_error_NO3_NO2"],  c(0,30), expression("Nitrate+Nitrite (" * mu * "M)"))+scale_y_continuous(expand=c(0,0), limits=c(0,30))+ theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+xlab(NULL)


gA <- ggplot_gtable(ggplot_build(temperature_full_dataset))
gB <- ggplot_gtable(ggplot_build(salinity_full_dataset))
gC <- ggplot_gtable(ggplot_build(viral_abundance_full_dataset_with_se))
gD <- ggplot_gtable(ggplot_build(bacterial_abundance_full_dataset_with_se ))
gE <- ggplot_gtable(ggplot_build(chlorophyll_a_full_dataset_with_se))
gE <- ggplot_gtable(ggplot_build(chlorophyll_with_small_one))
gF <- ggplot_gtable(ggplot_build(silicate_full_dataset_with_se))
gG <- ggplot_gtable(ggplot_build(phosphate_full_dataset_with_se))
gH <- ggplot_gtable(ggplot_build(nitrate_full_dataset_with_se ))
gI <- ggplot_gtable(ggplot_build(do_full_dataset  ))
gJ <- ggplot_gtable(ggplot_build(ph_full_dataset ))
#gI <- ggplot_gtable(ggplot_build(irradiation_full_dataset))

maxWidth = unit.pmax(gA$widths[2:3], gB$widths[2:3], 
                     gC$widths[2:3], gD$widths[2:3],
                     gE$widths[2:3], gF$widths[2:3], 
                     gG$widths[2:3], gH$widths[2:3],
                     gI$widths[2:3],gJ$widths[2:3])

# Set the widths
gA$widths[2:3] <- maxWidth
gB$widths[2:3] <- maxWidth
gC$widths[2:3] <- maxWidth
gD$widths[2:3] <- maxWidth
gE$widths[2:3] <- maxWidth
gF$widths[2:3] <- maxWidth
gG$widths[2:3] <- maxWidth
gH$widths[2:3] <- maxWidth
gI$widths[2:3] <- maxWidth
gJ$widths[2:3] <- maxWidth
# Arrange the charts
grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
             gI,gJ, 
             ncol=1)

### Print graph to file ###

pdf("../figures/Jericho_all_time_series_full_panel_chart_with_se.pdf", width = 15, height = 20, onefile = FALSE )
# grid.arrange(gA, gB, gC, gD,gE, gF, gG, gH,
#              gI,gJ,
#              ncol=1)

plot_grid(
 gA, gB, gC, gD,gE, gF, gG, gH,
 gI,gJ,
 ncol=1,
 labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
 align = "v")


dev.off()



pdf(
 "../figures/Jericho_all_time_series_env_with_se%03d.pdf", width = 20, height = 15, onefile = FALSE
)

environmental_plot(Jericho_enviro_and_bio_long, "Temperature_YSI", c(0,25), ylabel =  expression(Temperature~(degree ~ C)))

environmental_plot(Jericho_enviro_and_bio_long, "Salinity_ppt_YSI", c(0,30), "Salinity \n(psu)\n")

environmental_plot_with_se(Jericho_data,
                           Jericho_data[,"Average_viral_abundance"],
                           Jericho_data[,"Standard_error_viral_abundance"], c(1,60000000), "Viral \n abundance \n(viruses/mL)")+scale_y_continuous(labels=fancy_scientific,limits=c(1000000,80000000), expand=c(0,0))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_bacterial_abundance"],Jericho_data[,"Standard_error_bacterial_abundance"], c(1,70000000), "Bacterial \n abundance\n (cells/mL)")+scale_y_continuous(labels=fancy_scientific,limits=c(1,7000000), expand=c(0,0))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_chl_a"],Jericho_data[,"Standard_error_chl_a"], c(0,200), expression("Chlorophyll a (" * mu * "g/L)"))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_SiO2"],Jericho_data[,"Standard_error_SiO2"],  c(0,80), expression("Silicate (" * mu * "M)"))

environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_PO4"],Jericho_data[,"Standard_error_PO4"], c(0,4), expression("Phosphate (" * mu * "M)"))

nitrate_full_dataset_with_se <- environmental_plot_with_se(Jericho_data, Jericho_data[,"Average_NO3_NO2"],Jericho_data[,"Standard_error_NO3_NO2"],  c(0,30), expression("Nitrate+Nitrite (" * mu * "M)"))+scale_y_continuous(expand=c(0,0), limits=c(0,30))

environmental_plot(Jericho_enviro_and_bio_long, "Dissolved_oxygen_percent", c(50,150), "Dissolved \noxygen (%)")+scale_y_continuous(limits=c(50,150))
ph_full_dataset

dev.off()

## Plot nutrients on one graph together ####

target <- c("Average_PO4",
            "Average_SiO2",
            "Average_NO3_NO2")

target_se <- c( "Standard_error_PO4",  "Standard_error_SiO2", "Standard_error_NO3_NO2")
nuts_data <- filter(Jericho_enviro_and_bio_long, variable %in% target)


nutrients_plot_line <- ggplot( nuts_data, 
                               aes(Date, as.numeric(value), group=variable, colour=variable))+ 
 HAKA_bloom+
  spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(limits=c(0,80),expand=c(0,0))+ggtitle(NULL)+
 ylab("uM")+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())



pdf("../figures/Jericho_all_time_series_nutrients.pdf", width = 15, height = 11, onefile = FALSE )
nutrients_plot_line
dev.off()

## Plot chlorophyll as a small subsection ####
chla_data <- filter(Jericho_enviro_and_bio_long, variable == "Average_chl_a")


chla_plot_full_line <- ggplot( chla_data, 
                               aes(Date, as.numeric(value), group=variable, colour=variable))+ 
 HAKA_bloom+
  spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(limits=c(0,200))+ggtitle(NULL)+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())


chla_plot_mini_line <- ggplot( chla_data, 
                               aes(Date, as.numeric(value), group=variable, colour=variable))+ 
 HAKA_bloom+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(limits=c(0,10))+ggtitle(NULL)+
 theme_JAG_presentation()+
 theme(panel.grid.major.x = element_line(colour="grey20"),
       panel.grid.major.y=element_blank())


pdf("../figures/Jericho_all_time_series_chl_a.pdf", width = 15, height = 11, onefile = FALSE )
grid.arrange(chla_plot_full_line,
             chla_plot_mini_line,
             ncol=1)
dev.off()


pdf("../figures/Jericho_all_time_series_pH.pdf", width = 15, height = 11, onefile = FALSE )
ph_plot_full_line <- ggplot( droplevels(subset(Jericho_enviro_and_bio_long, variable=="pH")), 
                             aes(Date, as.numeric(value), group=variable))+ 
 HAKA_bloom+
  spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous()+ggtitle(NULL)+
 ylab("pH")+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())
ph_plot_full_line 
dev.off()

pdf("../figures/Jericho_all_time_series_DO.pdf", width = 15, height = 11, onefile = FALSE )
DO_plot_full_line <- ggplot( droplevels(subset(Jericho_enviro_and_bio_long, variable=="Dissolved_oxygen_percent")), 
                             aes(Date, as.numeric(value), group=variable))+ 
 HAKA_bloom+
  spring_bloom_line+
 season_line+
 geom_line()+
 geom_point(colour="grey50") +
 date_scaling +
 scale_y_continuous(limits=c(60,140))+ggtitle(NULL)+
 ylab("Dissolved oxygen(%)")+
 theme_JAG_presentation()+
 theme(panel.grid.major=element_blank())
DO_plot_full_line
dev.off()


pdf("../figures/Jericho_all_time_series_secchi.pdf", width = 15, height = 11, onefile = FALSE )

sechi_plot_full_line <- ggplot( droplevels(subset(Jericho_enviro_and_bio_long, (variable=="Secchi_disk_disappears"))), 
                                aes(x=Date,y=as.numeric(value), group=variable))+ 
 spring_bloom_line+
 season_line+
 geom_jitter(position = position_jitter(width = .1),size=6,colour=line_colour) +
 scale_y_reverse(limits=c(5,0))+
 ggtitle(NULL)+
 date_scaling +
 ylab("Secchi depth(m)")+
 theme_JAG_presentation()
sechi_plot_full_line

dev.off()


## Plot chlorophyll as a small subsection ####
chla_data <- filter(Jericho_enviro_and_bio_long,
                    variable == "Average_chl_a")

chla_plot_full_line <- ggplot(chla_data,
                              aes(Date, as.numeric(value),
                                  group = variable)) +
 HAKA_bloom+
 season_line+
 spring_bloom_line+
 geom_line(size=1.5) +
 geom_point(colour = "grey50") +
 date_scaling+
 scale_y_continuous(limits = c(0,200),
                    expand = c(0,0))+ 
 ggtitle(NULL)+
 ylab(expression("Chlorophyll a (" * mu * "g/L)"))+
 theme_JAG_presentation()+
 theme(panel.grid.major = element_blank())


chla_plot_mini_line <- ggplot(chla_data,
                              aes(Date, as.numeric(value),
                                  group = variable)) +
 HAKA_bloom+
 season_line+
 spring_bloom_line+
 geom_line(size=1.5) +
 geom_point(colour = "grey50") +
 date_scaling+
 #scale_y_continuous(limits = c(0,6),
 #                   expand = c(0,0))+
 coord_cartesian(ylim=c(0,10))+
 ggtitle(NULL)+
 xlab(NULL)+
 ylab(expression(mu * "g/L"))+
 theme_JAG_presentation()+
 theme(panel.grid.major = element_blank(),
       axis.text.x = element_text(size = 8, colour = "black"),
       plot.background = element_rect(size=2,linetype="solid",color="grey50"),
       axis.text.y = element_text(size = 8, colour = "black"))


vp <- grid::viewport(width = 0.6, 
                     height = 0.6,
                     x = 0.70,
                     y = 0.3, 
                     just = c("right",
                              "bottom"))
# 
# print(chla_plot_mini_line) 
# print(chla_plot_full_line, vp = vp)

print(chla_plot_full_line) 
print(chla_plot_mini_line, vp = vp)


pdf(
 "../figures/Jericho_all_time_series_chl_a%03d.pdf", width = 15, height = 11, onefile = FALSE
)
#grid.arrange(chla_plot_full_line,
#             chla_plot_mini_line,
#             ncol = 1)
print(chla_plot_full_line) 
print(chla_plot_mini_line, vp = vp)


dev.off()

