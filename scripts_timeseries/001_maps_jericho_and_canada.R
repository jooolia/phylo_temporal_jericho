
library(gridExtra)
library(reshape2)
library(plyr)

library(maps)
library(ggplot2)
library(mapdata)

args <- commandArgs(TRUE)
inputFile <- args[1]

inputFile <- "../../JAG_black_presentation.R"
## test to see if input file is given, so I can decide whether to use this argument or the orginal one. 
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

map.dat <- map_data("worldHires")
map.can <- droplevels(subset(map.dat, (region == "Canada" )) )
map.us <- droplevels(subset(map.dat, (region == "USA" )))
map.canus <- rbind(map.can, map.us)

y_coordi <- c(48.75, 51.5)
x_coordi <- c(-127.5,-122.5)

library(dplyr)

map.SOG <- map.canus %>%
 filter(., lat > y_coordi[1] & lat < y_coordi[2]) %>%
 filter(., long > x_coordi[1] & long < x_coordi[2])

## could add in better details here and then could use this...
station_map <- ggplot(aes(long,lat, group=group),data=map.SOG) + 
 #geom_polygon( fill="white")+
 geom_path(colour=path_colour)+
 theme_JAG_presentation()+
 theme(axis.text = element_blank(), axis.title=element_blank())
 station_map

Jericho_pier_lat <-  49.2768694
Jericho_pier_lon <- -123.2015

station_map+geom_text(aes(x=Jericho_pier_lon, y=Jericho_pier_lat),
          #data=vancouver,
          label="Jericho Pier", 
          colour="turquoise",
          size=9, vjust=-1)+ geom_point(aes(x=Jericho_pier_lon, y=Jericho_pier_lat), size=7, colour="turquoise")+
 theme(panel.grid.major =  element_blank())

 
 pdf(file="../figures/map_SOG.pdf")
 station_map
 dev.off()
 
 y_coordi <- c(25, 90)
x_coordi <- c(-155,-55)

map.Can_USA <- map.canus %>%
 filter(., lat > y_coordi[1] & lat < y_coordi[2]) %>%
 filter(., long > x_coordi[1] & long < x_coordi[2])

## could add in better details here and then could use this...
station_map <- ggplot(aes(long,lat, group=group), data=map.Can_USA)+
 #geom_polygon(fill="white")+
 geom_path(colour=path_colour, size=0.5)+
 theme_JAG_presentation() +
 theme(axis.text = element_blank(), axis.title=element_blank())
station_map

Van_lat <- 49.2827
Van_lon <- -123.1207

pdf(file="../figures/map_CAN_and_USA_with_Vancouver_labelled.pdf")
station_map + geom_text(aes(x=Van_lon, y=Van_lat),
                        #data=vancouver,
                        label="Vancouver", 
                        colour="turquoise",
                        size=9, vjust=-1)+ geom_point(aes(x=Van_lon, y=Van_lat), size=7, colour="turquoise")+
 theme(panel.grid.major =  element_blank())

dev.off()



### try out marmap

library(marmap)
library(ggplot2)

#  Fetch data on NOAA servers and write on disk
bat <- getNOAA.bathy(25, 90, -155,-55, res = 4, keep=TRUE)

# Create nice looking color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

# Plot
plot(bat, image = TRUE, land = TRUE, lwd = 0.1, bpal = list(c(0, max(bat), greys), c(min(bat), 0, blues)))
plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline

autoplot(bat, geom=c("raster", "contour"), colour=path_colour, size=0.1) + scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen")


## would also like to add vancouver here too. 

## SOG coordinates

SOG_y_coordi <- c(48.75, 51.5)
SOG_x_coordi <- c(-127.5,-122.5)

## Van coors

Van_lat <- 49.2827
Van_lon <- -123.1207

