## make globe with water coloured in OCE. 

library(oce)

data(topoWorld)
plot(topoWorld)
plot(topoWorld, clatitude=30, clongitude=370, span=9000)

library(maps)
library(maptools)
data(wrld_simpl)
plot(wrld_simpl)
plot(wrld_simpl, col="white")

## colour from Jericho Pres pallette #50666F
plot(wrld_simpl, col="#f2f2f2",bg="#50666F", lwd=0.1)



png("../figures/world_map_for_black_background.png", width=1500, height=1500)
map("world",col = "#50666F", fill = TRUE,resolution = 0,lty = 0, bg="black")
#map("world",col = "white",add=TRUE,lty=1,lwd=0.01, bg="black")
dev.off()