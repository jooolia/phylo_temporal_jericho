
## unzip file from nasa
## in folder home/h4tonccf

#bunzip2 A2014293.L3m_DAY_CHL_chlor_a_9km.bz2
#/home/julia/h4tonccf/h4tonccf_nc4 T2011174194000.L2_LAC_OC

## first ran >$ h4tonccf_nc4 INPUTFILE.hdf 
## this gave the .nc file

library(oce)
library(ncdf4)

args <- commandArgs(TRUE)
inputFile <- args[1]
plotFile <- args[2]

#bunzip2 ../data/M20111692011176.L3m_8D_CHL_chlor_a_4km.bz2 
#/home/julia/h4tonccf/h4tonccf_nc4 ../data/M20111692011176.L3m_8D_CHL_chlor_a_4km
#nc <- nc_open("../data/M20111692011176.nc")
## Finally works from here!! 
#http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MY1DMD_CHLORA&date=2011-06-18

inputFile <- "/home/julia/Downloads/A20151212015151.L3b_MO_CHL.nc"

nc <- nc_open("../data/M20111692011176.nc")

nc <- nc_open(inputFile) ## From June 18th

chl <- ncvar_get(nc, 'l3m_data')
chl[chl < 0] <- NA   # remove negative and missing values
lon <- ncvar_get(nc, "Number_of_Columns")
lat <- ncvar_get(nc, "Number_of_Lines")


data(coastlineWorld)
d <- read.oce("../data/ne_10m_admin_0_countries.shp")
## make a colormap that looks like the website
col <- colorRampPalette(c("purple", "#00007F", "blue", 
                          "#007FFF", "cyan","#7FFF7F", 
                          "yellow", "#FF7F00", "red", "#7F0000"))


pdf(plotFile, width = 15, height = 11, onefile = TRUE)
## in mg/m3, 
imagep(lon, lat, chl, col=col, zlim=c(0.01, 50), missingColor="white", zlab='Chl', axes=TRUE, xlim=c(-126,-122),ylim=c(48,50))
#polygon(coastlineWorld[['longitude']], coastlineWorld[['latitude']], col='grey')
#lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)
dev.off()

## works
# imagep(lon, lat, log10(chl), col=col, zlim=log10(c(0.01, 20)), missingColor=1, zlab='log10(Chl)', axes=TRUE, xlim=c(-126,-122),ylim=c(48,50))
# lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
# polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)



### Now do June 23rd ## should be day 174....
## download from nasa site
##bunzip2 ../data/A2011174.L3m_DAY_CHL_chlor_a_4km.bz2 
## /home/julia/h4tonccf/h4tonccf_nc4 ../data/A2011174.L3m_DAY_CHL_chlor_a_4km 

# nc <- nc_open("../data/A2011174.nc") ## From June 23rd

# chl <- ncvar_get(nc, 'l3m_data')
# chl[chl < 0] <- NA   # remove negative and missing values
# lon <- ncvar_get(nc, "Number_of_Columns")
# lat <- ncvar_get(nc, "Number_of_Lines")

# ## work has done above seeting this up
# ## but actually satellite did not pass by on that day...
# imagep(lon, lat, chl, col=col, zlim=c(0.01, 70), missingColor="white", zlab='Chl', axes=TRUE,xlim=c(-126,-122),ylim=c(48,50) )
# #lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
# polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)

# ### Try June 25th ### 
# ##bunzip2 ../data/A2011176.L3m_DAY_CHL_chlor_a_4km.bz2
# ## /home/julia/h4tonccf/h4tonccf_nc4 ../data/A2011176.L3m_DAY_CHL_chlor_a_4km
# nc <- nc_open("../data/A2011176.nc")
# chl <- ncvar_get(nc, 'l3m_data')
# chl[chl < 0] <- NA   # remove negative and missing values
# lon <- ncvar_get(nc, "Number_of_Columns")
# lat <- ncvar_get(nc, "Number_of_Lines")

# ## work has done above seeting this up
# ## but actually satellite did not pass by on that day...
# imagep(lon, lat, chl, col=col, zlim=c(0.01, 50), missingColor="white", zlab='Chl', axes=TRUE,xlim=c(-126,-122),ylim=c(48,50) )
# #lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
# polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)


# ### June 26th ### nothing!
# ## bunzip2 ../data/A2011177.L3m_DAY_CHL_chlor_a_4km.bz2
# ##/home/julia/h4tonccf/h4tonccf_nc4 ../data/A2011177.L3m_DAY_CHL_chlor_a_4km
# nc <- nc_open("../data/A2011177.nc")
# chl <- ncvar_get(nc, 'l3m_data')
# chl[chl < 0] <- NA   # remove negative and missing values
# lon <- ncvar_get(nc, "Number_of_Columns")
# lat <- ncvar_get(nc, "Number_of_Lines")

# ## work has done above seeting this up
# ## but actually satellite did not pass by on that day...
# imagep(lon, lat, chl, col=col, zlim=c(0.01, 50), missingColor="white", zlab='Chl', axes=TRUE,xlim=c(-126,-122),ylim=c(48,50) )
# #lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
# polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)

# ### June 27th #### nothing...

# ## bunzip2 ../data/A2011178.L3m_DAY_CHL_chlor_a_4km.bz2
# ##/home/julia/h4tonccf/h4tonccf_nc4 ../data/A2011178.L3m_DAY_CHL_chlor_a_4km
# nc <- nc_open("../data/A2011178.nc")
# chl <- ncvar_get(nc, 'l3m_data')
# chl[chl < 0] <- NA   # remove negative and missing values
# lon <- ncvar_get(nc, "Number_of_Columns")
# lat <- ncvar_get(nc, "Number_of_Lines")

# ## work has done above seeting this up
# ## but actually satellite did not pass by on that day...
# imagep(lon, lat, chl, col=col, zlim=c(0.01, 50), missingColor="white", zlab='Chl', axes=TRUE,xlim=c(-126,-122),ylim=c(48,50) )
# #lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
# polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)


# ## look at 8 Day average from June 18th-25th
# ## careful because gives wrong days!!!

# nc <- nc_open("../data/A20111692011176.nc")
# chl <- ncvar_get(nc, 'l3m_data')
# chl[chl < 0] <- NA   # remove negative and missing values
# lon <- ncvar_get(nc, "Number_of_Columns")
# lat <- ncvar_get(nc, "Number_of_Lines")

# ## work has done above seeting this up
# ## but actually satellite did not pass by on that day...
# imagep(lon, lat, chl, col=col, zlim=c(0.01, 50), missingColor="white", zlab='Chl', axes=TRUE,xlim=c(-126,-122),ylim=c(48,50) )
# #lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
# polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)


# ### Look at month of June


# nc <- nc_open("../data/A20111522011181.nc")
# chl <- ncvar_get(nc, 'l3m_data')
# chl[chl < 0] <- NA   # remove negative and missing values
# lon <- ncvar_get(nc, "Number_of_Columns")
# lat <- ncvar_get(nc, "Number_of_Lines")

# ## work has done above seeting this up
# ## but actually satellite did not pass by on that day...
# imagep(lon, lat, chl, col=col, zlim=c(0.01, 50), missingColor="white", zlab='Chl', axes=TRUE,xlim=c(-126,-122),ylim=c(48,50) )
# #lines(d[['longitude']], d[['latitude']], col='blue', lwd=1.4)
# polygon(d[['longitude']], d[['latitude']], col='black', lwd=1.4)

