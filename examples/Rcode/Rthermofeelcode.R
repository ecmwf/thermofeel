#################################################################################
# (C) Copyright 1996- ECMWF.                                                    #
#                                                                               #
# This software is licensed under the terms of the Apache Licence Version 2.0   #
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.          #
# In applying this licence, ECMWF does not waive the privileges and immunities  #
# granted to it by virtue of its status as an intergovernmental organisation    #
# nor does it submit to any jurisdiction.                                       #
#                                                                               #
# Code By Chloe Brimicombe                                                      #      
# January 2022                                                                  # 
# Acknowledgments Claudia Di Napoli                                             #
#################################################################################

#loading the libraries required
require(ncdf4)
require(raster)
require(rworldmap)
require(RColorBrewer)

#Remove all temporary files that are more than 24 hours old:
removeTmpFiles(h=24)
#Remove all temporary files currently in existence:
removeTmpFiles(h=0)
rm(list=ls())


#Read_In Function for when you require hours from a previous month
#data the month
#datapre the previous month data
read_in1<- function(data,datapre){
  ndata <- nlayers(data)
  npredata<- nlayers(datapre)
  data<-data[[1:(ndata-7)]]
  datapre <- datapre[[(npredata-6):npredata]] 
  tdata <- stack(datapre,data)
  return(tdata)
}
August2000 <- read_in1("month directory","previous month directory")
#Read_In Function for when you have a list of files to read in
read_in <- function(file){
  require(ncdf4)
  files <- list.files(file,pattern='*.nc',full.names=TRUE, recursive = TRUE)
  files.r <- apply(as.data.frame(files), FUN=brick, MARGIN=1)
  data <- stack(files.r)
  return(data)
}


#Calculating the Mean/Min/Max
#Yes you could do it as an if statement in 1 function but, I didn't do that!
daily_max<- function(month){
  nmonth<- nlayers(month)
  monthmax <- month[[1:(nmonth/24)]]
  i=1
  y=1
  while(i<= nlayers(month)){
    monthmax[[y]] <- calc(month[[i:(i+23)]], function(x){max(x, na.rm=T)})
    i<- i+24
    y<- y+1
  }
  return(monthmax)
}

#to calculate the sum of all the values in a day
daily_sum<- function(month){
  nmonth<- nlayers(month)
  monthsum <- month[[1:(nmonth/24)]]
  i=1
  y=1
  while(i<= nlayers(month)){
    monthsum[[y]] <- calc(month[[i:(i+23)]], function(x){sum(x, na.rm=T)})
    i<- i+24
    y<- y+1
  }
  return(monthsum)
}
#daily mean for monthly data
daily_mean<- function(month){
  nmonth<- nlayers(month)
  monthmean <- month[[1:(nmonth/24)]]
  i=1
  y=1
  while(i<= nlayers(month)){
    monthmean[[y]] <- calc(month[[i:(i+23)]], function(x){mean(x, na.rm=T)})
    i<- i+24
    y<- y+1
  }
  return(monthmean)
}






#Climatology Function Huge only for monthly
monthly_climatology(){
  #List all the files (min/max/mean) you want to include in the climatology, 
  #you have to do this to avoid storing the variables as global <<-
  DailyMeanAugust2000 <- brick("filename")
  #etc..
  
  #all<- stack(calc(DailyMeanAugust2000,function(x){max(x,na.rm=T)}),
        #etc.)
  
        return(all)
}


#masking examples using rworldmap

#masks
world<- getMap()
nigeria<- world[which(world$ISO3=="NGA"),]
mali<- world[which(world$ISO3=="MLI"),]
kenya <-world[which(world$ISO3 == "KEN"),]
burkina <- world[which(world$NAME =="Burkina Faso"),]
#to mask
my_data_masked <- mask(data,nigeria)


#plotting example with country lines using many different libraries
library(maps)
library(maptools)
countries <- map("world", plot=FALSE) 
countries <- map2SpatialLines(countries, proj4string = CRS("+proj=longlat"))
library(grid)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)

#It's important that ggplot2 is not used with this plotting method
# It overwrites a key element from latticeExtra

YourData = #reference a file of data here
cutpts = # a list of cuts for a scale i.e. c(1,2,3,4,5)

labs <- graticule::graticule_labels(lons, lats,
                                    xline = lons[2],
                                    yline = lats[2])
labsLon <- labs[labs$islon,]
labsLat <- labs[!labs$islon,]

lons <- pretty(c(xmin(YourData), xmax(YourData)))
lats <- pretty(c(ymin(YourData), ymax(YourData)))
grat=graticule::graticule(lons,lats,proj= projection(YourData))


plt1 <- rasterVis::levelplot(YourData, 
                             main = "My Plot", at=cutpts ,margin=F)
plt1<- plt1 + latticeExtra::layer(sp.lines(countries, col="black", lwd=1))
plt1+    layer(sp.lines(grat)) +
  layer(sp.text(coordinates(labsLon),
                txt = parse(text = labsLon$lab),
                adj = c(1.1, -0.25),
                cex = 0.6)) +
  layer(sp.text(coordinates(labsLat),
                txt = parse(text = labsLat$lab),
                adj = c(-0.25, -0.25),
                cex = 0.6))


#If you want to make a gif of your netcdf layers
require(animation)

saveGIF({
  for(i in c(1:nlayers(YOURSTACKHERE))){
    l <- levelplot(YOURSTACKHERE[[i]],margin=FALSE)
    l<- l + latticeExtra::layer(sp.lines(countries, col="black", lwd=1))
    plot(l)
  }
}, interval=0.4, movie.name="YOURNAMEHERE.gif")



