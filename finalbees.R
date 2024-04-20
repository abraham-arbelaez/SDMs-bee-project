# Load required packages

library(pacman)
p_load(sf,sp,raster,mgcv,plotrix,gstat,tidyverse,terra,geodata,predicts, RColorBrewer)



##############################################################
#
# DATA DOWNLOAD ----
#
##############################################################

#Assign CRS objects
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
albers <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"

# Download bpen & bimp occurence dataset
df1 <- read.csv("https://raw.githubusercontent.com/abraham-arbelaez/SDMs-bee-project/main/BBOccurenceRecords_Verified_KS.csv")
df1.bpen <- df1 %>%
  filter(YearObserved %in% c("2020", "2021", "2022", "2023")) %>%
  filter(ScientificName == "Bombus pensylvanicus (De Geer, 1773)") %>%
  dplyr::select(ID, YearObserved, Latitude, Longitude)# Keep only the data on current bpen records (2020-2023)

df1.bimp <- df1 %>%
  filter(YearObserved %in% c("2020", "2021", "2022", "2023")) %>%
  filter(ScientificName == "Bombus impatiens Cresson, 1863") %>%
  dplyr::select(ID, YearObserved, Latitude, Longitude)

# Download shapefile of Kansas from census.gov
#download.file("http://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_us_state_20m.zip", destfile = "states.zip")
#unzip("states.zip")
sf.us <- st_read("cb_2015_us_state_20m.shp",quiet = TRUE)
sf.kansas <- sf.us[48,6] # filter to KS
sf.kansas <- as(sf.kansas, 'Spatial') #reset as sf
plot(sf.kansas,main="",col="white")

# Download KS boundary alternative
#sf.ks.boundary <- st_read("KS_boundary.shp")
#sf.ks.boundary <- st_transform(sf.ks.boundary, crs = wgs84)
#sf.ks.boundary <- sf.ks.boundary[,7] #geometry column only

# Download PRISM climate data (monthly, annual values, 4km resolution)
precip <- raster("prism/PRISM_ppt_stable_4kmM3_2022_bil/PRISM_ppt_stable_4kmM3_2022_bil.bil")
tmax <- raster("prism/PRISM_tmax_stable_4kmM3_2022_bil/PRISM_tmax_stable_4kmM3_2022_bil.bil")
tmin <- raster("./prism/PRISM_tmin_stable_4kmM3_2022_bil/PRISM_tmin_stable_4kmM3_2022_bil.bil")
tmean <- raster("./prism/PRISM_tmean_stable_4kmM3_2022_bil/PRISM_tmean_stable_4kmM3_2022_bil.bil")

precip <- projectRaster(precip, crs = wgs84, method = "bilinear")
tmax <- projectRaster(tmax, crs = wgs84, method = "bilinear")
tmin <- projectRaster(tmin, crs = wgs84, method = "bilinear")
tmean <- projectRaster(tmean, crs = wgs84, method = "bilinear")

# Download land cover raster for 2022 from Cropscape
rl.nlcd2022 <- raster("CDL_2022_20.tif")
plot(rl.nlcd2022)

rl.nlcd2022 <- projectRaster(rl.nlcd2022, crs = wgs84, method = "ngb")



####################################################################
#
# Make presence / absence df
#
##################################################################

# Make SpatialPoints data frame
#pts.sample.bpen <- data.frame(long = df1.bpen$Longitude,lat = df1.bpen$Latitude)
#coordinates(pts.sample.bpen) =~ long + lat
#proj4string(pts.sample.bpen) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#pts.sample.bimp <- data.frame(long = df1.bimp$Longitude,lat = df1.bimp$Latitude)
#coordinates(pts.sample.bimp) =~ long + lat
#proj4string(pts.sample.bimp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Convert pts to sf, then create df with only X and Y columns, assign 1 for presence in pa column
pts.bpen <- st_as_sf(df1.bpen, coords = c("Longitude", "Latitude"), crs = wgs84)
df.bpen <- st_coordinates(pts.bpen) %>% as.data.frame
head(df.bpen) #X,Y
pres.bpen <- df.bpen #rename for clarity
pres.bpen$pa <- 1 #Create pa column, assign 1 for presence

pts.bimp <- st_as_sf(df1.bimp, coords = c("Longitude", "Latitude"), crs = wgs84)
df.bimp <- st_coordinates(pts.bimp) %>% as.data.frame
head(df.bpen) #X,Y
pres.bimp <- df.bimp #rename for clarity
pres.bimp$pa <- 1 #Create pa column, assign 1 for presence


# Generate df with pseudo-absences (similar number to presence points), assign 0 for absence in pa column
#pts.pseudo <- st_sample(sf.ks.boundary, 150)
#df.pseudo <- st_coordinates(pts.pseudo) %>% as.data.frame
#head(df.pseudo) #X,Y
#abs <- df.pseudo #rename for clarity
#abs$pa <- 0 #create pa column, assing 0 for absence

# Join into a single df with pres & abs
df.bpen.pa <- pres.bpen
df.bimp.pa <- pres.bimp

head(df.bpen.pa) #X,Y,pa
head(df.bimp.pa) #X,Y,pa

# Create sf object to plot points
sf.bpen.pa <- st_as_sf(df.bpen.pa, coords = c("X","Y"), crs = wgs84)
sf.bimp.pa <- st_as_sf(df.bimp.pa, coords = c("X", "Y"), crs = wgs84)

plot(sf.kansas, main = "KS BPEN pres & pseudo-abs locations between 2020-2023")
plot(sf.bpen.pa, add = T, pch = 17)

plot(sf.kansas, main = "KS BIMP pres & pseudo-abs locations between 2020-2023")
plot(sf.bimp.pa, add = T, pch = 19)




#####################################################################
#
# DATA WRANGLING: Land Cover ----
#
######################################################################


# Make SpatialPoints data frame
pts.sample.bpen <- data.frame(long = df.bpen.pa$X,lat = df.bpen.pa$Y)
coordinates(pts.sample.bpen) =~ long + lat
proj4string(pts.sample.bpen) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

pts.sample.bimp <- data.frame(long = df.bimp.pa$X,lat = df.bimp.pa$Y)
coordinates(pts.sample.bimp) =~ long + lat
proj4string(pts.sample.bimp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


# Assign 1 to grassland/pasture land cover categories
rl.nlcd.grass <- rl.nlcd2022
rl.nlcd.grass[] <- ifelse(rl.nlcd.grass[]==176,1,0) #0 if not grassland (water = 95)

plot(rl.nlcd.grass)
plot(sf.kansas, add= T, col = NA, lwd = 4)

df.bpen.pa$grass.perc <- unlist(lapply(raster::extract(rl.nlcd.grass,pts.sample.bpen,buffer=2000),mean))*100 #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$grass.perc <- unlist(lapply(raster::extract(rl.nlcd.grass,pts.sample.bimp,buffer=2000), mean)) *100



# Assign 1 to developed (all intensities) land cover categories
rl.nlcd.dev <- rl.nlcd2022
rl.nlcd.dev[] <- ifelse(rl.nlcd.dev[] %in% c(82, 121:124), 1,0)

plot(rl.nlcd.dev)
plot(sf.kansas, add= T, col = NA, lwd = 4)

df.bpen.pa$dev.perc <- unlist(lapply(raster::extract(rl.nlcd.dev,pts.sample.bpen,buffer=2000),mean))*100 #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$dev.perc <- unlist(lapply(raster::extract(rl.nlcd.dev,pts.sample.bimp,buffer=2000), mean)) *100



# Assign 1 to water/wetland land categories
rl.nlcd.wet <- rl.nlcd2022
rl.nlcd.wet[] <- ifelse(rl.nlcd.wet[] %in% c(83, 87, 92, 111, 129, 190), 1,0)

plot(rl.nlcd.wet)
plot(sf.kansas, add= T, col = NA, lwd = 4)

df.bpen.pa$wet.perc <- unlist(lapply(raster::extract(rl.nlcd.wet,pts.sample.bpen,buffer=2000),mean))*100 #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$wet.perc <- unlist(lapply(raster::extract(rl.nlcd.wet,pts.sample.bimp,buffer=2000), mean)) *100



# Assign 1 to forest/shrubland land categories
rl.nlcd.wood <- rl.nlcd2022
rl.nlcd.wood[] <- ifelse(rl.nlcd.wood[] %in% c(63:64,141:143, 152), 1,0)

plot(rl.nlcd.wood)
plot(sf.kansas, add= T, col = NA, lwd = 4)

df.bpen.pa$wood.perc <- unlist(lapply(raster::extract(rl.nlcd.wood,pts.sample.bpen,buffer=2000),mean))*100 #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$wood.perc <- unlist(lapply(raster::extract(rl.nlcd.wood,pts.sample.bimp,buffer=2000), mean)) *100




# Assign 1 to agriculture land categories
rl.nlcd.ag <- rl.nlcd2022
rl.nlcd.ag[] <- ifelse(rl.nlcd.ag[] %in% c(1:60,66:77,204:254), 1,0)

plot(rl.nlcd.ag)
plot(sf.kansas, add= T, col = NA, lwd = 4)

df.bpen.pa$ag.perc <- unlist(lapply(raster::extract(rl.nlcd.ag,pts.sample.bpen,buffer=2000),mean))*100 #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$ag.perc <- unlist(lapply(raster::extract(rl.nlcd.ag, pts.sample.bimp,buffer=2000), mean)) *100


#Inspect land cover histograms (NOTE: includes pres & abs)
hist(df.bpen.pa$grass.perc,col="grey", xlab = "% grassland", main="% grassland within \n2 km at BPEN P/A locations")
hist(df.bimp.pa$grass.perc,col="grey",xlab = "% grassland", main="% grassland within \n2 km at BIMP P/A locations")
hist(df.bpen.pa$dev.perc,col="grey",xlab = "% developed", main="% developed within \n2 km at BPEN P/A locations")
hist(df.bimp.pa$dev.perc,col="grey", xlab ="% developed", main="% developed within \n2 km at BIMP P/A locations")
hist(df.bpen.pa$wet.perc,col="grey",xlab = "% water/wetland", main="% water/wetland within \n2 km at BPEN P/A locations")
hist(df.bimp.pa$wet.perc,col="grey",xlab = "% water/wetland", main="% water/wetland within \n2 km at BIMP P/A locations")
hist(df.bpen.pa$wood.perc,col="grey",xlab = "% forest/shrubland", main="% forest/shrubland within \n2 km at BPEN P/A locations")
hist(df.bimp.pa$wood.perc,col="grey",xlab = "% forest/shrubland", main="% forest/shrubland \n2 km at BIMP P/A locations")
hist(df.bpen.pa$ag.perc,col="grey",xlab = "% agriculture", main ="% agriculture within \n2 km at BPEN P/A locations")
hist(df.bimp.pa$ag.perc,col="grey",xlab = "% agriculture", main="% agriculture within \n2 km at BIMP P/A locations")


#Inspect avg land cover percentage in pres vs abs
df.bpen.pa %>%
  group_by(pa) %>%
  summarize(avg.grass = mean(grass.perc),
            avg.dev = mean(dev.perc),
            avg.wet = mean(wet.perc),
            avg.wood = mean(wood.perc),
            avg.ag = mean(ag.perc)) # Abs = 42% grass & 40% ag; Pres = 23% grass & 44% developed...

df.bimp.pa %>%
  group_by(pa) %>%
  summarize(avg.grass = mean(grass.perc),
            avg.dev = mean(dev.perc),
            avg.wet = mean(wet.perc),
            avg.wood = mean(wood.perc),
            avg.ag = mean(ag.perc)) # Abs = 42% grass & 40% ag, Pres = 15% grass & 63% dev



head(df.bpen.pa) # X, Y, pa, grass.perc, dev.perc, wet.perc, wood.perc, ag.perc
head(df.bimp.pa) # X, Y, pa, grass.perc, dev.perc, wet.perc, wood.perc, ag.perc

################################################################################################
#
# Repeat for climate data (mean only) ----
#
################################################################################################

# Precip
df.bpen.pa$precip <- unlist(lapply(raster::extract(precip,pts.sample.bpen,buffer=2000),mean)) #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$precip <- unlist(lapply(raster::extract(precip, pts.sample.bimp,buffer=2000), mean))

# Maximum temp
df.bpen.pa$tmax <- unlist(lapply(raster::extract(tmax,pts.sample.bpen,buffer=2000),mean)) #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$tmax <- unlist(lapply(raster::extract(tmax, pts.sample.bimp,buffer=2000), mean))

# Minimum temp
df.bpen.pa$tmin <- unlist(lapply(raster::extract(tmin,pts.sample.bpen,buffer=2000),mean)) #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$tmin <- unlist(lapply(raster::extract(tmin, pts.sample.bimp,buffer=2000), mean))

# Annual mean temp
df.bpen.pa$tmean <- unlist(lapply(raster::extract(tmean,pts.sample.bpen,buffer=2000),mean)) #2km buffer around survey location, use maximum likelihood to estimate buffer parameter if you are concerned
df.bimp.pa$tmean <- unlist(lapply(raster::extract(tmean, pts.sample.bimp,buffer=2000), mean))



#Inspect climate data for pres vs abs
df.bpen.pa %>%
  group_by(pa) %>%
  summarize(avg.precip = mean(precip),
            avg.tmax = mean(tmax),
            avg.tmin = mean(tmin),
            avg.tmean = mean(tmean))

df.bimp.pa %>%
  group_by(pa) %>%
  summarize(avg.precip = mean(precip),
            avg.tmax = mean(tmax),
            avg.tmin = mean(tmin),
            avg.tmean = mean(tmean))


#Inspect updated df
head(df.bpen.pa) # X, Y, pa, grass.perc, dev.perc, wet.perc, wood.perc, ag.perc, precip, tmax, tmin, tmean
head(df.bimp.pa) # X, Y, pa, grass.perc, dev.perc, wet.perc, wood.perc, ag.perc, precip, tmax, tmin, tmean

df.bpen.final <- df.bpen.pa
df.bpen.final$grass.perc <- df.bpen.final$grass.perc*.01 
df.bpen.final$dev.perc <- df.bpen.final$dev.perc*.01 
df.bpen.final$wet.perc <- df.bpen.final$wet.perc*.01 
df.bpen.final$wood.perc <- df.bpen.final$wood.perc*.01 
df.bpen.final$ag.perc <- df.bpen.final$ag.perc*.01 

df.bimp.final <- df.bimp.pa
df.bimp.final$grass.perc <- df.bimp.final$grass.perc*.01 
df.bimp.final$dev.perc <- df.bimp.final$dev.perc*.01 
df.bimp.final$wet.perc <- df.bimp.final$wet.perc*.01 
df.bimp.final$wood.perc <- df.bimp.final$wood.perc*.01 
df.bimp.final$ag.perc <- df.bimp.final$ag.perc*.01 

#write.csv(df.bimp.final, file = "bimp.csv")
#write.csv(df.bpen.final, file = "bpen.csv")

#### FOLLOWING PACKAGE

library(ppmData)

?ppmData

library(ppmData)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path, pattern='*.tif',full.names = TRUE)
covariates <- rast(lst)
presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
npoints <- 1000
ppmdata1 <- ppmData(npoints = npoints, presences=presences,
                    window = covariates[[1]], covariates=covariates)


plot(ppmdata1)

ppp <- ppmdata1$ppmData
form  <- presence/weights ~ poly(annual_mean_precip, degree = 2) + 
  poly(annual_mean_temp, degree = 2) + 
  poly(distance_from_main_roads, degree = 2)

ft.ppm <- glm(formula = form, data = ppp,
              weights = as.numeric(ppp$weights),
              family = poisson())

pred <- predict(covariates,
                ft.ppm,
                type="response")
plot(pred*prod(res(pred))) # scale response by area of cell.
points(presences[!duplicated(presences),1:2],pch=16,cex=0.75,col=gray(0.2,0.75))


bimp <- df.bimp.final
names(bimp)[3] <- paste("SpeciesID") 
bimp$SpeciesID[1:171] <- "Bombus Pennsylvanicus"

window.ks <- rast(rl.nlcd2022)
class(window.ks)

ag.ks <- rast(rl.nlcd.ag)
dev.ks <- rast(rl.nlcd.dev)
grass.ks <- rast(rl.nlcd.grass)
wet.ks <- rast(rl.nlcd.wet)
wood.ks <- rast(rl.nlcd.wood)




lst <- list()

library(ppmData)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path, pattern='*.tif',full.names = TRUE)
covariates <- df.bimp.final[,4:12]
presences <- df.bimp.final[,c(1:3)]
npoints <- 1000
ppmdata1 <- ppmData(npoints = npoints, presences=presences,
                    window = window.ks, covariates=covariates)


plot(ppmdata1)

ppp <- ppmdata1$ppmData
form  <- presence/weights ~ poly(annual_mean_precip, degree = 2) + 
  poly(annual_mean_temp, degree = 2) + 
  poly(distance_from_main_roads, degree = 2)

ft.ppm <- glm(formula = form, data = ppp,
              weights = as.numeric(ppp$weights),
              family = poisson())

pred <- predict(covariates,
                ft.ppm,
                type="response")
plot(pred*prod(res(pred))) # scale response by area of cell.
points(presences[!duplicated(presences),1:2],pch=16,cex=0.75,col=gray(0.2,0.75))








































































## Notes from discussion w/ Trevor, Nikki, and Yoshi: 

#point process uses quadrature pts not pseudo absence pts

# Poisson point process solves pseudo-absence...paper 

# Include a clump of covariates that influence probability of presence AND a separate clump
#  that influence probability of reporting (e.g. distance to major cities)

# Model for presence-only, model for pa, model for pooled... compare differences
# Keep separate, then analyze slopes (one is high quality, one is low quality, do they show the same info??)
# Hill 1965 paper- number of datasets to include
