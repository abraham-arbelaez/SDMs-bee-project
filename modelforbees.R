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



# tweaking it with mgcv
library(mgcv)
form  <- presence/weights ~ s(annual_mean_precip, bs = "ts") + 
  s(annual_mean_temp, bs = "ts") + 
  s(distance_from_main_roads, bs = "ts")

ft.ppm <- gam(formula = form, data = ppp,
              weights = as.numeric(ppp$weights),
              family = poisson())

pred <- predict(covariates,
                ft.ppm,
                type="response")
plot(pred*prod(res(pred))) # scale response by area of cell.
points(presences[!duplicated(presences),1:2],pch=16,cex=0.75,col=gray(0.2,0.75))




# BEES

sf.us <- st_read("cb_2015_us_state_20m.shp",quiet = TRUE)
sf.kansas <- sf.us[48,6] # filter to KS

# ALL COVARIATES

precips <- rast(precip)
tmaxs <- rast(tmax)
tmins <- rast(tmin)
tmeans <- rast(tmean)
rl.nlcd <- rast(rl.nlcd2022)

precips <- rast(precip, crs(sf.kansas))
precips <- project(precips, crs(sf.kansas))
precips <- mask(precips, sf.kansas)
precips <- crop(precips, sf.kansas)
plot(precips)

tmaxs <- rast(tmax, crs(sf.kansas))
tmaxs <- project(tmaxs, crs(sf.kansas))
tmaxs <- mask(tmaxs, sf.kansas)
tmaxs <- crop(tmaxs, sf.kansas)
plot(tmaxs)

tmins <- rast(tmin, crs(sf.kansas))
tmins <- project(tmins, crs(sf.kansas))
tmins <- mask(tmins, sf.kansas)
tmins <- crop(tmins, sf.kansas)
plot(tmins)

tmeans <- rast(tmean, crs(sf.kansas))
tmeans <- project(tmeans, crs(sf.kansas))
tmeans <- mask(tmeans, sf.kansas)
tmeans <- crop(tmeans, sf.kansas)
plot(tmeans)

ag <- rast(rl.nlcd.ag, crs(sf.kansas))
ag <- project(ag, crs(sf.kansas))
ag <- mask(ag, sf.kansas)
ag <- crop(ag, sf.kansas)
plot(ag)

dev <- rast(rl.nlcd.dev, crs(sf.kansas))
dev <- project(dev, crs(sf.kansas))
dev <- mask(dev, sf.kansas)
dev <- crop(dev, sf.kansas)
plot(dev)

grass <- rast(rl.nlcd.grass, crs(sf.kansas))
grass <- project(grass, crs(sf.kansas))
grass <- mask(grass, sf.kansas)
grass <- crop(grass, sf.kansas)
plot(grass)

wet <- rast(rl.nlcd.wet, crs(sf.kansas))
wet <- project(wet, crs(sf.kansas))
wet <- mask(wet, sf.kansas)
wet <- crop(wet, sf.kansas)
plot(wet)

wood <- rast(rl.nlcd.wood, crs(sf.kansas))
wood <- project(wood, crs(sf.kansas))
wood <- mask(wood, sf.kansas)
wood <- crop(wood, sf.kansas)
plot(wood)

# resampling
# high res
precips <- resample(precips, ag)
tmaxs <- resample(tmaxs, ag)
tmins <- resample(tmins, ag)
tmeans <- resample(tmeans, ag)

# low res
rag <- resample(ag, tmeans)
rdev <- resample(dev, tmeans)
rgrass <- resample(grass, tmeans)
rwet <- resample(wet, tmeans)
rwood <- resample(wood, tmeans)

lst <- list(precips, tmaxs, tmins, tmeans, ag, dev, grass, wet, wood)
covariates <- rast(lst)

#sf.kansas <- as(sf.kansas, 'Spatial') #reset as sf
#sf.kansas <- as(sf.kansas, 'SpatRaster')

#window.ks <- rast(sf.kansas)
#ext(window.ks) <- ext()
#class(window.ks)

bimp <- df.bimp.final
names(bimp)[3] <- paste("SpeciesID") 
bimp$SpeciesID[1:171] <- "Bombus Impatiens"

bpen <- df.bpen.final
names(bpen)[3] <- paste("SpeciesID") 
bpen$SpeciesID[1:157] <- "Bombus Pennsylvanicus"

bees <- rbind(bimp, bpen)

# Pennsylvanicus

presences <- subset(bees,SpeciesID %in% "Bombus Pennsylvanicus")

npoints <- 100
ppmdata1 <- ppmData(npoints = npoints, presences=presences,
                    window = covariates[[1]], covariates=covariates)


#plot(ppmdata1)

ppp <- ppmdata1$ppmData

colnames(ppp) <- c("X", "Y", "precip", "tmax", "tmin", "tmean", "ag", "dev",
                   "grass", "wet", "wood", "presence", "weights")

form  <- presence/weights ~ s(precip, bs = "ts") + 
  s(tmax, bs = "ts") + 
  s(tmin, bs = "ts") +
  s(tmean, bs = "ts") +
  s(ag, bs = "ts") +
  s(dev, bs = "ts") +
  s(grass, bs = "ts") +
  s(wet, bs = "ts") +
  s(wood, bs = "ts")

ft.ppm <- gam(formula = form, data = ppp,
              weights = as.numeric(ppp$weights),
              family = poisson())

summary(ft.ppm)

names(covariates) <- c("precip", "tmax", "tmin", "tmean", "ag", "dev",
                   "grass", "wet", "wood")

pred <- predict(covariates,
                ft.ppm,
                type="response")

plot(pred*prod(res(pred))) # scale response by area of cell.
points(presences[!duplicated(presences),1:2],pch=16,cex=0.25,col="black")

#
## ggplot
#


library(viridis)

pred1_df <- as.data.frame(pred*prod(res(pred)), xy = TRUE)

pennsylvanicus_final <- ggplot() +
  geom_raster(data = pred1_df, aes(x = x, y = y, fill = lyr1)) +
  geom_point(data = bpen, aes(x = X, y = Y), color = "darkgreen", size = 0.5)+
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(x = "Longitude", y = "Latitude", fill = "Intensity") +
  theme_minimal() +
  ggtitle("Bombus Pennsylvanicus")


### Impatiens

presences <- subset(bees,SpeciesID %in% "Bombus Impatiens")

npoints <- 100000
ppmdata1 <- ppmData(npoints = npoints, presences=presences,
                    window = covariates[[1]], covariates=covariates)


#plot(ppmdata1)

ppp <- ppmdata1$ppmData

colnames(ppp) <- c("X", "Y", "precip", "tmax", "tmin", "tmean", "ag", "dev",
                   "grass", "wet", "wood", "presence", "weights")

form  <- presence/weights ~ s(precip, bs = "ts") + 
  s(tmax, bs = "ts") + 
  s(tmin, bs = "ts") +
  s(ag, bs = "ts") +
  s(dev, bs = "ts") +
  s(grass, bs = "ts") +
  s(wet, bs = "ts") +
  s(wood, bs = "ts")



form  <- presence/weights ~ s(precip, bs = "ts") + 
  s(tmin, bs = "ts") +
  s(dev, bs = "ts") +
  s(wet, bs = "ts") +
  s(wood, bs = "ts")


ft.ppm <- gam(formula = form, data = ppp,
              weights = as.numeric(ppp$weights),
              family = poisson())

summary(ft.ppm)

names(covariates) <- c("precip", "tmax", "tmin", "tmean", "ag", "dev",
                       "grass", "wet", "wood")

pred <- predict(covariates,
                ft.ppm,
                type="response")

plot(pred*prod(res(pred))) # scale response by area of cell.
points(presences[!duplicated(presences),1:2],pch=16,cex=0.75,col=gray(0.2,0.75))


pred2_df <- as.data.frame(pred*prod(res(pred)), xy = TRUE)

impatiens_final <- ggplot() +
  geom_raster(data = pred2_df, aes(x = x, y = y, fill = lyr1)) +
  geom_point(data = bimp, aes(x = X, y = Y), color = "darkgreen", size = 0.5)+
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(x = "Longitude", y = "Latitude", fill = "Intensity") +
  theme_minimal() +
  ggtitle("Bombus Impatiens")

library(patchwork)
pennsylvanicus_final + impatiens_final
