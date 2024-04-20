library(sf)
library(tidyverse)
library(spData)

ks <- st_read("KS_counties.shp") #WGS84

bba <- read.csv("BBW_KansasRecords_5April2024_Holthaus_KState.csv")

bba.bpen <- bba %>%
  filter(vernacularName == "American bumble bee" )

bpen_current <- read.csv(file = "ver_KS_bpen_current.csv")

sf.bba.bpen <- st_as_sf(bba.bpen, coords = c("decimalLongitude", "decimalLatitude"), crs = '+proj=longlat +datum=WGS84 +no_defs +type=crs')
sf.bpen_current <- st_as_sf(bpen_current, coords = c("Longitude", "Latitude"), crs = '+proj=longlat +datum=WGS84 +no_defs +type=crs')

ggplot() +
  geom_sf(data = ks) +
  geom_sf(data = sf.bpen_current, color = "black") +
  geom_sf(data = sf.bba, color = "tomato") +
  ggtitle("Current BPEN records from GPBBA (red) and GBIF (black)")

