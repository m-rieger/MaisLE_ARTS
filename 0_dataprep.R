#### data preparation ####

## packages 
library(tidyverse)
library(sf)

crs <- 31467 # Gauss Krüger in m
crsPlot <- 3857 # 
crsLL <- 4326 # lon lat in degree

#### 1) data for stations ------------------------------------------------------
df.stat <- read.csv(here("data", "station.csv"))

## reduce to stations included in analysis; remove columns
df.stat <- df.stat[df.stat$analysis == "yes",]
df.stat <- dplyr::select(df.stat, -c(analysis, station.hostname, station.wireguard_peer.id, 
                              station.wireguard_peer.address, station.shelter_id))

## get stations from analysis
df.stat <- st_as_sf(df.stat, coords = c("station.lon", "station.lat"), crs = crsLL)

df.stat$type <- "omni"
df.stat$type[df.stat$station.type == "VHF Directional Station"] <- "direct"

## save as gpgk file
st_write(df.stat, here("data", "Station.gpkg"), append = FALSE, 
         layer_options = "GEOMETRY_NAME=geometry")

#### 2) gpkg files for testtracks ----------------------------------------------
myfiles <- list.files(path = "../data/data_groundtruth", pattern=".+\\.gpx", full.names=F)

shp.p <- NULL
shp.l <- NULL

## loop through files
for(i in myfiles) {
  ## points
  gpx <- st_read(paste0("../data/data_groundtruth/", i), layer = "track_points")

  gpx <- gpx[, c("time", "ele")]
  gpx$track <- i

  shp.p <- rbind(shp.p, gpx)

  ## lines (only for length)
  gpxl <- st_read(paste0("../data/data_groundtruth/", i), layer = "tracks")

  gpxl <- gpxl[, c("name")]
  gpxl$track <- i

  shp.l <- rbind(shp.l, gpxl)

}

## separate filename to get info
shp.p <- separate(shp.p, track,
                     sep = " - ",
                     into = c("time_start", "area2", "stat", "comment"))

shp.p <- separate(shp.p, area2,
                     sep = " ",
                     into = c("area", "pos2"))
shp.p$date   <- as.Date(shp.p$time)
shp.p$site <- paste0("mais", shp.p$area)

## separate filename to get info
shp.l <- separate(shp.l, track,
                  sep = " - ",
                  into = c("time_start", "area2", "stat", "comment"))
shp.l <- separate(shp.l, time_start,
                  sep = "-",
                  into = c("date", "time"))
shp.l <- separate(shp.l, area2,
                  sep = " ",
                  into = c("area", "pos2"))
shp.l$date   <- as.POSIXct(shp.l$date, tryFormats = c("%Y%m%d"), tz = "CET")
shp.l$site <- paste0("mais", shp.l$area)

## save as gpkg
st_write(shp.p, here("data", "Testtracks_points.gpkg"), append = F, 
         layer_options = "GEOMETRY_NAME=geometry")
st_write(shp.l, here("data", "Testtracks_lines.gpkg"), append = F, 
         layer_options = "GEOMETRY_NAME=geometry")
