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


#### 3) merge raw data and remove unnecessary cols -----------------------------
## get all files and create folder for saved files
df.L <- list.files("../data/data_cali/", pattern = ".gpkg")

if(!dir.exists("./data/cali")) {
  dir.create("./data/cali")}

## separate files based on method:
# df.L1, antennabeams = direct antenna beams
# df.L2, antennabeams_s = omni antenna beams
# df.L3, multilateration = omni multilateration
# df.L4, intersections = direct angulation
df.L1 <- df.L[grep("antennabeams.gpkg", df.L)]
df.L2 <- df.L[grep("antennabeam_s", df.L)]
df.L3 <- df.L[grep("multilateration", df.L)]
df.L4 <- df.L[grep("intersections", df.L)]

## direct antenna beams
df1 <- data.frame()

for(i in df.L1) {
  tmp <- st_read(paste0("../data/data_cali/", i))
  
  ## transform crs to lon lat
  tmp <- st_transform(tmp, crs = crsLL)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)

  ## add cols
  tmp$r <- strsplit(i, "_")[[1]][2]
  tmp$ant <- strsplit(i, "_")[[1]][3]
  tmp$Distance..m..min <- NA
  tmp$meth <- "direct.ab"
  
  df1 <- rbind(df1, tmp)
  
}

## omni antenna beams
df2 <- data.frame()

for(i in df.L2) {
  tmp <- st_read(paste0("../data/data_cali/", i))
  
  ## transform crs to lon lat
  tmp <- st_transform(tmp, crs = crsLL)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## add cols
  tmp$r <- "no"
  tmp$Distance..m..min <- NA
  tmp$ant <- "1-10"
  tmp$meth <- "omni.ab"
  
  df2 <- rbind(df2, tmp)
  
}

## omni multilateration
df3 <- data.frame()

for(i in df.L3) {
  tmp <- st_read(paste0("../data/data_cali/", i))
  
  ## transform crs to lon lat
  tmp <- st_transform(tmp, crs = crsLL)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## add cols
  tmp$r <- "no"
  tmp$Weight <- NA
  tmp$ant <- "1-10"
  tmp$Antenna.Count <- tmp$Station.Count
  tmp$meth <- "omni.ml"
  
  df3 <- rbind(df3, tmp)
  
}

## direct angulation
df4 <- data.frame()

for(i in df.L4) {
  tmp <- st_read(paste0("../data/data_cali/", i))
  
  ## transform crs to lon lat
  tmp <- st_transform(tmp, crs = crsLL)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## add cols
  tmp$r <- "no"
  tmp$Weight <- NA
  tmp$Antenna.Count <- 2*tmp$Station.Count
  tmp$ant <- strsplit(i, "_")[[1]][3]
  tmp$meth <- "direct.an"
  
  df4 <- rbind(df4, tmp)
  
}

## combine data in one df
df <- do.call("rbind", list(df1, df2, df3, df4))
rm(df1); rm(df2); rm(df3); rm(df4)

## remove unnecessary IDs
df <- dplyr::select(df, -c("Planner", "Transmitter", "Transmitter.Id", 
                           "ant", "Distance..m..min", "Transmitter.Model"))

## remove "eurofins-" for anonymity
df$Stations <- str_replace_all(df$Stations, "eurofins-maisC-", "")
df$Stations <- str_replace_all(df$Stations, "eurofins-maisD-", "")

## save data as gpkg
dsn <- paste0("./data/cali/Data_cali_raw_maisC.gpkg")
st_write(df[df$Project == "maisC",], layer = 'unfiltered', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")
dsn <- paste0("./data/cali/Data_cali_raw_maisD.gpkg")
st_write(df[df$Project == "maisD",], layer = 'unfiltered', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")
