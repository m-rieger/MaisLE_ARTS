#### merge and summarize testtag positions for calibration ####

## ToDO
# integrate mean files to check whether distance error gets better

rm(list = ls())

## packages
library(ggdist)
library(tidyverse)
library(sf)
library(data.table)
library(here)
library(geosphere)
#library(spatstat) # density distribution for sf objects
library(raster) # density as raster
library(stars)
library(terra)
#library(zoo)

bs <- 14

## define variables (here, you might want to change something)
# chose data type (animal data or testtag data)
dtyp <- "cali" # "animal", "testtag" or "cali"
# chose time interval in Jupyter data (usually 2 sec)
t <- "2 sec" # "x sec"
# minimum number of stations per position
NStat <- 1
# minimum number of antennas per position
NAnt <- 1
# write gpkg for rolling mean and filtered data
gpkg <- TRUE

## additional variables
t2 <- as.numeric(gsub(" sec", "", t))

#### loop through gpkg-files ####

## get all files and create folder for saved files
df.L <- list.files(paste0("../data/data_", dtyp, "/"), pattern = ".gpkg")
if(!dir.exists(paste0("../data/data_", dtyp, "/savedFiles"))) dir.create(paste0("../data/data_", dtyp, "/savedFiles"))

## exclude -mean files
df.Lm <- df.L[grep("-mean", df.L)]
df.L <- df.L[-grep("-mean", df.L)]

df.L1 <- df.L[grep("antennabeams.gpkg", df.L)]
df.L2 <- df.L[grep("antennabeam_s", df.L)]
df.L3 <- df.L[grep("multilateration", df.L)]

df.Lm1 <- df.Lm[grep("antennabeams-mean.gpkg", df.Lm)]
df.Lm2 <- df.Lm[grep("antennabeam_s-mean", df.Lm)]
df.Lm3 <- df.Lm[grep("multilateration-mean", df.Lm)]

#### normal data ####
# antennabeam quadrologger
df1 <- data.frame()

for(i in df.L1) {
  tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  
  ## transform crs to lon lat and save crs
  tmp <- st_transform(tmp, crs = 4326)
  crs <- st_crs(tmp)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## get coordinates from geometry
  tmp <- tmp %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    st_drop_geometry()
  
  ## add cols
  tmp$detR <- strsplit(i, "_")[[1]][2]
  tmp$ant <- strsplit(i, "_")[[1]][3]
  tmp$Distance..m..min <- NA
  tmp$meth <- "ab_ql"
  
  df1 <- rbind(df1, tmp)
  
}

# antennabeam monologger
df2 <- data.frame()

for(i in df.L2) {
  tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  
  ## transform crs to lon lat and save crs
  tmp <- st_transform(tmp, crs = 4326)
  crs <- st_crs(tmp)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## get coordinates from geometry
  tmp <- tmp %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    st_drop_geometry()
  
  ## add cols
  tmp$detR <- "no"
  tmp$Distance..m..min <- NA
  tmp$ant <- "1-10"
  tmp$meth <- "ab_ml"
  
  df2 <- rbind(df2, tmp)
  
}

# multilateration monologger
df3 <- data.frame()

for(i in df.L3) {
  tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  
  ## transform crs to lon lat and save crs
  tmp <- st_transform(tmp, crs = 4326)
  crs <- st_crs(tmp)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## get coordinates from geometry
  tmp <- tmp %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    st_drop_geometry()
  
  ## add cols
  tmp$detR <- "no"
  tmp$Weight <- NA
  tmp$ant <- "1-10"
  tmp$Antenna.Count <- tmp$Station.Count
  tmp$meth <- "ml_ml"
  
  df3 <- rbind(df3, tmp)
  
}

df <- rbind(df1, df2)
df <- rbind(df, df3)
rm(df1); rm(df2); rm(df3)

#### mean data ####
# antennabeam quadrologger
df1 <- data.frame()

for(i in df.Lm1) {
  tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  
  ## transform crs to lon lat and save crs
  tmp <- st_transform(tmp, crs = 4326)
  crs <- st_crs(tmp)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## get coordinates from geometry
  tmp <- tmp %>%
    dplyr::mutate(lon.m = sf::st_coordinates(.)[,1],
                  lat.m = sf::st_coordinates(.)[,2]) %>%
    st_drop_geometry()
  
  ## add cols
  tmp$detR <- strsplit(i, "_")[[1]][2]
  tmp$ant <- strsplit(i, "_")[[1]][3]
  tmp$Distance..m..min <- NA
  tmp$meth <- "ab_ql"
  
  df1 <- rbind(df1, tmp)
  
}

# antennabeam monologger
df2 <- data.frame()

for(i in df.Lm2) {
  tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  
  ## transform crs to lon lat and save crs
  tmp <- st_transform(tmp, crs = 4326)
  crs <- st_crs(tmp)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## get coordinates from geometry
  tmp <- tmp %>%
    dplyr::mutate(lon.m = sf::st_coordinates(.)[,1],
                  lat.m = sf::st_coordinates(.)[,2]) %>%
    st_drop_geometry()
  
  ## add cols
  tmp$detR <- "no"
  tmp$Distance..m..min <- NA
  tmp$ant <- "1-10"
  tmp$meth <- "ab_ml"
  
  df2 <- rbind(df2, tmp)
  
}

# multilateration monologger
df3 <- data.frame()

for(i in df.Lm3) {
  tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  
  ## transform crs to lon lat and save crs
  tmp <- st_transform(tmp, crs = 4326)
  crs <- st_crs(tmp)
  
  ## get Individual and Transmitter
  ind <- unique(tmp$Individual)
  trans <- unique(tmp$Transmitter)
  
  ## get coordinates from geometry
  tmp <- tmp %>%
    dplyr::mutate(lon.m = sf::st_coordinates(.)[,1],
                  lat.m = sf::st_coordinates(.)[,2]) %>%
    st_drop_geometry()
  
  ## add cols
  tmp$detR <- "no"
  tmp$Weight <- NA
  tmp$ant <- "1-10"
  tmp$Antenna.Count <- tmp$Station.Count
  tmp$meth <- "ml_ml"
  
  df3 <- rbind(df3, tmp)
  
}

dfm <- rbind(df1, df2)
dfm <- rbind(dfm, df3)
rm(df1); rm(df2); rm(df3)


#### merge & add GT data ####
df <- left_join(unique(df), unique(dfm[, c("lon.m", "lat.m", "Individual", "X_time", "ant", "detR", "meth")]),
                by = c("Individual", "X_time", "ant", "detR", "meth"))

df <- select(df, -c("Planner", "Transmitter", "Transmitter.Id"))

## read groundtruth data to estimate position error
df.GPS  <- fread(here("..", "2023", "00_Telemetry_2023", "data_groundtruth", "01_GTdata-processed", "GTdata.csv"))

## change colnames
colnames(df.GPS)[colnames(df.GPS) == "lon"] <- "lon.true"
colnames(df.GPS)[colnames(df.GPS) == "lat"] <- "lat.true"
colnames(df.GPS)[colnames(df.GPS) == "Testtag"] <- "Individual"
colnames(df.GPS)[colnames(df.GPS) == "date_start"] <- "X_time"

## get moving data
df.GPS <- df.GPS[df.GPS$pos2 %in% c("grid", "Testtrack")]

## merge with df by transmitter and time
df <- left_join(df, 
                df.GPS[, c("Individual", "X_time", "lat.true", "lon.true", "pos1", "pos2")],
                by = c("Individual", "X_time"))
df <- unique(df)

## get distance error
df <- df[!is.na(df$lat.true),] %>% rowwise %>%
  mutate(distance = distm(x = c(lon, lat),
                          y = c(lon.true, lat.true)),
         distance.m = distm(x = c(lon.m, lat.m),
                          y = c(lon.true, lat.true))) %>%
  ungroup()


## write df
df <- st_as_sf(df[!is.na(df$lon) & !is.na(df$lat),], coords = c("lon", "lat"), crs = crs)
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raw.gpkg")
st_write(df, layer = 'unfiltered', append = F, dsn = dsn)

## read df
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raw.gpkg")
df <- st_read(layer = 'unfiltered', dsn = dsn)

## run markdown to visualize output
rmarkdown::render(input = here("R_pre", "plotCali.Rmd"), 
                  output_format = "bookdown::html_document2")


#### density ####
tmp <- st_read(here("data", "station.csv"))
df.p <- st_as_sf(tmp[tmp$station.id %in% c("c1l1", "c2l1", "c3l1", "c4l1", "c5l1", "c6l1", "c7l1", "c8l1", "c9l1", "c10l1"),], coords = c("station.lon", "station.lat"), crs = 4326)
df.p <- st_as_sf(tmp[tmp$station.id %in% c("c1l1", "c2l1"),], coords = c("station.lon", "station.lat"), crs = 4326)
df.p <- st_as_sf(tmp[tmp$station.id == "c1l1",], coords = c("station.lon", "station.lat"), crs = 4326)

crs <- 31467 # Gauss-Krüger (in m)
# 
# for(s in c("c1l1", "c2l1", "c3l1", "c4l1", "c5l1", "c6l1", "c7l1", "c8l1", "c9l1", "c10l1")) {
#   df.p <- st_as_sf(tmp[tmp$station.id == s,], coords = c("station.lon", "station.lat"), crs = 4326)
#   df.p <- st_transform(df.p, crs = crs)
#   b100 <- st_buffer(df.p, dist = 100)
#   b100$dens <- 1
#   b200 <- st_buffer(b100, dist = 100)
#   b200$dens <- 0.75
#   b300 <- st_buffer(b200, dist = 100)
#   b300$dens <- 0.5
#   b400 <- st_buffer(b300, dist = 100)
#   b400$dens <- 0.25
#   # shp.buffer <- rbind(b100, b200)
#   # shp.buffer <- rbind(shp.buffer, b300)
#   # shp.buffer <- rbind(shp.buffer, b400)
#   
#   df.rast1 <- st_rasterize(b100, st_as_stars(st_bbox(b400), nx = 40, ny = 40, values = b100$dens))
#   df.rast2 <- st_rasterize(b200, st_as_stars(st_bbox(b400), nx = 40, ny = 40, values = b200$dens))
#   df.rast3 <- st_rasterize(b300, st_as_stars(st_bbox(b400), nx = 40, ny = 40, values = b300$dens))
#   df.rast4 <- st_rasterize(b400, st_as_stars(st_bbox(b400), nx = 40, ny = 40, values = b400$dens))
#   
#   df.rast <- df.rast1 + df.rast2
#   df.rast <- df.rast + df.rast3
#   df.rast <- df.rast + df.rast4
#   
#   plot(df.rast, axes = TRUE)
#   
# }





# df.rast <- st_rasterize(shp.buffer, st_as_stars(st_bbox(shp.buffer), nx = 40, ny = 40, values = NA_real_))
# 
# df.rast$dens <- 1
# plot(df.rast, axes = TRUE)
# 
# df.rast2 <- df.rast + df.rast
# #df.rast <- raster(shp.buffer)
# # df.rast <- raster(shp.buffer)
# #plot(df.rast, breaks = "equal")
# 
# ggplot(shp.buffer) + geom_sf(aes(alpha = dens), fill = "blue")

#########################
stations <- c("c1l1", "c2l1", "c3l1", "c4l1", "c5l1", "c6l1", "c7l1", "c8l1", "c9l1", "c10l1")

b100 <- b200 <- b300 <- b400 <- shp.buffer <- list()
df.r <- NULL

## get raster dimension
df.p <- st_as_sf(tmp[tmp$station.id %in% stations,], coords = c("station.lon", "station.lat"), crs = 4326)
df.p <- st_transform(df.p, crs = crs)
r.dim <- st_buffer(df.p, dist = 500)

# Define extent and resolution
raster_extent <- extent(r.dim)  # Adjust according to your needs
raster_resolution <- 10  # Change resolution as needed

# Create an empty raster
r <- raster(raster_extent, resolution = raster_resolution)

for(i in 1:length(stations)) {
  s <- stations[i]
  
  df.p <- st_as_sf(tmp[tmp$station.id == s,], coords = c("station.lon", "station.lat"), crs = 4326)
  df.p <- st_transform(df.p, crs = crs)
  b100[[s]] <- st_buffer(df.p, dist = 100)
  b200[[s]] <- st_buffer(df.p, dist = 200)
  b300[[s]] <- st_buffer(df.p, dist = 300)
  b400[[s]] <- st_buffer(df.p, dist = 400)
  
  b400[[s]] <- st_difference(b400[[s]], b300[[s]])
  b300[[s]] <- st_difference(b300[[s]], b200[[s]])
  b200[[s]] <- st_difference(b200[[s]], b100[[s]])
  
  b100[[s]]$dens <- 1/length(stations)
  b200[[s]]$dens <- 0.75/length(stations)
  b300[[s]]$dens <- 0.5/length(stations)
  b400[[s]]$dens <- 0.25/length(stations)
  b100[[s]] <- b100[[s]][, "dens"]
  b200[[s]] <- b200[[s]][, "dens"]
  b300[[s]] <- b300[[s]][, "dens"]
  b400[[s]] <- b400[[s]][, "dens"]
  shp.buffer[[s]] <- rbind(b100[[s]], b200[[s]])
  shp.buffer[[s]] <- rbind(shp.buffer[[s]], b300[[s]])
  shp.buffer[[s]] <- rbind(shp.buffer[[s]], b400[[s]])
  
  r.tmp <- raster::rasterize(shp.buffer[[s]], r, field = "dens", fun = sum)
  r.tmp[is.na(r.tmp)] <- 0
  if(s == stations[1]) df.r <- r.tmp
  else df.r <- df.r + r.tmp

  
}

plot(df.r, main = "Rasterized Polygon")
