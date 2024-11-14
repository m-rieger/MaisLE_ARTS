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

## should data be read in new (TRUE) or loaded (FALSE)
READ <- FALSE

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

## stations to be included in stations density
stations <- c("c1l1", "c2l1", "c3l1", "c4l1", "c5l1", "c6l1", "c7l1", "c8l1", "c9l1", "c10l1",
              "d1l1", "d2l1", "d3l1", "d4l1", "d5l1", "d6l1", "d7l1", "d8l1")

## data for stations density
df.stat <- st_read(here("data", "station.csv"))
crs <- 31467 # Gauss-Krüger (in m)

# get raster dimension
df.stat <- st_as_sf(df.stat[df.stat$station.id %in% stations,], coords = c("station.lon", "station.lat"), crs = 4326)
df.stat <- st_transform(df.stat, crs = crs)
r.dim <- st_buffer(df.stat, dist = 800) # 800m buffer

# Define extent and resolution
raster_extent <- extent(r.dim)  # Adjust according to your needs
raster_resolution <- 10  # Change resolution as needed (in meter) for stat density raster
raster_incr <- 5  # Change resolution as needed (in meter) for point summary raster (raster_resoultion*raster_incr)

#### 1) loop through gpkg-files ####
if(READ) {
  
  ## get all files and create folder for saved files
  df.L <- list.files(paste0("../data/data_", dtyp, "/"), pattern = ".gpkg")
  if(!dir.exists(paste0("../data/data_", dtyp, "/savedFiles"))) dir.create(paste0("../data/data_", dtyp, "/savedFiles"))
  
  ## exclude -mean files
  df.Lm <- df.L[grep("-mean", df.L)]
  df.L <- df.L[-grep("-mean", df.L)]
  
  df.L1 <- df.L[grep("antennabeams.gpkg", df.L)]
  df.L2 <- df.L[grep("antennabeam_s", df.L)]
  df.L3 <- df.L[grep("multilateration", df.L)]
  df.L4 <- df.L[grep("intersections", df.L)]
  
  df.Lm1 <- df.Lm[grep("antennabeams-mean.gpkg", df.Lm)]
  df.Lm2 <- df.Lm[grep("antennabeam_s-mean", df.Lm)]
  df.Lm3 <- df.Lm[grep("multilateration-mean", df.Lm)]
  df.Lm4 <- df.L[grep("intersections-mean", df.Lm)]
  
  #### 1.1) normal data ####
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
  
  # intersection quadrologger
  df4 <- data.frame()
  
  for(i in df.L4) {
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
    tmp$Antenna.Count <- 2*tmp$Station.Count
    tmp$ant <- strsplit(i, "_")[[1]][3]
    tmp$meth <- "in_ql"
    
    df4 <- rbind(df4, tmp)
    
  }
  
  df <- rbind(df1, df2)
  df <- rbind(df, df3)
  df <- rbind(df, df4)
  rm(df1); rm(df2); rm(df3); rm(df4)
  
  #### 1.2) mean data ####
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
  
  # intersection quadrologger
  df4 <- data.frame()
  
  for(i in df.L4) {
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
    tmp$Antenna.Count <- 2*tmp$Station.Count
    tmp$ant <- strsplit(i, "_")[[1]][3]
    tmp$meth <- "in_ql"
    
    df4 <- rbind(df4, tmp)
    
  }
  dfm <- rbind(df1, df2)
  dfm <- rbind(dfm, df3)
  dfm <- rbind(dfm, df4)
  rm(df1); rm(df2); rm(df3); rm(df4)
  
  
  #### 1.3) merge & add GT data ####
  df <- left_join(unique(df), unique(dfm[, c("lon.m", "lat.m", "Individual", "X_time", "ant", "detR", "meth")]),
                  by = c("Individual", "X_time", "ant", "detR", "meth"))
  
  df <- dplyr::select(df, -c("Planner", "Transmitter", "Transmitter.Id"))
  
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
                  by = c("Individual", "X_time")) # many-to-many relationship?
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
  
  ## run markdown to visualize output
  rmarkdown::render(input = here("R_pre", "plotCali.Rmd"), 
                    output_format = "bookdown::html_document2")
}



#### 2) stations density ####

## read df
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raw.gpkg")
df <- st_read(layer = 'unfiltered', dsn = dsn)

b100 <- b200 <- b300 <- b400 <- b500 <- b600 <- b700 <- shp.buffer <- list()
df.r <- NULL

# Create an empty raster
r <- raster(raster_extent, resolution = raster_resolution)

for(i in 1:length(stations)) {
  s <- stations[i]
  
  df.p <- st_as_sf(df.stat[df.stat$station.id == s,], coords = c("station.lon", "station.lat"), crs = 4326)
  df.p <- st_transform(df.p, crs = crs)
  b100[[s]] <- st_buffer(df.p, dist = 100)
  b200[[s]] <- st_buffer(df.p, dist = 200)
  b300[[s]] <- st_buffer(df.p, dist = 300)
  b400[[s]] <- st_buffer(df.p, dist = 400)
  b500[[s]] <- st_buffer(df.p, dist = 500)
  b600[[s]] <- st_buffer(df.p, dist = 600)
  b700[[s]] <- st_buffer(df.p, dist = 700)
  
  b700[[s]] <- st_difference(b700[[s]], b600[[s]])
  b600[[s]] <- st_difference(b600[[s]], b500[[s]])
  b500[[s]] <- st_difference(b500[[s]], b400[[s]])
  b400[[s]] <- st_difference(b400[[s]], b300[[s]])
  b300[[s]] <- st_difference(b300[[s]], b200[[s]])
  b200[[s]] <- st_difference(b200[[s]], b100[[s]])
  
  b100[[s]]$dens <- 1
  b200[[s]]$dens <- 0.85
  b300[[s]]$dens <- 0.7
  b400[[s]]$dens <- 0.55
  b500[[s]]$dens <- 0.4
  b600[[s]]$dens <- 0.25
  b700[[s]]$dens <- 0.1
  b100[[s]] <- b100[[s]][, "dens"]
  b200[[s]] <- b200[[s]][, "dens"]
  b300[[s]] <- b300[[s]][, "dens"]
  b400[[s]] <- b400[[s]][, "dens"]
  b500[[s]] <- b500[[s]][, "dens"]
  b600[[s]] <- b600[[s]][, "dens"]
  b700[[s]] <- b700[[s]][, "dens"]
  shp.buffer[[s]] <- rbind(b100[[s]], b200[[s]])
  shp.buffer[[s]] <- rbind(shp.buffer[[s]], b300[[s]])
  shp.buffer[[s]] <- rbind(shp.buffer[[s]], b400[[s]])
  shp.buffer[[s]] <- rbind(shp.buffer[[s]], b500[[s]])
  shp.buffer[[s]] <- rbind(shp.buffer[[s]], b600[[s]])
  shp.buffer[[s]] <- rbind(shp.buffer[[s]], b700[[s]])
  
  r.tmp <- raster::rasterize(shp.buffer[[s]], r, field = "dens", fun = sum)
  r.tmp[is.na(r.tmp)] <- 0
  if(s == stations[1]) df.r <- r.tmp
  else df.r <- df.r + r.tmp

  
}

## transform to sf object and crs to lat-lon
df.r2 <- df.r
df.r <- st_make_valid(st_as_sf(st_as_stars(df.r), point = FALSE, merge = TRUE, connect8 = TRUE))
colnames(df.r) <- c("dens", "geometry")
df.r <- st_transform(df.r, crs = crs(df))

## intersect df.r with df
# e.g. try this https://gis.stackexchange.com/questions/271268/assigning-raster-values-to-spatial-point-using-r

## get df based on mean estimated positions
df.m <-  df %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()
df.m <- st_as_sf(df.m[!is.na(df.m$lon.m) & !is.na(df.m$lat.m),], coords = c("lon.m", "lat.m"), crs = crs(df))

## get df based on true positions
df.t <-  df %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()
df.t <- st_as_sf(df.t[!is.na(df.t$lon.true) & !is.na(df.t$lat.true),], coords = c("lon.true", "lat.true"), crs = crs(df))

# get raster values for each point
df$dens   <- extract(df.r2, df)
df.m$dens <- extract(df.r2, df.m)
df.t$dens <- extract(df.r2, df.t)

dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raster.gpkg")
st_write(df.r, layer = 'shp_raster', append = F, dsn = dsn)
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_density.gpkg")
st_write(df, layer = 'density', append = F, dsn = dsn)
st_write(df.m, layer = 'density_mean', append = F, dsn = dsn)
st_write(df.t, layer = 'density_gps', append = F, dsn = dsn)


## plot distance error depending on dens (indicator for how good the station cover is)
# ggplot(df.test) + 
#   # geom_point(aes(x = dens, y = distance, color = detR, group = detR),
#   #            pch = 1, alpha = 0.2,
#   #            position = position_dodge(width = 0.5)) +
#   geom_smooth(aes(x = dens, y = distance, group = detR, color = detR), method = "glm", formula = y ~ poly(x, 2)) +
#   geom_smooth(aes(x = dens, y = distance.m, group = detR, color = detR), method = "glm", formula = y ~ poly(x, 2), data = df.test.m, lty = "dashed") +
#   facet_wrap(~Individual, scales = "free_y") +
#   scale_color_viridis_d() +
#   theme_light()

# intersect df.r (raster as polygon) with df, df.m, df.t
## -> this does need a lot of time
# df.int <- sf::st_intersection(df.r, df[df$Individual == "TT090C" & df$meth == "ab_ql" & df$detR == "800m" & df$ant == "1-10",])

#### 3) raster per area with position error ####

if(READ) {
  
  ## create raster (raser_incr*bigger than first raster)
  r <- raster(extent(st_buffer(df.stat, dist = 800)), resolution = raster_resolution*raster_incr)
  
  ## change crs of df to match raster (Gauss-Krueger in m)
  df <- st_transform(df, crs = crs)
  df.m <- st_transform(df.m, crs = crs)
  df.t <- st_transform(df.t, crs = crs)
  
  ## rename distance
  df$dist <- df$distance # use distance or raw positions to true
  df.m$dist <- df.m$distance.m # use distance of mean positions to true
  df.t$dist <- df.t$distance.m # use distance of mean positions to true
  
  ## loop through Individuals, methods, detR, ant
  shp <- NULL
  
  for(s in c("pos", "pos.mean", "pos.gps")) {
    
    ## define source of data
    if(s == "pos")      dat <- df[!is.na(df$dist),]
    if(s == "pos.mean") dat <- df.m[!is.na(df.m$dist),]
    if(s == "pos.gps")  dat <- df.t[!is.na(df.t$dist),]
    
    for(p in unique(dat$Project)) {
      
      for(m in unique(dat$meth[dat$Project == p])) {
        
        for(d in unique(dat$detR[dat$Project == p & dat$meth == m])) {
          
          for(a in unique(dat$ant[dat$Project == p & dat$meth == m & dat$detR == d])) {
            
            ## raster all individuals
            tmp <- dat[dat$Project == p & dat$meth == m & dat$detR == d & dat$ant == a,]
            
            ## get attributes per raster cell (mean, median, sd, ..)
            r.mean   <- raster::rasterize(tmp, r, field = tmp$dist, fun = mean)
            r.median <- raster::rasterize(tmp, r, field = tmp$dist, fun = median)
            r.sd     <- raster::rasterize(tmp, r, field = tmp$dist, fun = sd)
            r.var    <- raster::rasterize(tmp, r, field = tmp$dist, fun = var)
            r.min    <- raster::rasterize(tmp, r, field = tmp$dist, fun = min)
            r.max    <- raster::rasterize(tmp, r, field = tmp$dist, fun = max)
            r.N      <- raster::rasterize(tmp, r, field = tmp$dist, fun = "count")
            ## stack raster files and add variable means
            r.stack <- stack(r.mean, r.median, r.sd, r.var, r.min, r.max, r.N)
            names(r.stack) <- c("mean", "median", "sd", "var", "min", "max", "N")
            ## transform to sf point object and add crs (in m)
            r.stack <- st_as_sf(rasterToPolygons(r.stack), point = FALSE, merge = TRUE, connect8 = TRUE, crs = crs)
            st_crs(r.stack) <- crs
            
            ## add looping variables
            r.stack$source     <- s
            r.stack$Project    <- p
            r.stack$Individual <- "all"
            r.stack$meth       <- m
            r.stack$detR       <- d
            r.stack$ant        <- a
            
            ## merge with previous data
            shp <- rbind(shp, r.stack)
            
            ## raster per individual
            for(i in unique(dat$Individual[dat$Project == p & dat$meth == m & dat$detR == d & dat$ant == a])) {
                
              tmp <- dat[dat$Project == p & dat$Individual == i & dat$meth == m & dat$detR == d & dat$ant == a,]
              
              ## get attributes per raster cell (mean, median, sd, ..)
              r.mean   <- raster::rasterize(tmp, r, field = tmp$dist, fun = mean)
              r.median <- raster::rasterize(tmp, r, field = tmp$dist, fun = median)
              r.sd     <- raster::rasterize(tmp, r, field = tmp$dist, fun = sd)
              r.var    <- raster::rasterize(tmp, r, field = tmp$dist, fun = var)
              r.min    <- raster::rasterize(tmp, r, field = tmp$dist, fun = min)
              r.max    <- raster::rasterize(tmp, r, field = tmp$dist, fun = max)
              r.N      <- raster::rasterize(tmp, r, field = tmp$dist, fun = "count")
              ## stack raster files and add variable means
              r.stack <- stack(r.mean, r.median, r.sd, r.var, r.min, r.max, r.N)
              names(r.stack) <- c("mean", "median", "sd", "var", "min", "max", "N")
              ## transform to sf point object and add crs (in m)
              r.stack <- st_as_sf(rasterToPolygons(r.stack), point = FALSE, merge = TRUE, connect8 = TRUE, crs = crs)
              st_crs(r.stack) <- crs
              
              ## add looping variables
              r.stack$source     <- s
              r.stack$Project    <- p
              r.stack$Individual <- i
              r.stack$meth       <- m
              r.stack$detR       <- d
              r.stack$ant        <- a
              
              ## merge with previous data
              shp <- rbind(shp, r.stack)
            
            } # end of i
            
          } # end of a
          
        } # end of d
        
      } # end of m
      
    } # end of p
    
  } # end of s
  
  ## save data (split by source of data)
  dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raster.gpkg")
  st_write(shp[shp$source == "pos",], layer = 'raster', append = F, dsn = dsn)
  st_write(shp[shp$source == "pos.mean",], layer = 'raster_mean', append = F, dsn = dsn)
  st_write(shp[shp$source == "pos.gps",], layer = 'raster_gps', append = F, dsn = dsn)
  
}

## read data
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raster.gpkg")
shp <- st_read(layer = 'raster', dsn = dsn)
shp.m <- st_read(layer = 'raster_mean', dsn = dsn)
shp.t <- st_read(layer = 'raster_gps', dsn = dsn)

## run markdown to visualize output
rmarkdown::render(input = here("R_pre", "plotError.Rmd"))
# 