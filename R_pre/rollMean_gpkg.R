#### apply rolling mean and NStat, NAnt filters to gpkg-files from Jupyter ####

## input
########-
## needed folder structure:
# data/: folder containing two subfolders:
 # data/data_animal/: this is where you store the gpkg-files for individuals/birds
 # data/data_testtag/: this is where you store the gpkg-files for testtags
# MaisLE_ARTS/: this is the complete repository from github

# data and MaisLE_ARTS should be neighboring folders

## output
#########-
# output is saved in
  # data/data_animal/savedFiles/ or
  # data/data_testtag/savedFiles/
# there are different types of output:
  # *.csv:  contains the unfiltered data for no rolling mean, rolling mean with 15s and 30s
  # *_filtered.csv:  as above but contains the filtered data (NStat, NAnt) 
  # *.gpkg: only if gpkg = TRUE, gpkg with 6 layers
    # unfiltered
    # unfiltered_rollM15
    # unfiltered_rollM30
    # filtered
    # filtered_rollM15
    # filtered_rollM30
  # *.jpeg: display of calculated positions per filter and rolling mean
  # *_Error.jpeg: (only for testtags) histogram of distance error per filter and rolling mean

## define variables (here, you might want to change something)
# chose data type (animal data or testtag data)
dtyp <- "testtag" # "animal" or "testtag"
# chose time interval in Jupyter data (usually 2 sec)
t <- "2 sec" # "x sec"
# minimum number of stations per position
NStat <- 2
# minimum number of antennas per position
NAnt <- 3
# write gpkg for rolling mean and filtered data
gpkg <- TRUE

## no further changes needed down here :) ##

## packages
library(ggdist)
library(tidyverse)
library(sf)
#library(data.table)
library(geosphere)
library(zoo)

## additional variables
t2 <- as.numeric(gsub(" sec", "", t))

#### loop through gpkg-files ####

## get all files
df.L <- list.files(paste0("../data/data_", dtyp, "/"), pattern = ".gpkg")

# create folder to store new files
if(!dir.exists(paste0("../data/data_", dtyp, "/savedFiles"))) dir.create(paste0("../data/data_", dtyp, "/savedFiles"))

for(i in df.L) {
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
  
  ## clone tmp for filtered data
  tmp.f <- tmp[tmp$Station.Count >= NStat & tmp$Antenna.Count >= NAnt,]
  
  ## get time sequence and order by time
  tmp <- full_join(tmp, 
                   data.frame(X_time = seq(from = min(tmp$X_time)-14, 
                                           to = max(tmp$X_time)+14, by = t)), 
                   by = c("X_time"))
  tmp <- tmp[order(tmp$X_time),]
  
  tmp.f <- full_join(tmp.f, 
                   data.frame(X_time = seq(from = min(tmp.f$X_time)-14, 
                                           to = max(tmp.f$X_time)+14, by = t)), 
                   by = c("X_time"))
  tmp.f <- tmp.f[order(tmp.f$X_time),]
  
  ## calculate rolling mean for 15 and 30 seconds
  tmp$lon15 <- zoo::rollapply(tmp$lon, width = 15/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  tmp$lat15 <- zoo::rollapply(tmp$lat, width = 15/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  tmp$lon30 <- zoo::rollapply(tmp$lon, width = 30/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  tmp$lat30 <- zoo::rollapply(tmp$lat, width = 30/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  
  tmp.f$lon15 <- zoo::rollapply(tmp.f$lon, width = 15/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  tmp.f$lat15 <- zoo::rollapply(tmp.f$lat, width = 15/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  tmp.f$lon30 <- zoo::rollapply(tmp.f$lon, width = 30/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  tmp.f$lat30 <- zoo::rollapply(tmp.f$lat, width = 30/t2, FUN = mean, fill = NA, na.rm = T, partial = TRUE)
  
  ## fill columns
  tmp$Individual <- tmp.f$Individual <- ind
  tmp$Transmitter <- tmp.f$Transmitter <- trans
  
  ## get real coordinates for comparison
  if(dtyp == "testtag") {
    
    tmp <- tmp %>% rowwise() %>%
      mutate(error = distm(x = c(lon, lat),
                              y = c(Longitude.Calibration, Latitude.Calibration)),
             error15 = distm(x = c(lon15, lat15),
                           y = c(Longitude.Calibration, Latitude.Calibration)),
             error30 = distm(x = c(lon30, lat30),
                           y = c(Longitude.Calibration, Latitude.Calibration))) %>%
      ungroup()
    
    tmp.f <- tmp.f %>% rowwise() %>%
      mutate(error = distm(x = c(lon, lat),
                           y = c(Longitude.Calibration, Latitude.Calibration)),
             error15 = distm(x = c(lon15, lat15),
                             y = c(Longitude.Calibration, Latitude.Calibration)),
             error30 = distm(x = c(lon30, lat30),
                             y = c(Longitude.Calibration, Latitude.Calibration))) %>%
      ungroup()
    
  }
  
  ## remove unnecessary columns
  column <- c("X_time", "Station.Count", "Antenna.Count", "Stations", "Weight", "Signal.max",
              "Transmitter", "Individual",
              "lon", "lat", "lon15", "lat15", "lon30", "lat30")
  if(dtyp == "testtag") column <- c(column, "Latitude.Calibration", "Longitude.Calibration", "error", "error15", "error30")

  tmp   <- tmp[, column]
  tmp.f <- tmp.f[, column]
  
  ## write data
  write.csv(tmp, paste0("../data/data_", dtyp, "/savedFiles/", gsub(".gpkg", "", i), ".csv"), row.names = F)
  write.csv(tmp.f, paste0("../data/data_", dtyp, "/savedFiles/", gsub(".gpkg", "", i), "_filtered.csv"), row.names = F)
  
  ## get shp per filter and rollMean
  tmp1 <- st_as_sf(tmp[!is.na(tmp$lon) & !is.na(tmp$lat),], coords = c("lon", "lat"), crs = crs)
  tmp15 <- st_as_sf(tmp[!is.na(tmp$lon15) & !is.na(tmp$lat15),], coords = c("lon15", "lat15"), crs = crs)
  tmp30 <- st_as_sf(tmp[!is.na(tmp$lon30) & !is.na(tmp$lat30),], coords = c("lon30", "lat30"), crs = crs)
  
  tmp.f1 <- st_as_sf(tmp.f[!is.na(tmp.f$lon) & !is.na(tmp.f$lat),], coords = c("lon", "lat"), crs = crs)
  tmp.f15 <- st_as_sf(tmp.f[!is.na(tmp.f$lon15) & !is.na(tmp.f$lat15),], coords = c("lon15", "lat15"), crs = crs)
  tmp.f30 <- st_as_sf(tmp.f[!is.na(tmp.f$lon30) & !is.na(tmp.f$lat30),], coords = c("lon30", "lat30"), crs = crs)
  
  ## get distance error
  if(dtyp == "testtag") {
    tmp1$error_m <- tmp1$error
    tmp15$error_m <- tmp15$error15
    tmp30$error_m <- tmp30$error30
    tmp.f1$error_m <- tmp.f1$error
    tmp.f15$error_m <- tmp.f15$error15
    tmp.f30$error_m <- tmp.f30$error30    
  }
  
  ## remove unecessary columns
  column <- c("X_time", "Station.Count", "Antenna.Count", "Stations", "Weight", "Signal.max",
                                   "Transmitter", "Individual")
  if(dtyp == "testtag") column <- c(column, "Latitude.Calibration", "Longitude.Calibration", "error_m")
  
  tmp1 <- tmp1[, column]
  tmp15 <- tmp15[, column]
  tmp30 <- tmp30[, column]
  tmp.f1 <- tmp.f1[, column]
  tmp.f15 <- tmp.f15[, column]
  tmp.f30 <- tmp.f30[, column]  
  
  ## write gpkg per filter and rollMean
  if(gpkg) {
    dsn <- paste0("../data/data_", dtyp, "/savedFiles/", i)
    st_write(tmp1, layer = 'unfiltered', append = F, dsn = dsn)
    st_write(tmp15, layer = 'unfiltered_rollM15', append = F, dsn = dsn)
    st_write(tmp30, layer = 'unfiltered_rollM30', append = F, dsn = dsn)
    st_write(tmp.f1, layer = 'filtered', append = F, dsn = dsn)
    st_write(tmp.f15, layer = 'filtered_rollM15', append = F, dsn = dsn)
    st_write(tmp.f30, layer = 'filtered_rollM30', append = F, dsn = dsn)
  }
  
  ## add rollMean column
  tmp1$rollM <- paste0(1, "S, ", 1, "A - rollM 1")
  tmp15$rollM <- paste0(1, "S, ", 1, "A - rollM 15")
  tmp30$rollM <- paste0(1, "S, ", 1, "A - rollM 30")
  tmp.f1$rollM <- paste0(NStat, "S, ", NAnt, "A - rollM 1")
  tmp.f15$rollM <- paste0(NStat, "S, ", NAnt, "A - rollM 15")
  tmp.f30$rollM <- paste0(NStat, "S, ", NAnt, "A - rollM 30")
    
  ## merge data for plotting
  column <- c("X_time", "Antenna.Count", "Station.Count", "rollM")
  if(dtyp == "testtag") column <- c(column, "error_m")
  
  df.plot <- rbind(tmp1[, column],
                   tmp15[, column])
  df.plot <- rbind(df.plot,
                   tmp30[, column])
  df.plot <- rbind(df.plot,
                   tmp.f1[, column])  
  df.plot <- rbind(df.plot,
                   tmp.f15[, column])
  df.plot <- rbind(df.plot,
                   tmp.f30[, column])  
  
  ## get number of points per method and merge with rollM
  sum.df <- st_drop_geometry(df.plot) %>% group_by(rollM) %>%
    summarize(nPos = n(), .groups = "drop")
  
  if(dtyp == "testtag") {
    sum.df <- st_drop_geometry(df.plot) %>% group_by(rollM) %>%
      summarize(nPos = n(), .groups = "drop",
                error_median = median(error_m, na.rm = T),
                error_25 = quantile(error_m, probs = 0.25, na.rm = T),
                error_75 = quantile(error_m, probs = 0.75, na.rm = T))
  }
  
  sum.df$group <- paste0(sum.df$rollM, " (n = ", sum.df$nPos, ")")
  
  df.plot <- left_join(df.plot, sum.df[, c("rollM", "group")], by = "rollM")
  
  ## plot
  jpeg(paste0("../data/data_", dtyp, "/savedFiles/", gsub(".gpkg", "", i), ".jpeg"), 
       height = 1200, width = 1500)
  print(ggplot(df.plot) +
    geom_sf(aes(color = X_time), pch = 1, alpha = 0.5, size = 3) +
    facet_wrap(~group, ncol = 3) +
    scale_color_viridis_c(end = 0.9) +
    ggtitle(gsub(".gpkg", "", i)) +
    theme_light(base_size = 25))
  dev.off()
  
  if(dtyp == "testtag") {
    ## plot error
    jpeg(paste0("../data/data_", dtyp, "/savedFiles/", gsub(".gpkg", "", i), "_Error.jpeg"), 
         height = 1000, width = 1000)
    print(df.plot %>%
        ggplot(aes(x = error_m, y = group)) +
        stat_dots() +
        ggtitle(gsub(".gpkg", "", i)) +
        stat_slabinterval(point_size = 5) +
        theme_light(base_size = 15))
    dev.off()
  }
  
}  
