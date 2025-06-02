#### merge and summarize testtag positions for calibration ---------------------


#### 0) create other files -----------------------------------------------------
## raster resolution (rast_buf buffer around stations)
rast_ext <- extent(st_buffer(shp.stat, dist = rast_buf*1.1))

## mean locations from shp.GPS (per time.int t based on Jupyter)
shp.GPS <- shp.GPS %>%
  dplyr::mutate(lon.true = sf::st_coordinates(.)[,1],
                lat.true = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()
# round to t
shp.GPS$X_time <- floor_date(shp.GPS$time, t)
# if uneven, add 1 (since Jupyter data is uneven)
if(second(shp.GPS$X_time[1]) %% 2 == 0) shp.GPS$X_time <- shp.GPS$X_time + 1
# summarize per X_time
shp.GPS <- shp.GPS %>% group_by(X_time, site) %>%
  summarise(lat.true = mean(lat.true, na.rm = T),
            lon.true = mean(lon.true, na.rm = T), .groups = "drop")

#### 1) merge estimated positions with true positions --------------------------
if(READ) {
  ## create folder for saved gpkg
  if(!dir.exists("./data/cali/savedFiles")) dir.create("./data/cali/savedFiles")
  
  #### 1.1) read in data (raw data output from position calculation) -----------
  df1 <- st_read(dsn = "./data/cali/Data_cali_raw_maisC.gpkg", layer = "unfiltered")
  df2 <- st_read(dsn = "./data/cali/Data_cali_raw_maisD.gpkg", layer = "unfiltered")
  
  df <- rbind(df1, df2); rm(df1); rm(df2)
  
  #### 1.2) merge & add GT data, expand time -----------------------------------
  ## add date (to filter each testtrack)
  shp.GPS$date <- as.Date(shp.GPS$X_time)
  df$date <- as.Date(df$X_time)
  colnames(df)[colnames(df) == "Project"] <- "site"
  
  ## merge GPS with df by site and time
  df <- left_join(df, 
                  shp.GPS[, c("X_time", "lon.true", "lat.true", "site")],
                  by = c("site", "X_time"))
  df <- df[!is.na(df$lat.true),] # due to gaps in test tracks
  
  #### 1.3) position error -----------------------------------------------------
  df <-  df %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    st_drop_geometry()
  
  ## here you need only data != na in lon lat (to compute distances)
  # NA or NaN values exist for positions with received signal(s) but a failed 
  # calculation, e.g. because intersection lines do not intersect (direct an)
  df <- df[!is.na(df$lon), ] %>% rowwise %>%
    mutate(PE = distm(x = c(lon, lat),
                            y = c(lon.true, lat.true))) %>%
    ungroup()  

  ## write df as gpkg
  df <- st_as_sf(df, coords = c("lon", "lat"), crs = crsLL)
  
  dsn <- "./data/cali/savedFiles/Data_cali_raw.gpkg"
  st_write(df, layer = 'unfiltered', append = F, dsn = dsn, 
           layer_options = "GEOMETRY_NAME=geometry")

}


#### 2) stations density -------------------------------------------------------

## read df
dsn <- "./data/cali/savedFiles/Data_cali_raw.gpkg"
df <- st_read(layer = 'unfiltered', dsn = dsn)

shp.buffer <- list()
r.direct <- r.omni <- NULL

#### 2.1) create station density raster ----------------------------------------
## Create an empty raster
r <- raster(rast_ext, resolution = rast_res)

radii <- seq(100, rast_buf, by = 100)
densities <- seq(1, 0.1, length.out = length(radii))

## calculate station density for directional stations
for (s in stat.direct) {
  df.p <- shp.stat[shp.stat$station.id == s,]
  buf <- list()
  
  ## create buffer circles
  for (i in seq_along(radii)) buf[[i]] <- st_buffer(df.p, dist = radii[i])
  
  ## create buffer rings (based on circles)
  for (i in length(buf):2) buf[[i]] <- st_difference(buf[[i]], buf[[i - 1]])
  
  ## assign density
  for (i in seq_along(buf)) {
    buf[[i]]$dens <- densities[i]
    buf[[i]] <- buf[[i]][, "dens"]
  }
  
  ## merge buffers
  shp.buffer[[s]] <- do.call(rbind, buf)
  
  ## rasterize
  r.tmp <- raster::rasterize(shp.buffer[[s]], r, field = "dens", fun = sum)
  r.tmp[is.na(r.tmp)] <- 0
  
  ## and sum up rasters
  if (s == stat.direct[1]) r.direct <- r.tmp
  else r.direct <- r.direct + r.tmp
}

## calculate station density for omnidirectional stations
for (s in stat.omni) {
  df.p <- shp.stat[shp.stat$station.id == s,]
  buf <- list()
  
  ## create buffer circles
  for (i in seq_along(radii)) buf[[i]] <- st_buffer(df.p, dist = radii[i])
  
  ## create buffer rings (based on circles)
  for (i in length(buf):2) buf[[i]] <- st_difference(buf[[i]], buf[[i - 1]])
  
  ## assign density
  for (i in seq_along(buf)) {
    buf[[i]]$dens <- densities[i]
    buf[[i]] <- buf[[i]][, "dens"]
  }
  
  ## merge buffers
  shp.buffer[[s]] <- do.call(rbind, buf)
  
  ## rasterize
  r.tmp <- raster::rasterize(shp.buffer[[s]], r, field = "dens", fun = sum)
  r.tmp[is.na(r.tmp)] <- 0
  
  ## and sum up rasters
  if (s == stat.omni[1]) r.omni <- r.tmp
  else r.omni <- r.omni + r.tmp
}

## merge density rings in one shape, label station type
shp.buffer <- do.call(rbind, Map(cbind, shp.buffer, station = names(shp.buffer)))
shp.buffer$type <- "direct"
shp.buffer$type[shp.buffer$station %in% stat.omni] <- "omni"

## save shp file
dsn <- paste0("./data/cali/savedFiles/Data_cali_raster.gpkg")
st_write(shp.buffer, layer = 'density_rings', append = F, dsn = dsn)

## transform raster to sf object and crs to lat-lon
r.direct2 <- st_make_valid(st_as_sf(st_as_stars(r.direct), point = FALSE, merge = TRUE, connect8 = TRUE))
r.direct2$type = "direct"

r.omni2 <- st_make_valid(st_as_sf(st_as_stars(r.omni), point = FALSE, merge = TRUE, connect8 = TRUE))
r.omni2$type = "omni"

r.shp <- rbind(r.direct2, r.omni2); rm(r.direct2); rm(r.omni2)

colnames(r.shp) <- c("dens", "type", "geometry")

r.shp <- st_transform(r.shp, crs = crs(df))
st_write(r.shp, layer = 'density_raster', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")

#### 2.2) intersect with data (estimated and true) -----------------------------
## add type to df
df$type <- "direct"
df$type[df$meth %in% c("omni.ab", "omni.ml")] <- "omni"

## get df based on true positions
df.t <-  df %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()
df.t <- st_as_sf(df.t[!is.na(df.t$lon.true) & !is.na(df.t$lat.true),], 
                 coords = c("lon.true", "lat.true"), crs = crs(df))

# get raster values for each point
df$dens[df$type == "direct"] <- extract(r.direct, df[df$type == "direct",])
df.t$dens[df.t$type == "direct"] <- extract(r.direct, df.t[df.t$type == "direct",])
df$dens[df$type == "omni"] <- extract(r.omni, df[df$type == "omni",])
df.t$dens[df.t$type == "omni"] <- extract(r.omni, df.t[df.t$type == "omni",])

dsn <- paste0("./data/cali/savedFiles/Data_cali_density.gpkg")
st_write(df, layer = 'density_est', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")
st_write(df.t, layer = 'density_true', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")

#### 3) raster per area with position error ------------------------------------

if(READ) {
  
  ## create raster (rast_inc*bigger than first raster)
  r <- raster(extent(st_buffer(shp.stat, dist = rast_buf*1.1)), resolution = rast_res*rast_inc)
  
  ## change crs of df to match raster (Gauss-Krueger in m)
  df <- st_transform(df, crs = crs)
  df.t <- st_transform(df.t, crs = crs)
  
  ## loop through Individuals, methods, r
  shp <- NULL
  
  for(s in c("pos", "pos.gps")) {
    
    ## define source of data
    if(s == "pos")      dat <- df[!is.na(df$PE),]
    if(s == "pos.gps")  dat <- df.t[!is.na(df.t$PE),]
    
    for(p in unique(dat$site)) {
      
      for(m in unique(dat$meth[dat$site == p])) {
        
        for(detr in unique(dat$r[dat$site == p & dat$meth == m])) {
          
          ## raster all individuals
          tmp <- dat[dat$site == p & dat$meth == m & dat$r == detr,]
          
          ## get attributes per raster cell (mean, median, sd, ..)
          r.mean   <- raster::rasterize(tmp, r, field = tmp$PE, fun = mean)
          r.median <- raster::rasterize(tmp, r, field = tmp$PE, fun = median)
          r.sd     <- raster::rasterize(tmp, r, field = tmp$PE, fun = sd)
          r.var    <- raster::rasterize(tmp, r, field = tmp$PE, fun = var)
          r.min    <- raster::rasterize(tmp, r, field = tmp$PE, fun = min)
          r.max    <- raster::rasterize(tmp, r, field = tmp$PE, fun = max)
          r.N      <- raster::rasterize(tmp, r, field = tmp$PE, fun = "count")
          ## stack raster files and add variable means
          r.stack <- stack(r.mean, r.median, r.sd, r.var, r.min, r.max, r.N)
          names(r.stack) <- c("mean", "median", "sd", "var", "min", "max", "N")
          ## transform to sf point object and add crs (in m)
          r.stack <- st_as_sf(rasterToPolygons(r.stack), point = FALSE, merge = TRUE, connect8 = TRUE, crs = crs)
          st_crs(r.stack) <- crs
          
          ## add looping variables
          r.stack$source     <- s
          r.stack$site       <- p
          r.stack$Individual <- "all"
          r.stack$meth       <- m
          r.stack$r          <- detr

          ## merge with previous data
          shp <- rbind(shp, r.stack)
          
          ## raster per individual
          for(i in unique(dat$Individual[dat$site == p & dat$meth == m & dat$r == detr])) {
              
            tmp <- dat[dat$site == p & dat$Individual == i & dat$meth == m & dat$r == detr,]
            
            ## get attributes per raster cell (mean, median, sd, ..)
            r.mean   <- raster::rasterize(tmp, r, field = tmp$PE, fun = mean)
            r.median <- raster::rasterize(tmp, r, field = tmp$PE, fun = median)
            r.sd     <- raster::rasterize(tmp, r, field = tmp$PE, fun = sd)
            r.var    <- raster::rasterize(tmp, r, field = tmp$PE, fun = var)
            r.min    <- raster::rasterize(tmp, r, field = tmp$PE, fun = min)
            r.max    <- raster::rasterize(tmp, r, field = tmp$PE, fun = max)
            r.N      <- raster::rasterize(tmp, r, field = tmp$PE, fun = "count")
            ## stack raster files and add variable means
            r.stack <- stack(r.mean, r.median, r.sd, r.var, r.min, r.max, r.N)
            names(r.stack) <- c("mean", "median", "sd", "var", "min", "max", "N")
            ## transform to sf point object and add crs (in m)
            r.stack <- st_as_sf(rasterToPolygons(r.stack), point = FALSE, merge = TRUE, connect8 = TRUE, crs = crs)
            st_crs(r.stack) <- crs
            
            ## add looping variables
            r.stack$source     <- s
            r.stack$site       <- p
            r.stack$Individual <- i
            r.stack$meth       <- m
            r.stack$r          <- detr

            ## merge with previous data
            shp <- rbind(shp, r.stack)
          
          } # end of i
            
        } # end of d
        
      } # end of m
      
    } # end of p
    
  } # end of s
  
  ## save data (split by source of data)
  dsn <- paste0("./data/cali/savedFiles/Data_cali_raster.gpkg")
  st_write(shp, layer = 'raster_site', append = F, dsn = dsn, 
           layer_options = "GEOMETRY_NAME=geometry")

  
}

## read data
dsn <- paste0("./data/cali/savedFiles/Data_cali_raster.gpkg")
shp <- st_read(layer = 'raster_site', dsn = dsn)

## run markdown to visualize output
rmarkdown::render(input = here("R", "plotError.Rmd"))
