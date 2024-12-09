#### merge and summarize testtag positions for calibration ---------------------


#### 0) create other files -----------------------------------------------------
## raster resolution (rast_buf buffer around stations)
rast_ext <- extent(st_buffer(shp.stat, dist = rast_buf))


## mean loctions from shp.GPS (per time.int t based on Jupyter)
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

#### 1) loop through gpkg-files ------------------------------------------------
if(READ) {
  
  ## get all files and create folder for saved files
  df.L <- list.files(paste0("../data/data_", dtyp, "/"), pattern = ".gpkg")
  # df.L <- df.L[85:136]
  if(!dir.exists(paste0("../data/data_", dtyp, "/savedFiles"))) {
    dir.create(paste0("../data/data_", dtyp, "/savedFiles"))}
  
  ## exclude -mean files
  # df.Lm <- df.L[grep("-mean", df.L)]
  df.L <- df.L[-grep("-mean", df.L)]
  
  df.L1 <- df.L[grep("antennabeams.gpkg", df.L)]
  df.L2 <- df.L[grep("antennabeam_s", df.L)]
  df.L3 <- df.L[grep("multilateration", df.L)]
  df.L4 <- df.L[grep("intersections", df.L)]
  
  # df.Lm1 <- df.Lm[grep("antennabeams-mean.gpkg", df.Lm)]
  # df.Lm2 <- df.Lm[grep("antennabeam_s-mean", df.Lm)]
  # df.Lm3 <- df.Lm[grep("multilateration-mean", df.Lm)]
  # df.Lm4 <- df.Lm[grep("intersections-mean", df.Lm)]
  
  #### 1.1) normal data --------------------------------------------------------
  ## antennabeam quadrologger
  df1 <- data.frame()
  
  for(i in df.L1) {
    tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
    
    ## transform crs to lon lat and save crs
    tmp <- st_transform(tmp, crs = crsLL)
    
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
  
  ## antennabeam monologger
  df2 <- data.frame()
  
  for(i in df.L2) {
    tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
    
    ## transform crs to lon lat and save crs
    tmp <- st_transform(tmp, crs = crsLL)
    
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
  
  ## multilateration monologger
  df3 <- data.frame()
  
  for(i in df.L3) {
    tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
    
    ## transform crs to lon lat and save crs
    tmp <- st_transform(tmp, crs = crsLL)
    
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
  
  ## intersection quadrologger
  df4 <- data.frame()
  
  for(i in df.L4) {
    tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
    
    ## transform crs to lon lat and save crs
    tmp <- st_transform(tmp, crs = crsLL)
    
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
  
  df <- do.call("rbind", list(df1, df2, df3, df4))
  rm(df1); rm(df2); rm(df3); rm(df4)
  
  # #### 1.2) mean data ----------------------------------------------------------
  # ## antennabeam quadrologger
  # df1 <- data.frame()
  # 
  # for(i in df.Lm1) {
  #   tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  #   
  #   ## transform crs to lon lat and save crs
  #   tmp <- st_transform(tmp, crs = crsLL)
  #   
  #   ## get Individual and Transmitter
  #   ind <- unique(tmp$Individual)
  #   trans <- unique(tmp$Transmitter)
  #   
  #   ## get coordinates from geometry
  #   tmp <- tmp %>%
  #     dplyr::mutate(lon.m = sf::st_coordinates(.)[,1],
  #                   lat.m = sf::st_coordinates(.)[,2]) %>%
  #     st_drop_geometry()
  #   
  #   ## add cols
  #   tmp$detR <- strsplit(i, "_")[[1]][2]
  #   tmp$ant <- strsplit(i, "_")[[1]][3]
  #   tmp$Distance..m..min <- NA
  #   tmp$meth <- "ab_ql"
  #   
  #   df1 <- rbind(df1, tmp)
  #   
  # }
  # 
  # ## antennabeam monologger
  # df2 <- data.frame()
  # 
  # for(i in df.Lm2) {
  #   tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  #   
  #   ## transform crs to lon lat and save crs
  #   tmp <- st_transform(tmp, crs = crsLL)
  #   
  #   ## get Individual and Transmitter
  #   ind <- unique(tmp$Individual)
  #   trans <- unique(tmp$Transmitter)
  #   
  #   ## get coordinates from geometry
  #   tmp <- tmp %>%
  #     dplyr::mutate(lon.m = sf::st_coordinates(.)[,1],
  #                   lat.m = sf::st_coordinates(.)[,2]) %>%
  #     st_drop_geometry()
  #   
  #   ## add cols
  #   tmp$detR <- "no"
  #   tmp$Distance..m..min <- NA
  #   tmp$ant <- "1-10"
  #   tmp$meth <- "ab_ml"
  #   
  #   df2 <- rbind(df2, tmp)
  #   
  # }
  # 
  # ## multilateration monologger
  # df3 <- data.frame()
  # 
  # for(i in df.Lm3) {
  #   tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  #   
  #   ## transform crs to lon lat and save crs
  #   tmp <- st_transform(tmp, crs = crsLL)
  #   
  #   ## get Individual and Transmitter
  #   ind <- unique(tmp$Individual)
  #   trans <- unique(tmp$Transmitter)
  #   
  #   ## get coordinates from geometry
  #   tmp <- tmp %>%
  #     dplyr::mutate(lon.m = sf::st_coordinates(.)[,1],
  #                   lat.m = sf::st_coordinates(.)[,2]) %>%
  #     st_drop_geometry()
  #   
  #   ## add cols
  #   tmp$detR <- "no"
  #   tmp$Weight <- NA
  #   tmp$ant <- "1-10"
  #   tmp$Antenna.Count <- tmp$Station.Count
  #   tmp$meth <- "ml_ml"
  #   
  #   df3 <- rbind(df3, tmp)
  #   
  # }
  # 
  # ## intersection quadrologger
  # df4 <- data.frame()
  # 
  # for(i in df.Lm4) {
  #   tmp <- st_read(paste0("../data/data_", dtyp, "/", i))
  #   
  #   ## transform crs to lon lat and save crs
  #   tmp <- st_transform(tmp, crs = crsLL)
  # 
  #   ## get Individual and Transmitter
  #   ind <- unique(tmp$Individual)
  #   trans <- unique(tmp$Transmitter)
  #   
  #   ## get coordinates from geometry
  #   tmp <- tmp %>%
  #     dplyr::mutate(lon.m = sf::st_coordinates(.)[,1],
  #                   lat.m = sf::st_coordinates(.)[,2]) %>%
  #     st_drop_geometry()
  #   
  #   ## add cols
  #   tmp$detR <- "no"
  #   tmp$Weight <- NA
  #   tmp$Antenna.Count <- 2*tmp$Station.Count
  #   tmp$ant <- strsplit(i, "_")[[1]][3]
  #   tmp$meth <- "in_ql"
  #   
  #   df4 <- rbind(df4, tmp)
  #   
  # }
  # dfm <- do.call("rbind", list(df1, df2, df3, df4))
  # rm(df1); rm(df2); rm(df3); rm(df4)
  
  df <- dplyr::select(df, -c("Planner", "Transmitter", "Transmitter.Id"))
  
  
  #### 1.2) add Ac and Sc ------------------------------------------------------
  
  # in_ql: Ac >= 4, Sc >= 2
  # ab_ql: Ac >= 1, Sc >= 1
  # ml_ml: Ac >= 1, Sc >= 1, Ac = Sc
  # ab_ml: Ac >= 1, Sc >= 1, Ac = Sc
  df$AcSc <- "ac1-sc1"
  df$AcSc[df$meth == "in_ql"] <- "ac4-sc2"
  
  df.list <- list("df11" = df)
  
  df.list[["df21"]] <- df[df$Antenna.Count >= 2 & df$Station.Count >= 1 & df$meth == "ab_ql",]
  df.list[["df21"]]$AcSc <- "ac2-sc1"
  df.list[["df31"]] <- df[df$Antenna.Count >= 3 & df$Station.Count >= 1 & df$meth == "ab_ql",]
  df.list[["df31"]]$AcSc <- "ac3-sc1"
  df.list[["df41"]] <- df[df$Antenna.Count >= 4 & df$Station.Count >= 1 & df$meth == "ab_ql", ]
  df.list[["df41"]]$AcSc <- "ac4-sc1"
  
  df.list[["df22"]] <- df[df$Antenna.Count >= 2 & df$Station.Count >= 2 & df$meth != "in_ql",]
  df.list[["df22"]]$AcSc <- "ac2-sc2"
  df.list[["df32"]] <- df[df$Antenna.Count >= 3 & df$Station.Count >= 2 & df$meth == "ab_ql",]
  df.list[["df32"]]$AcSc <- "ac3-sc2"
  df.list[["df42"]] <- df[df$Antenna.Count >= 4 & df$Station.Count >= 2 & df$meth == "ab_ql",]
  df.list[["df42"]]$AcSc <- "ac4-sc2"
  
  df.list[["df33"]] <- df[df$Antenna.Count >= 3 & df$Station.Count >= 3 & df$meth != "in_ql",]
  df.list[["df33"]]$AcSc <- "ac3-sc3"
  df.list[["df43"]] <- df[df$Antenna.Count >= 4 & df$Station.Count >= 3 & df$meth == "ab_ql",]
  df.list[["df43"]]$AcSc <- "ac4-sc3"
  
  df.list[["df44"]] <- df[df$Antenna.Count >= 4 & df$Station.Count >= 4 & df$meth != "in_ql",]
  df.list[["df44"]]$AcSc <- "ac4-sc4"
  
  #### 1.3) merge & add GT data, expand time -----------------------------------
  ## add date (to filter each testtrack)
  shp.GPS$date <- as.Date(shp.GPS$X_time)
  df$date <- as.Date(df$X_time)
  colnames(df)[colnames(df) == "Project"] <- "site"
  
  df.time <- data.frame()
  
  ## expand time per testtrack, site, Individual, detR, ant, meth
  for(d in unique(shp.GPS$date)) {
    tmp <- merge(unique(df[df$date == d, c("site", "Individual", "detR", "ant", "meth")]),
                       data.frame(X_time = seq(from = min(shp.GPS$X_time[shp.GPS$date == d]), 
                                    to = max(shp.GPS$X_time[shp.GPS$date == d]), by = t)))
    
    df.time <- rbind(df.time, tmp)
    
  }
  
  
  for(d in names(df.list)){
    
    ## merge with dfs
    df.list[[d]] <- left_join(df.time, df.list[[d]], 
                              by = c("site" = "Project", "X_time", "Individual", "detR", "ant", "meth"))
    
    ## merge GPS with df by site and time
    df.list[[d]] <- left_join(df.list[[d]], 
                    shp.GPS[, c("X_time", "lat.true", "lon.true", "site")],
                    by = c("site", "X_time")) # many-to-many relationship?
    df.list[[d]] <- unique(df.list[[d]])
    
    #### 1.4) rolling mean positions -------------------------------------------
    ## here you need the full dataset (including all timestamps with NA lat and lon)
    df.list[[d]] <- df.list[[d]] %>% group_by(site, Individual, detR, ant, meth, AcSc) %>%
      mutate(
        # lonT.m = rmean(lon.true, width = rollM/2), # not used in Jupyter
        # latT.m = rmean(lat.true, width = rollM/2), # not used in Jupyter   
        lon.m = rmean(lon, width = rollM/2),
        lat.m = rmean(lat, width = rollM/2)
        ) %>%
      ungroup()   
    
    #### 1.5) position error ---------------------------------------------------
    ## here you need only data != na in lon lat (to compute distances)
    df.list[[d]] <- df.list[[d]][!is.na(df.list[[d]]$lon), ] %>% rowwise %>%
      mutate(PE = distm(x = c(lon, lat),
                              y = c(lon.true, lat.true)),
             PE.m = distm(x = c(lon.m, lat.m),
                              y = c(lon.true, lat.true))) %>%
      ungroup()  
  }
  
  # merge dfs
  df <- do.call("rbind", df.list)
  
  ## write df
  df <- st_as_sf(df, coords = c("lon", "lat"), crs = crsLL)
  
  dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raw.gpkg")
  st_write(df, layer = 'unfiltered', append = F, dsn = dsn, 
           layer_options = "GEOMETRY_NAME=geometry")

  ## run markdown to visualize output
  rmarkdown::render(input = here("R_pre", "plotCali.Rmd"), 
                    output_format = "bookdown::html_document2")
}


#### 2) stations density -------------------------------------------------------

## read df
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raw.gpkg")
df <- st_read(layer = 'unfiltered', dsn = dsn)

b100 <- b200 <- b300 <- b400 <- b500 <- b600 <- b700 <- shp.buffer <- list()
df.r <- NULL

# Create an empty raster
r <- raster(rast_ext, resolution = rast_res)

for(i in 1:length(stations)) {
  s <- stations[i]
  
  df.p <- shp.stat[shp.stat$station.id == s,]
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
df.m <- st_as_sf(df.m[!is.na(df.m$lon.m) & !is.na(df.m$lat.m),], 
                 coords = c("lon.m", "lat.m"), crs = crs(df))

## get df based on true positions
df.t <-  df %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry()
df.t <- st_as_sf(df.t[!is.na(df.t$lon.true) & !is.na(df.t$lat.true),], 
                 coords = c("lon.true", "lat.true"), crs = crs(df))

# get raster values for each point
df$dens   <- extract(df.r2, df)
df.m$dens <- extract(df.r2, df.m)
df.t$dens <- extract(df.r2, df.t)

dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raster.gpkg")
st_write(df.r, layer = 'shp_raster', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_density.gpkg")
st_write(df, layer = 'density', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")
st_write(df.m, layer = 'density_mean', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")
st_write(df.t, layer = 'density_gps', append = F, dsn = dsn, 
         layer_options = "GEOMETRY_NAME=geometry")


## plot PE error depending on dens (indicator for how good the station cover is)
# ggplot(df.test) + 
#   # geom_point(aes(x = dens, y = PE, color = detR, group = detR),
#   #            pch = 1, alpha = 0.2,
#   #            position = position_dodge(width = 0.5)) +
#   geom_smooth(aes(x = dens, y = PE, group = detR, color = detR), method = "glm", formula = y ~ poly(x, 2)) +
#   geom_smooth(aes(x = dens, y = PE.m, group = detR, color = detR), method = "glm", formula = y ~ poly(x, 2), data = df.test.m, lty = "dashed") +
#   facet_wrap(~Individual, scales = "free_y") +
#   scale_color_viridis_d() +
#   theme_light()

# intersect df.r (raster as polygon) with df, df.m, df.t
## -> this does need a lot of time
# df.int <- sf::st_intersection(df.r, df[df$Individual == "TT090C" & df$meth == "ab_ql" & df$detR == "800m" & df$ant == "1-10",])

#### 3) raster per area with position error ------------------------------------

if(READ) {
  
  ## create raster (rast_inc*bigger than first raster)
  r <- raster(extent(st_buffer(shp.stat, dist = rast_buf)), resolution = rast_res*rast_inc)
  
  ## change crs of df to match raster (Gauss-Krueger in m)
  df <- st_transform(df, crs = crs)
  df.m <- st_transform(df.m, crs = crs)
  df.t <- st_transform(df.t, crs = crs)
  
  ## rename PE
  df$dist <- df$PE # use PE or raw positions to true
  df.m$dist <- df.m$PE.m # use PE of mean positions to true
  df.t$dist <- df.t$PE.m # use PE of mean positions to true
  
  ## loop through Individuals, methods, detR, AcSc
  shp <- NULL
  
  for(s in c("pos", "pos.mean", "pos.gps")) {
    
    ## define source of data
    if(s == "pos")      dat <- df[!is.na(df$dist),]
    if(s == "pos.mean") dat <- df.m[!is.na(df.m$dist),]
    if(s == "pos.gps")  dat <- df.t[!is.na(df.t$dist),]
    
    for(p in unique(dat$site)) {
      
      for(m in unique(dat$meth[dat$site == p])) {
        
        for(d in unique(dat$detR[dat$site == p & dat$meth == m])) {
          
          for(a in unique(dat$AcSc[dat$site == p & dat$meth == m & dat$detR == d])) {
            
            ## raster all individuals
            tmp <- dat[dat$site == p & dat$meth == m & dat$detR == d & dat$AcSc == a,]
            
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
            r.stack$site    <- p
            r.stack$Individual <- "all"
            r.stack$meth       <- m
            r.stack$detR       <- d
            r.stack$AcSc        <- a
            
            ## merge with previous data
            shp <- rbind(shp, r.stack)
            
            ## raster per individual
            for(i in unique(dat$Individual[dat$site == p & dat$meth == m & dat$detR == d & dat$AcSc == a])) {
                
              tmp <- dat[dat$site == p & dat$Individual == i & dat$meth == m & dat$detR == d & dat$AcSc == a,]
              
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
              r.stack$site    <- p
              r.stack$Individual <- i
              r.stack$meth       <- m
              r.stack$detR       <- d
              r.stack$AcSc        <- a
              
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
  st_write(shp[shp$source == "pos",], layer = 'raster', append = F, dsn = dsn, 
           layer_options = "GEOMETRY_NAME=geometry")
  st_write(shp[shp$source == "pos.mean",], layer = 'raster_mean', append = F, dsn = dsn, 
           layer_options = "GEOMETRY_NAME=geometry")
  st_write(shp[shp$source == "pos.gps",], layer = 'raster_gps', append = F, dsn = dsn, 
           layer_options = "GEOMETRY_NAME=geometry")
  
}

## read data
dsn <- paste0("../data/data_", dtyp, "/savedFiles/Data_cali_raster.gpkg")
shp <- st_read(layer = 'raster', dsn = dsn)
shp.m <- st_read(layer = 'raster_mean', dsn = dsn)
shp.t <- st_read(layer = 'raster_gps', dsn = dsn)

## run markdown to visualize output
rmarkdown::render(input = here("R_pre", "plotError.Rmd"))
# 