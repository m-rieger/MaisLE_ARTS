#### help functions ####

## Bayes P (P>0)
BayesP <- function(x) sum(x > 0)/length(x)

## rolling mean
rmean <- function(x, width) {
  zoo::rollapply(x,             # column to apply roll mean
                 width = width, # window size of roll mean (e.g. 10, 20, 30)
                 FUN = mean, fill = NA, na.rm = T, partial = TRUE)
}

## rolling max
rmax <- function(x, width) {
  zoo::rollapply(x,             # column to apply roll max
                 width = width, # window size of roll max (e.g. 10, 20, 30)
                 FUN = max, fill = NA, na.rm = T, partial = TRUE)
}

## function to expand to all timestamps to apply rmean and get new position error
rmeanPE <- function(data = NULL, # dataframe with signals and est. raw positions (lon, lat)
                    c.time = "X_time", # timestamp (datetime, POSIXct)
                    c.group = c("site", "tagID", "detR", "meth"), # columns used for grouping
                    c.lon = "lon", c.lat = "lat", # estimated raw positions to average
                    c.add = NULL, # additional columns you want to apply roll mean
                    GT = TRUE, # do you have ground truth data to estimate position errors?
                    c.lonT = "lon.true", c.latT = "lat.true", # groundtruth (GPS) positions, only if GT = T
                    c.PE = "PE", # position error of est. raw positions, only if GT = T
                    w.size = 30) { # window size for roll mean
  
  if(is.null(data)) stop("You need to add a dataframe.")
  
  # get dataframe
  if(GT) {
    d <- data.frame(X_time = as.POSIXct(data[, c.time]),
                    lon = data[, c.lon],
                    lat = data[, c.lat],
                    lonT = data[, c.lonT],
                    latT = data[, c.latT])    
  }
  if(!GT) {
    d <- data.frame(X_time = as.POSIXct(data[, c.time]),
                    lon = data[, c.lon],
                    lat = data[, c.lat])
  }
  
  if(!is.null(c.add)) d[, c.add] <- data[, c.add]

  d$date <- as.Date(d$X_time)
  
  d[, c.group] <- data[, c.group]
  
  df.time <- data.frame()
  
  for(i in unique(d$date)) {
    tmp <- merge(unique(d[d$date == i, c.group]),
                 data.frame(X_time = seq(from = min(d$X_time[d$date == i]), 
                                         to = max(d$X_time[d$date == i]), by = "1 sec")))
    
    df.time <- rbind(df.time, tmp)
    
  }
  
  ## merge with dfs
  d <- left_join(df.time, d, 
                     by = c(c.group, "X_time"))
  
  ## rolling mean positions
  ## here you need the full dataset (including all timestamps with NA lat and lon)
  d <- d %>% group_by(across(all_of(c.group))) %>%
    mutate(
      lon.m = rmean(lon, width = w.size),
      lat.m = rmean(lat, width = w.size),
      Weight.m = if("Weight" %in% names(d)) rmean(Weight, width = w.size) else NA_real_,
      maxSig.m = if("maxSig" %in% names(d)) rmean(maxSig, width = w.size) else NA_real_,
      maxSig.max = if("maxSig" %in% names(d)) rmax(maxSig, width = w.size) else NA_real_,
      Ac.m = if("Ac" %in% names(d)) rmean(Ac, width = w.size) else NA_real_,
      Sc.m = if("Sc" %in% names(d)) rmean(Sc, width = w.size) else NA_real_,
    ) %>%
    ungroup()   
  
  ## Position Error (PE) based on mean positions
  ## here you need only data != na in lon lat (to compute distances)
  d <- d[!is.na(d$lon),]
  
  if(GT) {
    d <- d %>% rowwise %>%
      mutate(PE = distm(x = c(lon.m, lat.m),
                        y = c(lonT, latT))) %>%
      ungroup()    
  }
 
  ## remove or rename duplicate columns to generate identical df to data (but with mean lon lat and new PE)
  d[, c.lon]<- d$lon.m
  d[, c.lat]<- d$lat.m
  colnames(d)[colnames(d) == "X_time"] <- c.time
  colnames(d)[colnames(d) == "PE"] <- c.PE
  
  if(GT) {
    data <- dplyr::select(data, -all_of(c(c.lon, c.lat, c.PE)))
    d <- dplyr::select(d, -c("lon.m", "lat.m", "lonT", "latT", "date"))
  }
  if(!GT) {
    data <- dplyr::select(data, -all_of(c(c.lon, c.lat)))
    d <- dplyr::select(d, -c("lon.m", "lat.m", "date"))
  }
  
  data <- left_join(data, d, by = c(c.group, c.time))
  
  return(data)
}
