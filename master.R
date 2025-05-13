#### master script ####

rm(list = ls())
#### 1) packages ---------------------------------------------------------------
pckgs <- c(
  "ggdist",
  "tidyverse",
  "sf",
  "data.table",
  "here",
  "geosphere",
  "raster", # density as raster
  "stars",
  "terra",
  "ggspatial", # for annotation bar
  "ggnewscale", # adding several color and fill scales to ggplot
  "ggeffects", # to predict gams (in a fast way)
  "parallel", 
  "glmmTMB", # models
  "zoo" # roll mean
  )

## install packages if not already installed
for (i in pckgs) {
  if (!i %in% installed.packages()) {
    install.packages(i, dependencies = TRUE)
  }
}
## load packages
sapply(pckgs, require, character.only = TRUE); rm(pckgs)

####  2) variables -------------------------------------------------------------
## 2a) data --------------------------------------------------------------------
## should data be read in new (TRUE) or loaded (FALSE)
READ <- TRUE # read in or load data
# chose time interval in Jupyter data (usually 2 sec)
t <- "2 sec" # "x sec"

## 2b) spatial data ------------------------------------------------------------
crs     <- 31467 # Gauss Krüger in m
crsPlot <- 3857 # 
crsLL   <- 4326 # lon lat in degree

## 2c) raster data -------------------------------------------------------------
## stations to be included in station cover
stat.direct <- c("c1l1", "c2l1", "c3l1", "c4l1", "c5l1", "c6l1", "c7l1", "c8l1", "c9l1", "c10l1",
              "d1l1", "d2l1", "d3l1", "d4l1", "d5l1", "d6l1", "d7l1", "d8l1")
stat.omni <- c("s1l1", "s2l1", "s3l1", "s4l1", "s5l1", "s6l1", "s7l1", "s8l1", "s9l1", "s10l1")

rast_buf <- 700 # buffer around stations to make raster
rast_res <- 10  # Change resolution as needed (in meter) for stat density raster
rast_inc <- 5  # Change resolution as needed (in meter) for point summary raster (raster_resoultion*raster_incr)

## 2d) plotting ----------------------------------------------------------------
bs <- 14 # basesize

#### 3) data -------------------------------------------------------------------

## stations
shp.stat <- st_read(here("data", "Station.gpkg"), layer = "Station")
shp.stat <- st_transform(shp.stat, crs = crs)
shp.stat$plot <- "station" # to add the shape to legend

## groundtruth data
shp.GPS  <- st_read(here("data", "Testtracks_points.gpkg")) 

#### 4) functions --------------------------------------------------------------
source("plot_functions.R")
source("help_functions.R")

#### 5) run analysis
source("./R_pre/mergeCali.R")
source("./R_model/Model_PE.R") ## check whether all helper functions work (might need some more packages)
