#### master script ####

rm(list = ls())
#install.packages("raster")
#### 1) packages ---------------------------------------------------------------
library(ggdist)
library(tidyverse)
library(sf)
library(data.table)
library(here)
library(geosphere)
library(lubridate)
#library(spatstat) # density distribution for sf objects
library(raster) # density as raster
library(stars)
library(terra)
library(ggspatial) # for annotation bar
library(ggnewscale) # adding several (color, fill) scales to ggplot
library(zoo) # roll mean


####  2) variables -------------------------------------------------------------
## 2a) data --------------------------------------------------------------------
## should data be read in new (TRUE) or loaded (FALSE)
READ <- TRUE # read in or load data
# chose data type (animal data or testtag data)
dtyp <- "cali" # "animal", "testtag" or "cali"
# chose time interval in Jupyter data (usually 2 sec)
t <- "2 sec" # "x sec"
# minimum number of stations per position
NStat <- 1
# minimum number of antennas per position
NAnt <- 1
# rolling mean for mean positions
rollM <- 15

## 2b) spatial data ------------------------------------------------------------
crs <- 31467 # Gauss Krüger in m
crsPlot <- 3857 # 
crsLL <- 4326 # lon lat in degree

## 2c) raster data -------------------------------------------------------------
## stations to be included in stations density
stations <- c("c1l1", "c2l1", "c3l1", "c4l1", "c5l1", "c6l1", "c7l1", "c8l1", "c9l1", "c10l1",
              "d1l1", "d2l1", "d3l1", "d4l1", "d5l1", "d6l1", "d7l1", "d8l1")


rast_buf <- 800 # buffer around stations to make raster
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

rmean <- function(x, width) {
  zoo::rollapply(x, width = width, 
                 FUN = mean, fill = NA, na.rm = T, partial = TRUE)
}

#### 5) run analysis
source("./R_pre/mergeCali.R")
