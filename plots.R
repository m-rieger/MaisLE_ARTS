## plots for paper

rm(list = ls())

## packages
library(tidyverse)
library(glmmTMB)
library(sf)
library(basemaps) # for loading OSM background maps
library(here)
library(ggspatial) # for annotation bar and annotation_map_tiles
library(prettymapr) # needed for annotation_map_tiles
library(gridExtra)
library(grid)
library(ggnewscale) # adding several (color, fill) scales to ggplot
library(patchwork) # for combining plots
library(ggdist) # plot densities
# devtools::install_github("psyteachr/introdataviz")
library(introdataviz)
library(spatstat) # kernel density of sf object

source("./R_model/Linear modelling workflow_support functions.R") 
source("./help_functions.R") 

## variables
crs <- 31467 # Gauss Krüger in m
crsPlot <- 3857
crsLL <- 4326 # lon lat in degree

#### 0) load data --------------------------------------------------------------
## station
df.stat <- read.csv("./data/station.csv", encoding = "latin1")
df.stat <- st_as_sf(x = df.stat[df.stat$analysis == "yes",],                         
                    coords = c("station.lon", "station.lat"),
                    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    remove = F)
df.stat$plot <- "station" # to add the shape to legend
df.stat$type <- "omni"
df.stat$type[df.stat$station.type == "VHF Directional Station"] <- "direct"
df.stat <- st_transform(df.stat, crs = crsPlot) 

# df.stat <- df.stat %>% group_by(station.project_id, type) %>%
#   mutate(minDist = min(st_distance(geometry))) %>% ungroup()

# table for sites
tab.sites <- read.csv("./data/tabSites.csv", encoding = "latin1")

## get dimensions for plotting
## get plotting dimensions
dimC3 <- st_bbox(st_buffer(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisC",])), dist = 1300))
dimC2 <- st_bbox(st_buffer(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisC",])), dist = 2000))
dimC <- st_bbox(st_buffer(df.stat[df.stat$station.id == "c5l1",], dist = 1600))
dimD3 <- st_bbox(st_buffer(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisD",])), dist = 1300))
dimD2 <- st_bbox(st_buffer(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisD",])), dist = 2000))
dimD <- st_bbox(st_buffer(df.stat[df.stat$station.id == "d6l1",], dist = 1600))
dimA <- st_bbox(st_buffer(df.stat, dist = 200000))

## two study area
df.area <- rbind(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisC",])),
                 st_centroid(st_union(df.stat[df.stat$station.project_id == "maisD",])))

## isolines
# shp.iso <- read_sf("../data/Isolines_BB.shp")
# shp.iso <- st_zm(st_transform(shp.iso, crs = crsPlot)) # change crs and XYM to XY 
# isoC <- st_intersection(shp.iso, st_as_sfc(dimC))
# isoD <- st_intersection(shp.iso, st_as_sfc(dimD))
# isoA <- st_intersection(shp.iso, st_as_sfc(dimA))
# isoA <- isoA[isoA$z %in% seq(0, 100, 5),]
# rm(shp.iso)
# st_write(isoC, here("data", "Isolines_C.shp"), append = F) 
# st_write(isoD, here("data", "Isolines_D.shp"), append = F) 
# st_write(isoA, here("data", "Isolines_A.shp"), append = F) 

isoC <- read_sf(here("data", "Isolines_C.shp")) 
isoD <- read_sf(here("data", "Isolines_D.shp")) 
# isoA <- read_sf(here("data", "Isolines_A.shp")) 
shp.ger <- read_sf(here("data", "vg2500_bld_ganz.shp"))
shp.ger <- st_transform(shp.ger, crs = crsPlot) 
dimA <- st_bbox(st_buffer(shp.ger, dist = 200000))

## testtracks (points, lines)
shp.p <- st_read(here("data", "Testtracks_points.gpkg")) 
shp.l <- st_read(here("data", "Testtracks_lines.gpkg")) 

## raster
shp.r <- st_read(dsn = "../data/data_cali/savedFiles/Data_cali_raster.gpkg", layer = "shp_raster")
shp.r <- st_transform(shp.r, crs = crsPlot)

## circles (for raster)
shp.circ <- st_read(dsn = "../data/data_cali/savedFiles/Data_cali_circles.gpkg")

# model results and dfs
df.pred1 <- read.csv(here("output_model", "model-predictions_m1.csv"))
df.sim1  <- read.csv(here("output_model", "model-simulations_m1.csv"))
df.coef1 <- read.csv(here("output_model", "model-coefficients_m1.csv"))

df.pred4 <- read.csv(here("output_model", "model-predictions_m4.csv"))
df.sim4  <- read.csv(here("output_model", "model-simulations_m4.csv"))
df.coef4 <- read.csv(here("output_model", "model-coefficients_m4.csv"))

df.pred5ab <- read.csv(here("output_model", "model-predictions_m5_direct.ab.csv"))
df.sim5ab  <- read.csv(here("output_model", "model-simulations_m5_direct.ab.csv"))
df.coef5ab <- read.csv(here("output_model", "model-coefficients_m5_direct.ab.csv"))

df.pred5in <- read.csv(here("output_model", "model-predictions_m5_direct.in.csv"))
df.sim5in  <- read.csv(here("output_model", "model-simulations_m5_direct.in.csv"))
df.coef5in <- read.csv(here("output_model", "model-coefficients_m5_direct.in.csv"))

df.pred5ml <- read.csv(here("output_model", "model-predictions_m5_omni.ml.csv"))
df.sim5ml  <- read.csv(here("output_model", "model-simulations_m5_omni.ml.csv"))
df.coef5ml <- read.csv(here("output_model", "model-coefficients_m5_omni.ml.csv"))

df.pred5oab <- read.csv(here("output_model", "model-predictions_m5_omni.ab.csv"))
df.sim5oab  <- read.csv(here("output_model", "model-simulations_m5_omni.ab.csv"))
df.coef5oab <- read.csv(here("output_model", "model-coefficients_m5_omni.ab.csv"))

df.pred5ab$meth <- df.sim5ab$meth <- df.coef5ab$meth <- "direct ab"
df.pred5in$meth <- df.sim5in$meth <- df.coef5in$meth <- "direct an"
df.pred5ml$meth <- df.sim5ml$meth <- df.coef5ml$meth <- "omni ml"
df.pred5oab$meth <- df.sim5oab$meth <- df.coef5oab$meth <- "omni ab"

df.pred5 <- do.call(bind_rows, list(df.pred5ab, df.pred5in, df.pred5ml, df.pred5oab))
df.sim5 <- do.call(bind_rows, list(df.sim5ab, df.sim5in, df.sim5ml, df.sim5oab))
df.coef5 <- do.call(bind_rows, list(df.coef5ab, df.coef5in, df.coef5ml, df.coef5oab))

## model/method validation
df.meth <- data.frame(meth = c("direct.ab", "direct.ab", "direct.in", "direct.in", "omni.ab", "omni.ml"),
                      site = c("maisC", "maisD", "maisC", "maisD", "maisC", "maisC"),
                      R2 = NA)
df.meth$ID <- paste0(df.meth$site, "_", df.meth$meth)

Lk <- Lm <- Ld <- list()
for(i in df.meth$ID) {
  Lk[[i]] <- read.csv(here("output_model", paste0("kfold_m5_", i, ".csv"))) 
  Lm[[i]] <- readRDS( here("output_model", paste0("model_m5_", i, ".RDS"))) 
  Ld[[i]] <- read.csv(here("output_model", paste0("data_m5_", i, ".csv"))) 
}

# fileL <- list.files(here("output_model"), pattern = "data_m2")
# df.m2 <- data.frame()
# for(f in fileL) df.m2 <- rbind(df.m2, read.csv(here("output_model", f)))

fileL <- list.files(here("output_model"), pattern = "data_m4")
df.m4 <- data.frame()
for(f in fileL) df.m4 <- rbind(df.m4, read.csv(here("output_model", f)))

fileL <- list.files(here("output_model"), pattern = "data_test_m5")
df.t5 <- data.frame()
for(f in fileL) df.t5 <- rbind(df.t5, read.csv(here("output_model", f)))

#### 1) plot sites -------------------------------------------------------------

## using ggspatial (better)

gC <- ggplot() + 
  annotation_map_tile(type = "osm") +
  geom_sf(data = isoC, alpha = 0.3, aes(color = z)) +
  scale_color_viridis_c("elevation \n(m asl)", direction = -1, option = "rocket",
                        guide = guide_colorbar(order = 1)) +
  new_scale_color() +
  geom_sf(data = df.stat[df.stat$station.project_id == "maisC",], 
          aes(pch = type, fill = type), color = "black", stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("direct" = 23, "omni" = 21),
                     guide = guide_legend(order = 2)) +
  scale_fill_manual("station type", values = c("direct" = "white", "omni" = "grey60"),
                    guide = guide_legend(order = 2)) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimC3[1], dimC3[3]), ylim = c(dimC3[2], dimC3[4]), expand = F) +
  facet_grid(~station.project_id) +
  theme_void(base_size = 15) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 25, colour="white")) +
  guides(shape = guide_legend(position = "right"),
         fill = guide_legend(position = "right"))

gD <- ggplot() + 
  annotation_map_tile(type = "osm") +
  geom_sf(data = isoD, alpha = 0.3, aes(color = z)) +
  scale_color_viridis_c("elevation \n(m asl)", direction = -1, option = "rocket",
                        guide = guide_colorbar(order = 1)) +
  new_scale_color() +
  geom_sf(data = df.stat[df.stat$station.project_id == "maisD",], 
          aes(pch = type, fill = type), color = "black", stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("direct" = 23, "omni" = 21),
                     guide = guide_legend(order = 2)) +
  scale_fill_manual("station type", values = c("direct" = "white", "omni" = "grey60"),
                    guide = guide_legend(order = 2)) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimD3[1], dimD3[3]), ylim = c(dimD3[2], dimD3[4]), expand = F) +
  facet_grid(~station.project_id) +
  theme_void(base_size = 15) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 25, colour="white")) +
  guides(shape = "none",
         fill = "none")

gA <- ggplot() + 
  # annotation_map_tile(type = "osm") +
  # geom_sf(data = isoA, alpha = 0.3, aes(color = z)) +
  # scale_color_viridis_c("elevation \n(m asl)", direction = -1, option = "rocket",
  #                       guide = guide_colorbar(order = 1)) +
  geom_sf(data = shp.ger, color = "grey20", fill = "grey60") +
  geom_sf(data = st_centroid(st_union(df.stat)), 
                size = 7) +
  # geom_sf_label(data = st_centroid(st_union(df.stat)), 
  #               label = "study area", size = 7) +
  # geom_sf_label(data = st_centroid(st_union(df.stat[df.stat$station.project_id == "maisC",])), 
  #         label = "maisC", size = 7) +
  # geom_sf_label(data = st_centroid(st_union(df.stat[df.stat$station.project_id == "maisD",])), 
  #         label = "maisD", size = 7) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimA[1], dimA[3]), ylim = c(dimA[2], dimA[4]), expand = F) +
  theme_void(base_size = 15) +
  # facet_wrap(~"study area") +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 25, colour="white"))

# g <- gA | (gC/gD)  # + plot_layout(guides = 'collect')
g <- gA | gC | gD  # + plot_layout(guides = 'collect')

# ggsave("./plots/plotSite.pdf", plot = g, width = 30, height = 18, device = "pdf", units = "cm")
ggsave("./plots/plotSite.pdf", plot = g, width = 30, height = 12, device = "pdf", units = "cm")

## get min distance between stations
directC <- df.stat[df.stat$station.project_id == "maisD" & df.stat$type == "direct",]
directC <- st_transform(directC, crs = crs)

dist <- st_distance(directC, directC)
dist <- as.data.frame(t(dist))
dist <- as.data.frame(lapply(dist, units::drop_units))
dist[dist == 0] <- NA
dist$distance_m <- apply(dist, 1, FUN = min, na.rm = T)
min(dist$distance_m); mean(dist$distance_m); max(dist$distance_m)


#### 2) testtracks -------------------------------------------------------------
## plot testtracks and generate table
shp.l$length_km <- round(as.vector(st_length(shp.l))/1000, 2)

shp.p$time <- round(shp.p$time, 0)
shp.p <- shp.p[second(shp.p$time) %% 2 != 0, ]

# summarize shp.p
sum.p <- shp.p %>% group_by(date, area, site) %>%
  summarise("N points" = n(),
            "duration_h" = round(as.numeric(max(time) - min(time), units = "hours"), 2), 
            .groups = "drop")
sum.p$trackID <- paste0(sum.p$area, c(1, "test", 2, 3, 1, 2, 3, "test"))

sum.p <- left_join(st_drop_geometry(sum.p), 
                   st_drop_geometry(shp.l), by = c("date", "area", "site"))

## add trackID to shp.p
shp.p <- left_join(shp.p, sum.p[, c("date", "trackID")], by = "date")

sum.p <- sum.p[, c("trackID", "date", "duration_h", "length_km", "N points", "site")]
sum.p <- sum.p[order(sum.p$trackID),]
sum.p[, "N tags"] <- c(3, 4, 4, 4, 3, 3, 3, 3)

write.csv(sum.p, "./data/tabTesttracks.csv", row.names = F)            


gC <- ggplot() + 
  annotation_map_tile(type = "osm") +
  geom_sf(data = shp.p[shp.p$site == "maisC",], 
          alpha = 0.2, aes(color = trackID), pch = 1) +
  scale_color_viridis_d("track", option = "viridis", begin = 0.1, end = 0.9,
                        guide = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2, stroke = 1))) +
  new_scale_color() +
  geom_sf(data = df.stat[df.stat$station.project_id == "maisC",], 
          aes(pch = type, fill = type), color = "black", stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("direct" = 23, "omni" = 21),
                     guide = guide_legend(order = 2)) +
  scale_fill_manual("station type", values = c("direct" = "white", "omni" = "grey60"),
                     guide = guide_legend(order = 2)) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimC[1], dimC[3]), ylim = c(dimC[2], dimC[4]), expand = F) +
  facet_grid(~"maisC") +
  theme_void(base_size = 15) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 30, colour="white"))

gD <- ggplot() + 
  annotation_map_tile(type = "osm") +
  geom_sf(data = shp.p[shp.p$site == "maisD",], 
          alpha = 0.2, aes(color = trackID), pch = 1) +
  scale_color_viridis_d("track", option = "viridis", begin = 0.1, end = 0.9,
                        guide = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2, stroke = 1))) +
  new_scale_color() +
  geom_sf(data = df.stat[df.stat$station.project_id == "maisD",], 
          aes(pch = type, fill = type), color = "black", stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("direct" = 23, "omni" = 21),
                     guide = guide_legend(order = 2)) +
  scale_fill_manual("station type", values = c("direct" = "white", "omni" = "grey60"),
                    guide = guide_legend(order = 2)) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimD[1], dimD[3]), ylim = c(dimD[2], dimD[4]), expand = F) +
  facet_grid(~"maisD") +
  theme_void(base_size = 15) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 30, colour="white"))

g <- gC + gD  # + plot_layout(guides = 'collect')

g
ggsave("./plots/plotTrack.pdf", plot = g, width = 30, height = 15, device = "pdf", units = "cm")

#### 4) nextmatch lastmatch
df.filter <- data.frame(signal = c("s0", "s1", "s2", NA, NA, NA),
                        signal2 = c(0, 1, 2, NA, NA, NA),
                        nm = c(NA, "(s0 nextmatch)", "(s1 nextmatch)", NA, NA, NA),
                        y = 150.100,
                        x = c(0, 1.1, 2.15, 0.7, 1.3, 2.65),
                        color = c("darkblue", "darkorange", "grey30", "grey", "grey", "grey"))
texp <- 1 ## expected time interval
x0 <- df.filter$x[which(df.filter$signal == "s0")]
c0 <- df.filter$color[which(df.filter$signal == "s0")]
x1 <- df.filter$x[which(df.filter$signal == "s1")]
c1 <- df.filter$color[which(df.filter$signal == "s1")]
x2 <- df.filter$x[which(df.filter$signal == "s2")]

g <- ggplot() + 
  ## area for s0
  geom_line(aes(x = x0+texp, y = c(150.050, 150.200)), color = c0, lty = "dashed") +
  geom_rect(aes(xmin = x0+texp-0.5*texp, xmax = x0+texp+0.5*texp, ymin = 150.050, ymax = 150.200),
            fill = c0, alpha = 0.2) +
  geom_text(aes(x = x0+texp, y = 150.205), label = "window to search for s0 nextmatch", color = c0, size = 5) +
  ## area for s1
  geom_line(aes(x = x1+texp, y = c(150.050, 150.200)), color = c1, lty = "dashed") +
  geom_rect(aes(xmin = x1+texp-0.5*texp, xmax = x1+texp+0.5*texp, ymin = 150.050, ymax = 150.200),
            fill = c1, alpha = 0.2) +
  geom_text(aes(x = x1+texp, y = 150.205), label = "window to search for s1 nextmatch", color = c1, size = 5) +
  
  ## lines s0
  geom_line(aes(x = c(x0, x0+texp), y = 150.125), lwd = 1) +
  geom_text(aes(x = mean(c(x0, x0+texp)), y = 150.13), label = expression(italic(paste("t"["expect"]))), size = 5) +

  geom_line(aes(x = c(x0, x1), y = 150.135), color = c0, lwd = 1) +
  geom_text(aes(x = mean(c(x0, x1)), y = 150.14), color = c0, label = expression(italic(paste("t"["actual"]))), size = 5) +
  
  geom_line(aes(x = c(x1, x0+texp), y = 150.145), color = c0, lwd = 1) +
  geom_text(aes(x = mean(c(x1, x0+texp)), y = 150.15), color = c0, label = "Nextmatch delta", size = 5) +
  
  ## lines s1
  geom_line(aes(x = c(x1, x1+texp), y = 150.125), lwd = 1) +
  geom_text(aes(x = mean(c(x1, x1+texp)), y = 150.13), label = expression(italic(paste("t"["expect"]))), size = 5) +
  
  geom_line(aes(x = c(x1, x2), y = 150.135), color = c1, lwd = 1) +
  geom_text(aes(x = mean(c(x1, x2)), y = 150.14), color = c1, label = expression(italic(paste("t"["actual"]))), size = 5) +
  
  geom_line(aes(x = c(x2, x1+texp), y = 150.145), color = c1, lwd = 1) +
  geom_text(aes(x = mean(c(x2, x1+texp)), y = 150.15), color = c1, label = "Nextmatch delta", size = 5) +
  
  ## points
  geom_point(aes(x = x, y = y), size = 5, color = df.filter$color, data = df.filter) +
  geom_text(aes(label = signal, x = x, y = y-0.005), data = df.filter, size = 5, color = df.filter$color) +
  geom_text(aes(label = nm, x = x, y = y-0.01), data = df.filter[!is.na(df.filter$nm),], size = 5, color = c(c0, c1)) +
  
  xlab("time [s]") +
  ylab("frequency [MHz]") +
  ylim(150.05, 150.21) +
  theme_light(base_size = 15)
  
ggsave("./plots/plotNextmatch.pdf", plot = g, width = 30, height = 18, device = "pdf", units = "cm")


#### 4) plot results -----------------------------------------------------------
## m1-4

## add column to indicate whether position was estimated for all methods
df.m4 <- df.m4 %>% group_by(X_time, tagID, site) %>%
  mutate(Nmeth = length(unique(meth)),
         nP = n()) %>%
  ungroup()

df.m4$allM <- "no"
df.m4$allM[df.m4$site == "maisC" & df.m4$Nmeth == 4] <- "yes"
df.m4$allM[df.m4$site == "maisD" & df.m4$Nmeth == 2] <- "yes"

df.m4$allM <- factor(df.m4$allM, levels = c("no", "yes"))
df.m4$allM2 <- factor(df.m4$allM, levels = c("yes", "no"))
df.m4$allM.site <- factor(paste0(df.m4$allM, " ", df.m4$site), levels = c("yes maisC", "no maisC", "yes maisD", "no maisD"))

## add column to indicate whether position is kept due to thresholds
df.m4$thresh <- "yes"
df.m4$thresh[df.m4$meth == "direct.ab" & df.m4$site == "maisC" & df.m4$Ac == df.m4$Sc] <- "no"
df.m4$thresh[df.m4$meth == "direct.ab" & df.m4$site == "maisD" & df.m4$Ac <= 3 & df.m4$Sc %in% c(1, 2, 3)] <- "no"
df.m4$thresh[df.m4$meth == "omni.ab" & df.m4$site == "maisC" & df.m4$Ac == 1] <- "no"

## add several IDs for plotting and merging
df.m4$group <- paste0(df.m4$site, " \n", df.m4$meth)
df.m4$ID <- paste0(df.m4$site, "_", df.m4$meth)
df.m4$meth2 <- df.m4$meth
df.m4$meth2[df.m4$meth2 == "direct.ab"] <- "direct ab"
df.m4$meth2[df.m4$meth2 == "direct.in"] <- "direct an"
df.m4$meth2[df.m4$meth2 == "omni.ab"] <- "omni ab"
df.m4$meth2[df.m4$meth2 == "omni.ml"] <- "omni ml"
df.m4$meth3 <- df.m4$meth
df.m4$meth3[df.m4$meth3 == "direct.ab"] <- "direct \nab"
df.m4$meth3[df.m4$meth3 == "direct.in"] <- "direct \nan"
df.m4$meth3[df.m4$meth3 == "omni.ab"] <- "omni \nab"
df.m4$meth3[df.m4$meth3 == "omni.ml"] <- "omni \nml"
df.m4$group2 <- paste0(df.m4$meth3, " \n", df.m4$site)

# get n per group
df.m4 <- df.m4 %>% group_by(site, allM2, meth) %>%
  mutate(N = n()) %>% ungroup()

df.sim1$detR <- factor(df.sim1$detR, ordered = T)

tmp.m1 <- subset(df.sim1, model == "m1" & meth == "direct.ab")
#tmp.m3 <- subset(df.sim, model == "m3" & meth == "ab_ql")
tmp.m4 <- subset(df.sim4, model == "m4")
tmp.m4$meth2 <- tmp.m4$meth
tmp.m4$meth2[tmp.m4$meth2 == "direct.ab"] <- "direct ab"
tmp.m4$meth2[tmp.m4$meth2 == "direct.in"] <- "direct an"
tmp.m4$meth2[tmp.m4$meth2 == "omni.ab"] <- "omni ab"
tmp.m4$meth2[tmp.m4$meth2 == "omni.ml"] <- "omni ml"
tmp.m4$meth3 <- tmp.m4$meth
tmp.m4$meth3[tmp.m4$meth3 == "direct.ab"] <- "direct \nab"
tmp.m4$meth3[tmp.m4$meth3 == "direct.in"] <- "direct \nan"
tmp.m4$meth3[tmp.m4$meth3 == "omni.ab"] <- "omni \nab"
tmp.m4$meth3[tmp.m4$meth3 == "omni.ml"] <- "omni \nml"
tmp.m4$group2 <- paste0(tmp.m4$meth3, " \n", tmp.m4$site)
#tmp.m2 <- subset(df.pred, model == "m2" & meth == "ab_ql")

df.pred5ab$diffCI2 <- ifelse(df.pred5ab$diffCI <= 10, "0-10",
                          ifelse(df.pred5ab$diffCI <= 25, "10-25",
                                ifelse(df.pred5ab$diffCI <= 50, "25-50",
                                       ifelse(df.pred5ab$diffCI <= 75, "50-75", paste0("75-", floor(max(df.pred5ab$diffCI)))))))
df.pred5ab$diffCIq502 <- ifelse(df.pred5ab$diffCIq50 <= 10, "0-10",
                          ifelse(df.pred5ab$diffCIq50 <= 25, "10-25",
                                 ifelse(df.pred5ab$diffCIq50 <= 50, "25-50",
                                        ifelse(df.pred5ab$diffCIq50 <= 75, "50-75", paste0("75-", floor(max(df.pred5ab$diffCIq50)))))))
df.pred5ab$diffCIq652 <- ifelse(df.pred5ab$diffCIq65 <= 10, "0-10",
                            ifelse(df.pred5ab$diffCIq65 <= 25, "10-25",
                                   ifelse(df.pred5ab$diffCIq65 <= 50, "25-50",
                                          ifelse(df.pred5ab$diffCIq65 <= 75, "50-75", paste0("75-", floor(max(df.pred5ab$diffCIq65)))))))

df.pred5in$diffCI2 <- ifelse(df.pred5in$diffCI <= 10, "0-10",
                             ifelse(df.pred5in$diffCI <= 25, "10-25",
                                    ifelse(df.pred5in$diffCI <= 50, "25-50",
                                           ifelse(df.pred5in$diffCI <= 75, "50-75", paste0("75-", floor(max(df.pred5in$diffCI)))))))
df.pred5in$diffCIq502 <- ifelse(df.pred5in$diffCIq50 <= 10, "0-10",
                                ifelse(df.pred5in$diffCIq50 <= 25, "10-25",
                                       ifelse(df.pred5in$diffCIq50 <= 50, "25-50",
                                              ifelse(df.pred5in$diffCIq50 <= 75, "50-75", paste0("75-", floor(max(df.pred5in$diffCIq50)))))))
df.pred5in$diffCIq652 <- ifelse(df.pred5in$diffCIq65 <= 10, "0-10",
                                ifelse(df.pred5in$diffCIq65 <= 25, "10-25",
                                       ifelse(df.pred5in$diffCIq65 <= 50, "25-50",
                                              ifelse(df.pred5in$diffCIq65 <= 75, "50-75", paste0("75-", floor(max(df.pred5in$diffCIq65)))))))

## get prop Points per method
df.points <- df.m4 %>% group_by(site, meth, group2, allM2) %>%
  summarise(nP = n(), .groups = "drop")
df.points <- left_join(df.points,
                     sum.p %>% group_by(site) %>% 
                       summarize(nPtot = sum(`N points` * `N tags`), .groups = "drop"), 
                     by = c("site"))

df.points$propP <- round(df.points$nP/df.points$nPtot*100, 0)

## model output horizontal combined
g1 <- ggplot(tmp.m1) +
  stat_halfeye(aes(x = detR, y = sim.m, group = site, color = site, fill = site,
                   linewidth = after_stat(.width)), # needed for linewidth
               .width = c(0.5, 0.95),
               #color = "black",
               fatten_point = 3,
               point_color = "white",
               normalize = "groups", scale = 0.8, # "groups" ?
               slab_alpha = 0.45, side = "left", pch = 24) +
  scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
  scale_color_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_fill_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  xlab("detection range [m]") +
  ylab("mean pPE [m]") +
  # ylim(15, 100) +
  ylim(35, 128) +
  theme_light(base_size = 16) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank())

g4 <- ggplot() + 
  stat_halfeye(aes(x = meth3, y = sim.m, group = site, color = site, fill = site,
                   linewidth = after_stat(.width)), # needed for linewidth
               .width = c(0.5, 0.95),
               #color = "black",
               fatten_point = 3,
               point_color = "white",
               normalize = "groups", scale = 0.8, # "groups" ?
               slab_alpha = 0.6, side = "left",
               data = tmp.m4, pch = 24) +
  scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
  scale_color_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_fill_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  xlab("method") +
  ylab("mean pPE [m]") +
  # ylim(15, 100) +
  ylim(35, 128) +
  
  ## legend
  geom_rect(aes(xmin = 0.6, xmax = 4.2, ymin = 112, ymax = 128), color = "grey", fill = "white") +
  geom_text(aes(x = c(0.9, 2.4), y = 125, label = c("site", "all meth.")), size = 6,
            hjust = 0) +
  geom_text(aes(x = c(0.9), y = c(120, 115), label = c("maisC", "maisD")), color = c("#440154", "#23898e"), size = 5,
            hjust = 0) +
  geom_text(aes(x = c(2.4), y = c(120, 115), label = c("yes", "yes")), color = c("#440154", "#23898e"), size = 5,
            hjust = 0) +
  geom_text(aes(x = c(3.3), y = c(120, 115), label = c("no", "no")), color = c("#440154", "#23898e"), alpha = 0.5, size = 5,
            hjust = 0) +
  
  theme_light(base_size = 16) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank())

g4c <- ggplot(data = df.m4, aes(x = group2, y = cover, fill = site, color = site, alpha = allM2)) +
  stat_eye(.width = NA, adjust = 3, 
               #point_interval = NULL,
               aes(side = ifelse(allM2 == "yes", "left", "right"),
                   thickness = after_stat(pdf*n)), # scales to counts
           point_color = "white", pch = 21,
           normalize = "all", scale = 0.6) +
  
  scale_color_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_fill_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_alpha_manual("all meth.", values = c("yes" = 0.6, "no" = 0.3)) +
  ylim(0, max(df.m4$cover)) +
  xlab("method") +
  ylab("station cover") +
  theme_light(base_size = 16) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())


g4PE <- ggplot() +
  stat_eye(.width = NA, adjust = 3, 
           #point_interval = NULL,
           data = df.m4, aes(x = group2, y = PE, fill = site, color = site, alpha = allM2,
                             side = ifelse(allM2 == "yes", "left", "right"),
                             thickness = after_stat(pdf*n)), # scales to counts
           point_color = "white", pch = 21,
           normalize = "all", scale = 0.6)  +  
  scale_color_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_fill_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_alpha_manual("all meth.", values = c("yes" = 0.6, "no" = 0.3)) +
  xlab("method") +
  ylab("PE [m]") +
  geom_text(data = df.points[df.points$allM2 == "no",], 
            aes(x = group2, y = 0.4, label = paste0(propP, "%"), color = site), alpha = 0.5) +
  geom_text(data = df.points[df.points$allM2 == "yes",], 
            aes(x = group2, y = 0.7, label = paste0(propP, "%"), color = site), alpha = 0.8) +
  scale_y_continuous(trans = "log10", breaks = c(5, 10, 50, 100, 500, 1000), limits = c(0.4, NA)) +
  scale_x_discrete(labels = c("direct \nab", "direct \nab", "direct \nan", "direct \nan", "omni \nab", "omni \nml")) +
  theme_light(base_size = 16) +
  theme(legend.position = "none")

layout <- "
AABBCCCC
AABBDDDD
"

g <- g1 + g4 + g4c + g4PE + plot_layout(design = layout)
g
ggsave("./plots/plotModel_horizontal.pdf", plot = g, width = 30, height = 15, device = "pdf", units = "cm")

# g <- g1 + g4
# g
# ggsave("./plots/plotModel_horizontal.pdf", plot = g, width = 30, height = 18, device = "pdf", units = "cm")

g5 <- ggplot(df.pred5ab) +
  geom_point(aes(y = Sc, x = Ac, color = pred_mod, fill = pred_mod, size = diffCI2), pch = 22) +
  scale_color_viridis_c("mean \npPE [m]", option = "rocket", limits = c(1, NA), na.value = "#FAEBDDFF") +
  scale_fill_viridis_c("mean \npPE [m]", option = "rocket", limits = c(1, NA), na.value = "#FAEBDDFF") +
  scale_size_discrete("range CI") +
  xlab("number of antennas") +
  ylab("number of stations") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))

g5.2 <- ggplot(df.pred5ab) + 
  # geom_point(aes(x = Ac, y = PE, group = as.factor(Sc), color = as.factor(Sc)),
  #            data = df.m4[df.m4$meth2 == "direct \nab",],
  #            size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.3)) +
  geom_line(aes(x = Ac, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.5, position = position_dodge(width = 0.4)) + 
  geom_linerange(aes(x = Ac, ymin = lwr, ymax = upr,
                      group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 0.5, alpha = 0.7, position = position_dodge(width = 0.4)) + 
  geom_linerange(aes(x = Ac, ymin = lwr25, ymax = upr75,
                     group = as.factor(Sc), color = as.factor(Sc)),
                 lwd = 1, alpha = 0.7, position = position_dodge(width = 0.4)) + 
  # geom_point(aes(x = Ac, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
  #                size = 1, alpha = 0.5, position = position_dodge(width = 0.3)) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("mean pPE [m]") +
  # ylim(0, 500) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))

# ggplot(df.m4[df.m4$meth3 == "direct \nan",]) + geom_histogram(aes(x = PE, fill = site)) + facet_wrap(~Sc)

ggsave("./plots/plotPA.pdf", plot = g5, width = 34, height = 12, device = "pdf", units = "cm")
ggsave("./plots/plotPA2.pdf", plot = g5.2, width = 34, height = 12, device = "pdf", units = "cm")

## add dummy sites to have violins one same side
# tmp.sim <- unique(df.sim5[, c("site", "meth", "Ac")])
# tmp.sim$Sc <- 0
# df.sim5 <- bind_rows(df.sim5, tmp.sim)

df.sim5$Sc <- as.factor(df.sim5$Sc)
df.pred5$Sc <- as.factor(df.pred5$Sc)

## correlations between coefficients

gPE <- ggplot(df.pred5ab) + 
  geom_point(aes(x = Ac, y = PE, group = as.factor(Sc), color = as.factor(Sc)),
             data = df.m4[df.m4$meth3 == "direct \nab",],
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Ac, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.7) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("position error (raw and pred.)") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))
ggsave("./plots/plotAcSc_PE.pdf", plot = gPE, width = 34, height = 12, device = "pdf", units = "cm")

gcov <- ggplot(df.pred5ab) + 
  geom_point(aes(x = Ac, y = cover, group = as.factor(Sc), color = as.factor(Sc)),
             data = df.m4[df.m4$meth3 == "direct \nab",],
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Ac, y = cover, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.7) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("station cover") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))
ggsave("./plots/plotAcSc_cover.pdf", plot = gcov, width = 34, height = 12, device = "pdf", units = "cm")

gmaxS <- ggplot(df.pred5ab) + 
  geom_point(aes(x = Ac, y = maxSig, group = as.factor(Sc), color = as.factor(Sc)),
             data = df.m4[df.m4$meth3 == "direct \nab",],
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Ac, y = maxSig, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.7) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("maximum recorded signal") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))
ggsave("./plots/plotAcSc_maxSig.pdf", plot = gmaxS, width = 34, height = 12, device = "pdf", units = "cm")


gwgt <- ggplot(df.pred5ab) + 
  geom_point(aes(x = Ac, y = Weight, group = as.factor(Sc), color = as.factor(Sc)),
             data = df.m4[df.m4$meth3 == "direct \nab",],
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Ac, y = Weight, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.7) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("weight (sum of normalized signals") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))
ggsave("./plots/plotAcSc_weight.pdf", plot = gwgt, width = 34, height = 12, device = "pdf", units = "cm")


gcw <- ggplot(df.m4[df.m4$meth3 == "direct \nab",]) + 
  geom_point(aes(x = cover, y = Weight, group = as.factor(Sc), color = as.factor(Sc)),
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("station cover") +
  ylab("weight (sum of normalized signals)") +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))
ggsave("./plots/plotcover_weight.pdf", plot = gcw, width = 34, height = 12, device = "pdf", units = "cm")

gsw <- ggplot(df.m4[df.m4$meth3 == "direct \nab",]) + 
  geom_point(aes(x = maxSig, y = Weight, group = as.factor(Sc), color = as.factor(Sc)),
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("maximum recorded signal") +
  ylab("weight (sum of normalized signals)") +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))

ggsave("./plots/plotmaxSig_weight.pdf", plot = gsw, width = 34, height = 12, device = "pdf", units = "cm")

gsc <- ggplot(df.m4[df.m4$meth3 == "direct \nab",]) + 
  geom_point(aes(x = maxSig, y = cover, group = as.factor(Sc), color = as.factor(Sc)),
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("maximum recorded signal") +
  ylab("station cover") +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))

ggsave("./plots/plotmaxSig_cover.pdf", plot = gsc, width = 34, height = 12, device = "pdf", units = "cm")

#### predictive performance test testtrack -------------------------------------
## check whether X_time is present in all methods
df.t5 <- df.t5 %>% group_by(X_time, tagID, site) %>% 
  mutate(allM2 = ifelse(site == "maisC" & n() == 4, "yes", 
                          ifelse(site == "maisD" & n() == 2, "yes", "no"))) %>% 
  ungroup()

## add site_meth groups
df.t5$group <- paste0(df.t5$site, "_", df.t5$meth)
df.t5$group2 <- paste0(df.t5$meth, "_", df.t5$site)

## add column to indicate whether position is kept due to thresholds
df.t5$thresh <- "yes"
df.t5$thresh[df.t5$meth == "direct.ab" & df.t5$site == "maisC" & df.t5$Ac == df.t5$Sc] <- "no"
df.t5$thresh[df.t5$meth == "direct.ab" & df.t5$site == "maisD" & df.t5$Ac <= 3 & df.t5$Sc %in% c(1, 2, 3)] <- "no"
df.t5$thresh[df.t5$meth == "omni.ab" & df.t5$site == "maisC" & df.t5$Ac <= 2] <- "no"
df.t5$thresh[df.t5$meth == "omni.ml" & df.t5$site == "maisC" & df.t5$Ac <= 2] <- "no"

df.t5 <- st_as_sf(df.t5, coords = c("lon", "lat"), crs = crsLL)

## get core area (convex hull of stations)
coreCd <- st_convex_hull(st_union(df.stat[df.stat$station.project_id == "maisC" & df.stat$type == "direct",]))
coreCo <- st_convex_hull(st_union(df.stat[df.stat$station.project_id == "maisC" & df.stat$type == "omni",]))
coreDd <- st_convex_hull(st_union(df.stat[df.stat$station.project_id == "maisD" & df.stat$type == "direct",]))

df.t5 <- st_transform(df.t5, crs = crsPlot) 

df.t5 <- df.t5 %>%
  mutate(keepCd = lengths(st_within(df.t5, coreCd)) > 0,
         keepCo = lengths(st_within(df.t5, coreCo)) > 0,
         keepDd = lengths(st_within(df.t5, coreDd)) > 0)

df.t5$core <- NA
df.t5$core[df.t5$meth %in% c("direct.ab", "direct.in") & df.t5$site == "maisC"] <- df.t5$keepCd[df.t5$meth %in% c("direct.ab", "direct.in") & df.t5$site == "maisC"]
df.t5$core[df.t5$meth %in% c("omni.ab", "omni.ml") & df.t5$site == "maisC"] <- df.t5$keepCo[df.t5$meth %in% c("omni.ab", "omni.ml") & df.t5$site == "maisC"]
df.t5$core[df.t5$meth %in% c("direct.ab", "direct.in") & df.t5$site == "maisD"] <- df.t5$keepDd[df.t5$meth %in% c("direct.ab", "direct.in") & df.t5$site == "maisD"]

g1 <- ggplot() +
  geom_hline(yintercept = 0, lty = "dashed", color = "grey") +
  stat_eye(.width = NA, adjust = 3, 
           point_interval = NULL,
           mapping = aes(x = group2, y = pred_mod, fill = site, color = site, alpha = allM2,
                         side = "right",
               thickness = after_stat(pdf*n)), # scales to counts
           normalize = "all", scale = 0.6,
           point_color = "white",
           data = df.t5, pch = 24) +
  
  stat_eye(.width = NA, adjust = 3, 
           point_interval = NULL,
           mapping = aes(x = group2, y = PE, fill = site, color = site, alpha = allM2,
                         side = "left",
                         thickness = after_stat(pdf*n)), # scales to counts
           normalize = "all", scale = 0.6,
           point_color = "white",
           data = df.t5, pch = 21) +
  
  stat_pointinterval(.width = NA,
                     mapping = aes(x = group2, y = pred_mod, fill = site, color = site, alpha = allM2,
                                   side = "right"),
                     point_color = "white", pch = 24,
                     data = df.t5) +
  
  stat_pointinterval(.width = NA,
                     mapping = aes(x = group2, y = PE, fill = site, color = site, alpha = allM2,
                                   side = "left"),
                     point_color = "white", pch = 21,
                     data = df.t5) +
  
  scale_color_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_fill_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_alpha_manual("all meth.", values = c("yes" = 0.6, "no" = 0.3)) +
  #ylim(-300, 500) +
  scale_x_discrete(labels = c("direct \nab", "direct \nab", "direct \nan", "direct \nan", "omni \nab", "omni \nml")) +
  scale_y_continuous(trans = "log10", breaks = c(5, 10, 50, 100, 500, 1000), limits = c(1, NA)) +
  xlab("method") +
  ylab("PE (left), pPE (right) [m]") +
  
  ## legend
  geom_rect(aes(xmin = 1.25, xmax = 3.75, ymin = 600, ymax = 2500), color = "grey", fill = "white") +
  geom_text(aes(x = c(1.5, 2.5), y = 1900, label = c("site", "all meth.")), size = 6,
            hjust = 0) +
  geom_text(aes(x = c(1.5), y = c(850, 1200), label = c("maisC", "maisD")), color = c("#440154", "#23898e"), size = 5,
            hjust = 0) +
  geom_text(aes(x = c(2.5), y = c(850, 1200), label = c("yes", "yes")), color = c("#440154", "#23898e"), size = 5,
            hjust = 0) +
  geom_text(aes(x = c(3), y = c(850, 1200), label = c("no", "no")), color = c("#440154", "#23898e"), alpha = 0.4, size = 5,
            hjust = 0) +
  
  theme_light(base_size = 16) +
  theme(legend.position = "none"); g1

g2 <- ggplot() +
  geom_hline(yintercept = 0, lty = "dashed", color = "grey") +
  stat_eye(.width = NA, adjust = 3, 
           #point_interval = NULL,
           mapping = aes(x = group2, y = diff, fill = site, color = site, alpha = allM2,
                         side = ifelse(allM2 == "yes", "left", "right"),
                         thickness = after_stat(pdf*n)), # scales to counts
           point_color = "white", pch = 21,
           normalize = "all", scale = 0.6,
           data = df.t5) +
  
  scale_color_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_fill_manual("site", values = c("maisC" = "#440154", "maisD" = "#23898e")) +
  scale_alpha_manual("all meth.", values = c("yes" = 0.7, "no" = 0.3)) +
  ylim(-300, 500) +
  scale_x_discrete(labels = c("direct \nab", "direct \nab", "direct \nan", "direct \nan", "omni \nab", "omni \nml")) +
  xlab("method") +
  ylab("PE - pPE [m]") +
  
  # ## legend
  # geom_rect(aes(xmin = 2, xmax = 5, ymin = 300, ymax = 500), color = "grey", fill = "white") +
  # geom_text(aes(x = c(2.5, 3.5), y = 450, label = c("site", "all meth.")), size = 6,
  #           hjust = 0) +
  # geom_text(aes(x = c(2.5), y = c(400, 350), label = c("maisC", "maisD")), color = c("#440154", "#23898e"), size = 5,
  #           hjust = 0) +
  # geom_text(aes(x = c(3.5), y = c(400, 350), label = c("yes", "yes")), color = c("#440154", "#23898e"), size = 5,
  #           hjust = 0) +
  # geom_text(aes(x = c(4), y = c(400, 350), label = c("no", "no")), color = c("#440154", "#23898e"), alpha = 0.4, size = 5,
  #           hjust = 0) +
  
  theme_light(base_size = 16) +
  theme(legend.position = "none")

g <- g1 + g2
ggsave("./plots/plotTest.pdf", plot = g, width = 34, height = 12, device = "pdf", units = "cm")
ggsave("./plots/plotTest.pdf", plot = g1, width = 18, height = 13, device = "pdf", units = "cm")


df.t5 <- as.data.frame(df.t5)
df.t5$Sc <- as.factor(df.t5$Sc)

## summarize df.t5
sum.t5 <- df.t5 %>% group_by(site, meth, Sc, Ac) %>%
  summarise(nPoints = n(),
            mdiff = median(diff),
            lwr = quantile(diff, probs = 0.025),
            upr = quantile(diff, probs = 0.975),
            lwr25 = quantile(diff, probs = 0.25),
            upr75 = quantile(diff, probs = 0.75), .groups = "drop"
  )

sum.t5 <- as.data.frame(sum.t5)

sum.t5.2 <- df.t5 %>% group_by(site, meth, Sc, Ac) %>%
  summarise(nPoints = n(),
            mPE = median(pred_mod),
            lwr = quantile(pred_mod, probs = 0.025),
            upr = quantile(pred_mod, probs = 0.975),
            lwr25 = quantile(pred_mod, probs = 0.25),
            upr75 = quantile(pred_mod, probs = 0.75), .groups = "drop"
  )

sum.t5.2 <- as.data.frame(sum.t5.2)

gdiff <- ggplot(sum.t5[sum.t5$meth == "direct.ab",]) + 

  geom_hline(yintercept = 0, lty = "dashed", color = "grey40") +
  geom_line(aes(x = Ac, y = mdiff, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.5, position = position_dodge(width = 0.4)) + 
  geom_linerange(aes(x = Ac, ymin = lwr, ymax = upr,
                     group = as.factor(Sc), color = as.factor(Sc)),
                 lwd = 0.5, alpha = 0.7, position = position_dodge(width = 0.4)) + 
  geom_linerange(aes(x = Ac, ymin = lwr25, ymax = upr75,
                     group = as.factor(Sc), color = as.factor(Sc)),
                 lwd = 1, alpha = 0.7, position = position_dodge(width = 0.4)) + 
  
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("PE - pPE [m]") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32), limits = c(0, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

gpPE <- ggplot(sum.t5.2[sum.t5.2$meth == "direct.ab",]) + 
  # geom_point(aes(x = Ac, y = PE, group = as.factor(Sc), color = as.factor(Sc)),
  #            data = df.m4[df.m4$meth2 == "direct \nab",],
  #            size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.3)) +
  geom_line(aes(x = Ac, y = mPE, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.5, position = position_dodge(width = 0.4)) + 
  geom_linerange(aes(x = Ac, ymin = lwr, ymax = upr,
                     group = as.factor(Sc), color = as.factor(Sc)),
                 lwd = 0.5, alpha = 0.7, position = position_dodge(width = 0.4)) + 
  geom_linerange(aes(x = Ac, ymin = lwr25, ymax = upr75,
                     group = as.factor(Sc), color = as.factor(Sc)),
                 lwd = 1, alpha = 0.7, position = position_dodge(width = 0.4)) + 
  # geom_point(aes(x = Ac, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
  #                size = 1, alpha = 0.5, position = position_dodge(width = 0.3)) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("pPE [m]") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32), limits = c(0, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"),
        axis.title.x = element_blank(),
        legend.position = c(0.9, 0.65),
        legend.background = element_rect(fill = "white",
                                         colour = "grey")) +
  guides(color = guide_legend(ncol = 2))

layout <- "
          A
          B
          "
g <- gpPE + gdiff + plot_layout(design = layout)
ggsave("./plots/plotPEdiff.pdf", plot = g, width = 34, height = 20, device = "pdf", units = "cm")



####5) table kfold -------------------------------------------------------------
## function for pseudo R2
pR2 <- function(model, resp) {
  1 - sum((resp - predict(model, type = "response"))^2)/
    sum((resp - mean(resp))^2)
  
}

df.meth <- data.frame(meth = c("direct.ab", "direct.ab", "direct.in", "direct.in", "omni.ab", "omni.ml"),
                      site = c("maisC", "maisD", "maisC", "maisD", "maisC", "maisC"),
                      R2 = NA)
df.meth$ID <- paste0(df.meth$site, "_", df.meth$meth)

# for(i in df.meth$ID) {
  # df.meth$R2[df.meth$ID == i] <- round(pR2(model = Lm[[i]], resp = Lm[[i]]$frame$PE)*100, 2)
  # df.meth$MAE[df.meth$ID == i] <- round(mean(Lk[[i]]$mfull2), 2)
  # Ld[[i]]$ID <- paste0(Ld[[i]]$site, "_", Ld[[i]]$meth)
  # Ld[[i]] <- left_join(Ld[[i]], df.m4[df.m4$ID == i, c("X_time", "tagID", "allM", "thresh", "ID")], by = c("X_time", "tagID", "ID"))
  # df.meth$MAE2[df.meth$ID == i] <- mean(abs(predict(Lm[[i]], newdata = Ld[[i]], type = "response")-Ld[[i]]$PE))
  # df.meth$MAE3[df.meth$ID == i] <- mean(abs(predict(Lm[[i]], newdata = Ld[[i]][Ld[[i]]$allM == "yes",], type = "response")-Ld[[i]]$PE[Ld[[i]]$allM == "yes"]))
  # df.meth$MAE4[df.meth$ID == i] <- mean(abs(predict(Lm[[i]], newdata = Ld[[i]][Ld[[i]]$thresh == "yes",], type = "response")-Ld[[i]]$PE[Ld[[i]]$thresh == "yes"]))
  # df.meth$ME2[df.meth$ID == i] <- -mean(predict(Lm[[i]], newdata = Ld[[i]], type = "response")-Ld[[i]]$PE)
  # df.meth$ME3[df.meth$ID == i] <- -mean(predict(Lm[[i]], newdata = Ld[[i]][Ld[[i]]$allM == "yes",], type = "response")-Ld[[i]]$PE[Ld[[i]]$allM == "yes"])
  # df.meth$ME4[df.meth$ID == i] <- -mean(predict(Lm[[i]], newdata = Ld[[i]][Ld[[i]]$thresh == "yes",], type = "response")-Ld[[i]]$PE[Ld[[i]]$thresh == "yes"])
  # df.meth$m2[df.meth$ID == i] <- median(predict(Lm[[i]], newdata = Ld[[i]], type = "response"))
  # df.meth$m3[df.meth$ID == i] <- median(predict(Lm[[i]], newdata = Ld[[i]][Ld[[i]]$allM == "yes",], type = "response"))
  # df.meth$m4[df.meth$ID == i] <- median(predict(Lm[[i]], newdata = Ld[[i]][Ld[[i]]$thresh == "yes",], type = "response"))
# }

df.t5 %>% group_by(group) %>% 
  summarize(MAE = mean(abs(diff)), MAEq50 = mean(abs(diffq50)), MAEq65 = mean(abs(diffq65)), 
            meanPA = mean(pred_mod), meanPA50 = mean(q50), meanPA65 = mean(q65),
            meanPE = mean(PE), meanPE50 = median(PE), meanPE65 = quantile(PE, probs = 0.65))

df.meth <- left_join(df.meth,
                     df.t5 %>% group_by(group) %>% 
                       summarize(MAE = mean(abs(diff)), 
                                 meanPA = mean(pred_mod),
                                 meanPE = mean(PE),
                                 meanDiff = mean(diff)),
                     by = c("ID" = "group"))

df.meth <- left_join(df.meth,
                     df.t5[df.t5$allM2 == "yes",] %>% group_by(group) %>% 
                       summarize(MAE2 = mean(abs(diff)), 
                                 meanPA2 = mean(pred_mod),
                                 meanPE2 = mean(PE),
                                 meanDiff2 = mean(diff)),
                     by = c("ID" = "group"))

df.meth <- left_join(df.meth,
                     df.t5[df.t5$thresh == "yes",] %>% group_by(group) %>% 
                       summarize(MAE3 = mean(abs(diff)), 
                                 meanPA3 = mean(pred_mod),
                                 meanPE3 = mean(PE),
                                 meanDiff3 = mean(diff)),
                     by = c("ID" = "group"))

df.meth <- left_join(df.meth,
                     df.t5[df.t5$core == TRUE,] %>% group_by(group) %>% 
                       summarize(MAE4 = mean(abs(diff)), 
                                 meanPA4 = mean(pred_mod),
                                 meanPE4 = mean(PE),
                                 meanDiff4 = mean(diff)),
                     by = c("ID" = "group"))

df.meth <- left_join(df.meth,
                     df.t5 %>% group_by(site, meth) %>% summarize(nPfull = n(), .groups = "drop"), 
                 by = c("site", "meth"))

df.meth <- left_join(df.meth,
                     df.t5[df.t5$allM2 == "yes",] %>% group_by(site, meth) %>% 
                       summarize(nPallM = n(), .groups = "drop"), 
                     by = c("site", "meth"))

df.meth <- left_join(df.meth,
                     df.t5[df.t5$thresh == "yes",] %>% group_by(site, meth) %>% 
                       summarize(nPthresh = n(), .groups = "drop"), 
                     by = c("site", "meth"))

df.meth <- left_join(df.meth,
                     df.t5[df.t5$core == T,] %>% group_by(site, meth) %>% 
                       summarize(nPcore = n(), .groups = "drop"), 
                     by = c("site", "meth"))

df.meth <- left_join(df.meth,
                     sum.p[sum.p$trackID %in% c("Ctest", "Dtest"),] %>% group_by(site) %>% 
                       summarize(nPtot = sum(`N points` * `N tags`), .groups = "drop"), 
                     by = c("site"))

df.meth$propPfull <- round(df.meth$nPfull/df.meth$nPtot*100, 0)
df.meth$propPallM <- round(df.meth$nPallM/df.meth$nPtot*100, 0)
df.meth$propPthresh <- round(df.meth$nPthresh/df.meth$nPtot*100, 0)
df.meth$propPcore <- round(df.meth$nPcore/df.meth$nPtot*100, 0)

colL <- c("meanPE", "meanPE2", "meanPE3", "meanPE4", "meanPA", "meanPA2", "meanPA3", "meanPA4", 
          "MAE", "MAE2", "MAE3", "MAE4", "meanDiff", "meanDiff2", "meanDiff3", "meanDiff4",
          "propPfull", "propPallM", "propPthresh", "propPcore")
colL2 <- c("meanPE",  "meanPA",  "MAE",  "meanDiff",  "propPfull", 
           "meanPE2", "meanPA2", "MAE2", "meanDiff2", "propPallM", 
           "meanPE3", "meanPA3", "MAE3", "meanDiff3", "propPthresh",
           "meanPE4", "meanPA4", "MAE4", "meanDiff4", "propPcore")

df.meth[, colL] <- round(df.meth[, colL], 0)

df.meth <- df.meth[order(df.meth$site),]

df.meth <- df.meth[, c("site", "meth", colL2)]

write.csv(df.meth, "./plots/tabMeth.csv", row.names = F)


## plot test data --------------------------------------------------------------
#dsn = "../data/data_test/GT13_R2_S1_Test.antennabeams.gpkg"
# df.t1 <- st_read(dsn = "../data/data_test/GT07_R1_S1.antennabeams.gpkg", layer = "PA")
# df.t2 <- st_read(dsn = "../data/data_test/GT07_R1_S3.antennabeams.gpkg", layer = "PA")
# df.t1 <- st_read(dsn = "../data/data_test/GT11_S_S1.antennabeams.gpkg", layer = "PA")
# df.t2 <- st_read(dsn = "../data/data_test/GT11_S_S2.antennabeams.gpkg", layer = "PA")
# df.t3 <- st_read(dsn = "../data/data_test/GT11_S_S3.antennabeams.gpkg", layer = "PA")
# df.t <- rbind(df.t1, df.t2)
# df.t <- rbind(df.t, df.t3)
# 
# print(ggplot(df.t[df.t$Ac > df.t$Sc,]) + 
#         annotation_map_tile(type = "osm") +
#         geom_sf(aes(color = X_time), alpha = 0.2, pch = 1, size = 5) +
#         geom_sf(data = df.stat[df.stat$station.project_id == "maisC" & df.stat$type == "direct",], 
#                 color = "black", fill = "grey", stroke = 1.2, size = 5, alpha = 0.7, aes(shape = "21")) +
#         scale_shape_manual(name = "Station", values = c(21), labels = c(" ")) +
#         # geom_sf(aes(fill = date_start, color = date_start, shape = "23"), alpha = 0.5, size = 10, color = "black",
#         #         data = shp.tmp) +
#         # scale_size_continuous("PA est. points",range = c(2, 10)) +
#         # coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
#         annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
#         scale_color_viridis_c("Zeit", option = "inferno", end = 0.8, begin = 0.1,
#                               # labels = c("08:00", "12:00", "16:00"), breaks = c(1697004000, 1697018400, 1697032800)) +
#                               labels = c("08:00", "12:00", "16:00"), breaks = c(1697090400, 1697104800, 1697119200)) +
#   # scale_fill_viridis_c("time", option = "inferno", end = 0.8, begin = 0.1,
#         #                      limits = c(min(df.t$X_time, shp.tmp$date_start), max(df.t$X_time, shp.tmp$date_start))) +
#         theme_void(base_size = 15))

df.t1 <- st_read(dsn = "../data/data_test/GT11_S_S1.antennabeams.gpkg", layer = "PA")
df.t2 <- st_read(dsn = "../data/data_test/GT11_S_S2.antennabeams.gpkg", layer = "PA")
df.t3 <- st_read(dsn = "../data/data_test/GT11_S_S3.antennabeams.gpkg", layer = "PA")
df.t <- rbind(df.t1, df.t2)
df.t <- rbind(df.t, df.t3)
df.t1 <- st_read(dsn = "../data/data_test/RO08_R1_S1.antennabeams.gpkg", layer = "PA")
df.t2 <- st_read(dsn = "../data/data_test/RO08_R1_S2.antennabeams.gpkg", layer = "PA")
df.t3 <- st_read(dsn = "../data/data_test/RO08_R1_S3.antennabeams.gpkg", layer = "PA")
df.tr <- rbind(df.t1, df.t2)
df.tr <- rbind(df.tr, df.t3)

df.t$time <- as.numeric(format(df.t$X_time, "%H")) + as.numeric(format(df.t$X_time, "%M"))/60
df.tr$time <- as.numeric(format(df.tr$X_time, "%H")) + as.numeric(format(df.tr$X_time, "%M"))/60

df.bird <- read.csv("../data/data_test/GTdata_full_bird.csv")
df.bird$date_start <- as.POSIXct(df.bird$date_start, tz = "UTC")
df.bird$date_start <- with_tz(df.bird$date_start, "CET")
df.bird$date_stop <- as.POSIXct(df.bird$date_stop, tz = "UTC")
df.bird$date_stop <- with_tz(df.bird$date_stop, "CET")
colnames(df.bird)[colnames(df.bird) %in% c("lon", "lat")] <- c("lon.true", "lat.true")

## summarize df.bird by ID and transform to sf
df.bird <- df.bird %>% group_by(lon.true, lat.true, ID, Testtag) %>%
  summarize(date_start = mean(date_start), .groups = "drop")
df.bird <- st_as_sf(df.bird, coords = c("lat.true", "lon.true"), crs = crsLL)

shp.tmp <- df.bird[df.bird$Testtag == unique(df.t$Individual) &
                     df.bird$date_start >= min(df.t$X_time) &
                     df.bird$date_start <= max(df.t$X_time),]
shp.tmpr <- df.bird[df.bird$Testtag == unique(df.tr$Individual) &
                       df.bird$date_start >= min(df.tr$X_time) &
                       df.bird$date_start <= max(df.tr$X_time),]
shp.tmp$time <- as.numeric(format(shp.tmp$date_start, "%H")) + as.numeric(format(shp.tmp$date_start, "%M"))/60
shp.tmpr$time <- as.numeric(format(shp.tmpr$date_start, "%H")) + as.numeric(format(shp.tmpr$date_start, "%M"))/60


g2 <- ggplot() + 
  geom_sf(aes(shape = "B"), alpha = 0.7, size = 4, stroke = 1, data = df.t[1,]) +
  geom_sf(aes(shape = "C"), alpha = 0.7, size = 4, stroke = 1, data = df.tr[1,]) +        
  annotation_map_tile(type = "osm") +
  # est. data
  geom_sf(aes(color = time, shape = "B", size = PA), alpha = 0.2, df.t[df.t$Ac > df.t$Sc,]) +
  geom_sf(aes(color = time, shape = "C", size = PA), alpha = 0.2, df.tr[df.tr$Ac > df.tr$Sc,]) +
  # handheld data
  geom_sf(aes(fill = time), shape = 21, alpha = 0.6, size = 5, color = "black",
          data = shp.tmp, stroke = 1) +
  geom_sf(aes(fill = time), shape = 24, alpha = 0.6, size = 5, color = "black",
          data = shp.tmpr, stroke = 1) +
  # stations
  geom_sf(data = df.stat[df.stat$station.project_id == "maisC" & df.stat$type == "direct",], 
          color = "black", fill = "white", stroke = 1.2, size = 5, alpha = 0.7, aes(shape = "A")) +
  scale_shape_manual(name = " ", values = c(23, 1, 2), labels = c("directional \nstation", "Great Tit", "European \nRobin")) +
  scale_size_continuous(name = "pPE", range = c(1, 7)) +
  # scale_size_continuous("PA est. points",range = c(2, 10)) +
  # coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1.5, width_hint = 0.4, bar_cols = c("grey60", "white"), location = "br") +
  scale_color_viridis_c("time", option = "inferno", end = 0.8, begin = 0.1, direction = -1,
                       labels = c("08:00 am", "12:00 pm", "04:00 pm"), breaks = c(8, 12, 16), limits = c(6.9, 18.5)) +
  scale_fill_viridis_c("time", option = "inferno", end = 0.8, begin = 0.1, direction = -1,
                        labels = c("08:00 am", "12:00 pm", "04:00 pm"), breaks = c(8, 12, 16), limits = c(6.9, 18.5)) +
  theme_void(base_size = 20)

# ggsave("./plots/plotAnimal1.pdf", plot = g1, width = 34, height = 18, device = "pdf", units = "cm")
ggsave("./plots/plotAnimal.pdf", plot = g2, width = 34, height = 18, device = "pdf", units = "cm")

## run markdown for supplement -------------------------------------------------
rmarkdown::render(input = here("supplement.Rmd"))

