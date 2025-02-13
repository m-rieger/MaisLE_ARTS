## plots for paper

rm(list = ls())

## packages
library(tidyverse)
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

## variables
crs <- 31467 # Gauss Krüger in m
crsPlot <- 3857

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

## model/method validation
df.meth <- data.frame(meth = c("direct.ab", "direct.ab", "direct.in", "direct.in", "omni.ab", "omni.ml"),
                      site = c("maisC", "maisD", "maisC", "maisD", "maisC", "maisC"),
                      R2 = NA, MAE = NA)
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

#### 1) plot sites -------------------------------------------------------------
## using basemaps
# set defaults for the basemap
# set_defaults(map_service = "osm", map_type = "topographic")
# 
# ggplot() + 
#   basemap_gglayer(st_buffer(df.stat[df.stat$station.project_id == "maisC",], dist = 800)) +
#   scale_fill_identity() + 
#   geom_sf(data = df.stat[df.stat$station.project_id == "maisC",], 
#           aes(pch = station.type, color = station.type), stroke = 1.2) +
#   coord_sf() +
#   scale_shape_manual("station type", values = c(3, 16)) +
#   scale_color_manual("station type", values = c("black", "grey50")) +
#   # annotation_scale(width_hint = 0.2, bar_cols = c("grey60", "white"), color = "grey60") +
#   theme_classic()
# 
# ggplot() + 
#   basemap_gglayer(st_buffer(df.stat[df.stat$station.project_id == "maisD",], dist = 800)) +
#   scale_fill_identity() + 
#   geom_sf(data = df.stat[df.stat$station.project_id == "maisD",], 
#           aes(pch = station.type, color = station.type), stroke = 1.2) +
#   coord_sf() +
#   scale_shape_manual("station type", values = c(3, 16)) +
#   scale_color_manual("station type", values = c("black", "grey50")) +
#   # annotation_scale(width_hint = 0.2, bar_cols = c("grey60", "white"), color = "grey60") +
#   theme_classic()

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
         fill = guide_legend(position = "right")); gC

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

g <- grid.arrange(gC, gD, nrow = 1)

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

#### 2) testtracks -------------------------------------------------------------
## plot testtracks and generate table
shp.l$length_km <- round(as.vector(st_length(shp.l))/1000, 2)

shp.p <- shp.p[second(shp.p$time) %% 2 != 0, ]

# summarize shp.p
sum.p <- shp.p %>% group_by(date, area, site) %>%
  summarise("N points" = n(),
            "duration_h" = round(as.numeric(max(time) - min(time), units = "hours"), 2), 
            .groups = "drop")
sum.p$trackID <- paste0(sum.p$area, c(1, 2, 3, 4, 1, 2, 3, 4))

sum.p <- left_join(st_drop_geometry(sum.p), 
                   st_drop_geometry(shp.l), by = c("date", "area", "site"))

## add trackID to shp.p
shp.p <- left_join(shp.p, sum.p[, c("date", "trackID")], by = "date")

sum.p <- sum.p[, c("site", "trackID", "date", "duration_h", "length_km", "N points")]
sum.p <- sum.p[order(sum.p$trackID),]
sum.p[, "N tag"] <- c(3, 4, 4, 4, 3, 3, 3, 3)

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
gC

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
gD

g <- gC + gD  # + plot_layout(guides = 'collect')

g
ggsave("./plots/plotTrack.pdf", plot = g, width = 30, height = 15, device = "pdf", units = "cm")

#### 3) raster (station cover) -------------------------------------------------


## plot raster (= merged polygons)
g <- ggplot(shp.r) +
  annotation_map_tile(type = "osm") +
  geom_sf(aes(fill = dens, color = dens), lwd = 0.1) +
  scale_fill_viridis_c("station cover", direction = -1, option = "rocket", na.value = NA, limits = c(0.1, NA)) +
  scale_color_viridis_c("station cover", direction = -1, option = "rocket", na.value = NA, limits = c(0.1, NA)) +

  new_scale_fill() +
  
  geom_sf(data = df.stat[df.stat$station.project_id == "maisC" & df.stat$type == "direct",], 
          aes(pch = type, fill = type), color = "black", stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("direct" = 23, "omni" = 21),
                     guide = guide_legend(order = 2)) +
  scale_fill_manual("station type", values = c("direct" = "white", "omni" = "grey60"),
                    guide = guide_legend(order = 2)) +
  # geom_sf(data = df[100000:100100,], pch = 1, color = "grey80", alpha = 0.5) + # check whether dimensions are correct
  coord_sf(xlim = c(dimC2[1], dimC2[3]), ylim = c(dimC2[2], dimC2[4]), expand = F) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  theme_void(base_size = 15) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 30, colour="white"))

g
ggsave("./plots/plotCover.pdf", plot = g, width = 18, height = 15, device = "pdf", units = "cm")

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

df.m4$allM <- factor(df.m4$allM, levels = c("yes", "no"))

## add column to indicate whether position is kept due to thresholds
df.m4$thresh <- "yes"
df.m4$thresh[df.m4$meth == "direct.ab" & df.m4$site == "maisC" & df.m4$Ac == df.m4$Sc] <- "no"
df.m4$thresh[df.m4$meth == "direct.ab" & df.m4$site == "maisD" & df.m4$Ac <= 3 & df.m4$Sc %in% c(1, 2, 3)] <- "no"
df.m4$thresh[df.m4$meth == "omni.ab" & df.m4$site == "maisC" & df.m4$Ac == 1] <- "no"

## add several IDs for plotting and merging
df.m4$group <- paste0(df.m4$site, " \n", df.m4$meth)
df.m4$ID <- paste0(df.m4$site, "_", df.m4$meth)
df.m4$meth2 <- df.m4$meth
df.m4$meth2[df.m4$meth2 == "direct.ab"] <- "direct \nab"
df.m4$meth2[df.m4$meth2 == "direct.in"] <- "direct \nin"
df.m4$meth2[df.m4$meth2 == "omni.ab"] <- "omni \nab"
df.m4$meth2[df.m4$meth2 == "omni.ml"] <- "omni \nml"

df.sim1$detR <- factor(df.sim1$detR, ordered = T)

tmp.m1 <- subset(df.sim1, model == "m1" & meth == "direct.ab")
#tmp.m3 <- subset(df.sim, model == "m3" & meth == "ab_ql")
tmp.m4 <- subset(df.sim4, model == "m4")
tmp.m4$meth2 <- tmp.m4$meth
tmp.m4$meth2[tmp.m4$meth2 == "direct.ab"] <- "direct \nab"
tmp.m4$meth2[tmp.m4$meth2 == "direct.in"] <- "direct \nin"
tmp.m4$meth2[tmp.m4$meth2 == "omni.ab"] <- "omni \nab"
tmp.m4$meth2[tmp.m4$meth2 == "omni.ml"] <- "omni \nml"
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


g1 <- ggplot(tmp.m1) +
  stat_halfeye(aes(x = detR, y = sim.m,
                   linewidth = after_stat(.width)), # needed for linewidth
               .width = c(0.5, 0.95),
               #color = "black",
               fatten_point = 3,
               normalize = "panels", scale = 0.7) +
  scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
 # scale_color_viridis_d("site", end = 0.5) +
  xlab("detection range (m)") +
  ylab("est. mean PE (m)") +
  ylim(30, 130) +
  #ggtitle("m1") +
  facet_wrap(~site, ncol = 1, scales = "free_y", strip.position = "left") +
  theme_light(base_size = 16) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank())


# g3 <- ggplot(tmp.m3) + 
#   stat_halfeye(aes(x = pos, y = sim,
#                    linewidth = after_stat(.width), color = site), # needed for linewidth
#                .width = c(0.5, 0.95),
#                #color = "black",
#                fatten_point = 3,
#                normalize = "panels", scale = 0.7, strip.position = "right") +
#   scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
#   scale_color_viridis_d("site", end = 0.5) +
#   xlab("position") +
#   ylab("est. mean PE (m)") +
#   #ggtitle("m3") +
#   facet_wrap(~site, ncol = 1, scales = "free_y", strip.position = "right") +
#   theme_light(base_size = 14) +
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         strip.text = element_blank(),
#         strip.background = element_blank()
#   )

g4 <- ggplot(tmp.m4) + 
  stat_halfeye(aes(x = meth2, y = sim.m,
                   linewidth = after_stat(.width)), # needed for linewidth
               .width = c(0.5, 0.95),
               #color = "black",
               fatten_point = 3,
               normalize = "panels", scale = 0.7) +
  scale_linewidth_continuous(range = c(15, 5)) + # Define range of linewidths (reverse!!)
  # scale_color_viridis_d("site", end = 0.5) +
  xlab("method") +
  ylab("est. mean PE (m)") +
  ylim(30, 130) +
  #ggtitle("m4") +
  facet_wrap(~site, ncol = 1, scales = "free_y", strip.position = "right") +
  theme_light(base_size = 16) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank())

g4c <- ggplot(data = df.m4, aes(x = meth2, y = cover, fill = allM, color = allM)) +
  geom_split_violin(alpha = 0.4, trim = TRUE, scale = "count", width = 1.1) +
  scale_color_viridis_d("include", end = 0.5) +
  scale_fill_viridis_d("include", end = 0.5) +
  facet_wrap(~site, ncol = 1, strip.position = "right") +
  ylim(0, max(df.m4$cover)) +
  #coord_flip() +
  xlab("method") +
  theme_light(base_size = 16) +
  theme(legend.position = c(0.85, 0.25),
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))

g5 <- ggplot(df.pred5ab) +
  geom_point(aes(y = Sc, x = Ac, color = pred_mod, fill = pred_mod, size = diffCI2), pch = 22) +
  scale_color_viridis_c("est. mean \nPE (m)", option = "rocket", limits = c(0, NA), na.value = "#FAEBDDFF") +
  scale_fill_viridis_c("est. mean \nPE (m)", option = "rocket", limits = c(0, NA), na.value = "#FAEBDDFF") +
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
            lwd = 1, alpha = 0.5, position = position_dodge(width = 0.3)) + 
  geom_linerange(aes(x = Ac, y = pred_mod, ymin = lwr, ymax = upr,
                      group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.5, position = position_dodge(width = 0.3)) + 
  # geom_point(aes(x = Ac, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
  #                size = 1, alpha = 0.5, position = position_dodge(width = 0.3)) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("est. mean PE (m) and 95% CI") +
  # ylim(0, 500) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))

g5in <- ggplot(df.pred5in) +
  geom_point(aes(y = Sc, x = Ac, color = pred_mod, size = diffCI2)) +
  scale_color_viridis_c("est. mean \nPE (m)", option = "rocket", limits = c(0, NA), na.value = "#FAEBDDFF") +
  scale_size_discrete("range CI") +
  xlab("number of antennas") +
  ylab("number of stations") +
  facet_wrap(~site, ncol = 1, strip.position = "right") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10)) +
  theme_light(base_size = 16)

ggplot(df.m4[df.m4$meth2 == "direct \nin",]) + geom_histogram(aes(x = PE, fill = site)) + facet_wrap(~Sc)
g5in.2 <- ggplot(df.m4[df.m4$meth2 == "direct \nin",]) + 
  # geom_split_violin(aes(x = as.factor(Sc), y = PE, fill = site, color = site),
  #                   alpha = 0.4, trim = F, scale = "area") +
  stat_halfeye(aes(x = as.factor(Sc), y = PE, fill = site, thickness = after_stat(pdf*n)),
                    alpha = 0.4, scale = 0.7) +
  geom_linerange(aes(x = as.factor(Sc), y = pred_mod, ymin = lwr, ymax = upr,
                     group = as.factor(Sc), color = as.factor(Sc)),
                 lwd = 1, alpha = 0.5, data = df.pred5in,
                 position = position_dodge(width = 0.5)) +
  # geom_point(aes(x = Ac+0.2, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
  #           size = 3, alpha = 0.5) + 
  # geom_point(aes(x = Ac, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
  #                size = 1, alpha = 0.5, position = position_dodge(width = 0.3)) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  scale_fill_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("est. mean PE (m) and 95% CI") +
  # ylim(0, 500) +
  # scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white")); g5in.2

g5in.2 <- ggplot(df.m4[df.m4$meth2 == "direct \nab",]) + 
  # geom_split_violin(aes(x = as.factor(Sc), y = PE, fill = site, color = site),
  #                   alpha = 0.4, trim = F, scale = "area") +
  stat_halfeye(aes(x = as.factor(Ac), y = PE, thickness = after_stat(pdf*n)),
               alpha = 0.4, scale = 0.7) +
  geom_linerange(aes(x = as.factor(Ac), y = pred_mod, ymin = lwr, ymax = upr,
                     group = as.factor(Sc), color = as.factor(Sc)),
                 lwd = 1, alpha = 0.5, data = df.pred5ab,
                 position = position_dodge(width = 0.5)) +
  # geom_point(aes(x = Ac+0.2, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
  #           size = 3, alpha = 0.5) + 
  # geom_point(aes(x = Ac, y = pred_mod, group = as.factor(Sc), color = as.factor(Sc)),
  #                size = 1, alpha = 0.5, position = position_dodge(width = 0.3)) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  scale_fill_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("est. mean PE (m) and 95% CI") +
  # ylim(0, 500) +
  # scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16, 18)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white")); g5in.2


layout <- "
AABCC
DDDDD
"

layout <- "
ABCC
ABCC
"

g <- g1 + g4 + g4c + plot_layout(design = layout)
g
ggsave("./plots/plotModel.pdf", plot = g, width = 30, height = 18, device = "pdf", units = "cm")

ggsave("./plots/plotPA.pdf", plot = g5, width = 34, height = 12, device = "pdf", units = "cm")
ggsave("./plots/plotPA2.pdf", plot = g5.2, width = 34, height = 12, device = "pdf", units = "cm")

gPE <- ggplot(df.pred5ab) + 
  geom_point(aes(x = Ac, y = PE, group = as.factor(Sc), color = as.factor(Sc)),
             data = df.m4[df.m4$meth2 == "direct \nab",],
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
             data = df.m4[df.m4$meth2 == "direct \nab",],
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
             data = df.m4[df.m4$meth2 == "direct \nab",],
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
             data = df.m4[df.m4$meth2 == "direct \nab",],
             size = 0.5, alpha = 0.2, pch = 1, position = position_dodge(width = 0.5)) +
  geom_line(aes(x = Ac, y = Weight, group = as.factor(Sc), color = as.factor(Sc)),
            lwd = 1, alpha = 0.7) + 
  scale_color_viridis_d("Sc", option = "viridis") +
  facet_wrap(~site, ncol = 2, strip.position = "top") +
  xlab("number of antennas") +
  ylab("weight (sum of normalized signals)") +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20, 24, 28, 32)) +
  theme_light(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 22, colour="white"))
ggsave("./plots/plotAcSc_weight.pdf", plot = gwgt, width = 34, height = 12, device = "pdf", units = "cm")


gcw <- ggplot(df.m4[df.m4$meth2 == "direct \nab",]) + 
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

gsw <- ggplot(df.m4[df.m4$meth2 == "direct \nab",]) + 
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

gsc <- ggplot(df.m4[df.m4$meth2 == "direct \nab",]) + 
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

####5) table kfold -------------------------------------------------------------
pR2 <- function(model, resp) {
  1 - sum((resp - predict(model, type = "response"))^2)/
    sum((resp - mean(resp))^2)
  
}

for(i in df.meth$ID) {
  df.meth$R2[df.meth$ID == i] <- round(pR2(model = Lm[[i]], resp = Lm[[i]]$frame$PE)*100, 2)
  df.meth$MAE[df.meth$ID == i] <- round(mean(Lk[[i]]$mfull2), 2)
  Ld[[i]]$ID <- paste0(Ld[[i]]$site, "_", Ld[[i]]$meth)
  Ld[[i]] <- left_join(Ld[[i]], df.m4[df.m4$ID == i, c("X_time", "tagID", "allM", "ID")], by = c("X_time", "tagID", "ID"))
  df.meth$MAE2[df.meth$ID == i] <- mean(abs(predict(Lm[[i]], newdata = Ld[[i]], type = "response")-Ld[[i]]$PE))
  df.meth$MAE3[df.meth$ID == i] <- mean(abs(predict(Lm[[i]], newdata = Ld[[i]][Ld[[i]]$allM == "yes",], type = "response")-Ld[[i]]$PE[Ld[[i]]$allM == "yes"]))
}

df.meth <- left_join(df.meth,
                 df.m4 %>% group_by(site, meth) %>% summarize(nPfull = n(), .groups = "drop"), 
                 by = c("site", "meth"))

df.meth <- left_join(df.meth,
                     df.m4[df.m4$allM == "yes",] %>% group_by(site, meth) %>% summarize(nP = n(), .groups = "drop"), 
                     by = c("site", "meth"))

df.meth <- left_join(df.meth,
                     sum.p %>% group_by(site) %>% summarize(nPtot = sum(`N points` * `N tag`), .groups = "drop"), 
                     by = c("site"))

df.meth$propPfull <- round(df.meth$nPfull/df.meth$nPtot*100, 2)
df.meth$propP <- round(df.meth$nP/df.meth$nPtot*100, 2)

df.meth[, c("MAE2", "MAE3", "MAE")] <- round(df.meth[, c("MAE2", "MAE3", "MAE")], 0)

write.csv(df.meth, "./plots/tabMeth.csv", row.names = F)


### convex hull for each method ------------------------------------------------
dist_C <- 0.5*800
dist_D <- 0.5*900

tmp <- st_transform(df.stat[df.stat$station.project_id == "maisC" & df.stat$type == "direct",], crs = crs)
pointCd <- tmp %>%
  mutate(geometry = st_geometry(.) + c(dist_C, 0))
pointCd <- rbind(pointCd,
                 tmp %>% mutate(geometry = st_geometry(.) + c(-dist_C, 0)))
pointCd <- rbind(pointCd,
                 tmp %>% mutate(geometry = st_geometry(.) + c(0, dist_C)))
pointCd <- rbind(pointCd,
                 tmp %>% mutate(geometry = st_geometry(.) + c(0, -dist_C)))
pointCd <- st_set_crs(pointCd, st_crs(tmp))

polyCd <- st_convex_hull(st_union(pointCd))


pointCo <- df.stat[df.stat$station.project_id == "maisC" & df.stat$type == "omni",]
polyCo <- st_convex_hull(st_union(pointCo))


gC <- ggplot() + 
  annotation_map_tile(type = "osm") +
  geom_sf(data = df.stat[df.stat$station.project_id == "maisC",], 
          aes(pch = type, fill = type), color = "black", stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("direct" = 23, "omni" = 21),
                     guide = guide_legend(order = 2)) +
  scale_fill_manual("station type", values = c("direct" = "white", "omni" = "grey60"),
                    guide = guide_legend(order = 2)) +
  
  geom_sf(data = polyCo, alpha = 0.2, fill = "grey60") +
  geom_sf(data = polyCd, alpha = 0.2, fill = "white") +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimC[1], dimC[3]), ylim = c(dimC[2], dimC[4]), expand = F) +
  facet_grid(~"maisC") +
  theme_void(base_size = 15) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 30, colour="white"))
gC
