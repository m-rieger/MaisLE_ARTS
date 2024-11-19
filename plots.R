## plots for paper
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

## variables
crs <- 31467 # Gauss Krüger in m
crsPlot <- 3857

#### 0) load data ####
## station
df.stat <- read.csv("./data/station.csv", encoding = "latin1")
df.stat <- st_as_sf(x = df.stat[df.stat$analysis == "yes",],                         
                    coords = c("station.lon", "station.lat"),
                    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    remove = F)
df.stat$plot <- "station" # to add the shape to legend
df.stat$type <- "mono"
df.stat$type[df.stat$station.type == "VHF Directional Station"] <- "quadro"
df.stat <- st_transform(df.stat, crs = crsPlot) 

## get dimensions for plotting
## get plotting dimensions
dimC <- st_bbox(st_buffer(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisC",])), dist = 1300))
dimD <- st_bbox(st_buffer(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisD",])), dist = 1300))
dimA <- st_bbox(st_buffer(df.stat, dist = 3000))

## two study area
df.area <- rbind(st_centroid(st_union(df.stat[df.stat$station.project_id == "maisC",])),
                 st_centroid(st_union(df.stat[df.stat$station.project_id == "maisD",])))

## isolines
shp.iso <- read_sf("../data/Isolines_BB.shp")
shp.iso <- st_zm(st_transform(shp.iso, crs = crsPlot)) # change crs and XYM to XY 
isoC <- st_intersection(shp.iso, st_as_sfc(dimC))
isoD <- st_intersection(shp.iso, st_as_sfc(dimD))
isoA <- st_intersection(shp.iso, st_as_sfc(dimA))
isoA <- isoA[isoA$z %in% seq(0, 100, 5),]

#### 1) plot sites ####
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
          aes(pch = type, color = type), stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("quadro" = 3, "mono" = 1),
                     guide = guide_legend(order = 2)) +
  scale_color_manual("station type", values = c("quadro" = "black", "mono" = "grey30"),
                     guide = guide_legend(order = 2)) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimC[1], dimC[3]), ylim = c(dimC[2], dimC[4]), expand = F) +
  facet_grid(~station.project_id) +
  theme_void(base_size = 15) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 30, colour="white"))


gD <- ggplot() + 
  annotation_map_tile(type = "osm") +
  geom_sf(data = isoD, alpha = 0.3, aes(color = z)) +
  scale_color_viridis_c("elevation \n(m asl)", direction = -1, option = "rocket",
                        guide = guide_colorbar(order = 1)) +
  new_scale_color() +
  geom_sf(data = df.stat[df.stat$station.project_id == "maisD",], 
          aes(pch = type, color = type), stroke = 1.2, size = 3, alpha = 0.7) +
  scale_shape_manual("station type", values = c("quadro" = 3, "mono" = 1),
                     guide = guide_legend(order = 2)) +
  scale_color_manual("station type", values = c("quadro" = "black", "mono" = "grey30"),
                     guide = guide_legend(order = 2)) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimD[1], dimD[3]), ylim = c(dimD[2], dimD[4]), expand = F) +
  facet_grid(~station.project_id) +
  theme_void(base_size = 15) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 30, colour="white"))

g <- grid.arrange(gC, gD, nrow = 1)

gA <- ggplot() + 
  annotation_map_tile(type = "osm") +
  geom_sf(data = isoA, alpha = 0.3, aes(color = z)) +
  scale_color_viridis_c("elevation \n(m asl)", direction = -1, option = "rocket",
                        guide = guide_colorbar(order = 1)) +
  geom_sf_label(data = st_centroid(st_union(df.stat[df.stat$station.project_id == "maisC",])), 
          label = "maisC", size = 7) +
  geom_sf_label(data = st_centroid(st_union(df.stat[df.stat$station.project_id == "maisD",])), 
          label = "maisD", size = 7) +
  annotation_scale(height = unit(0.4, "cm"), text_cex = 1, width_hint = 0.4, bar_cols = c("grey60", "white")) +
  coord_sf(xlim = c(dimA[1], dimA[3]), ylim = c(dimA[2], dimA[4]), expand = F) +
  theme_void(base_size = 15) +
  facet_wrap(~"study area") +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey60", color = "grey60"),
        strip.text = element_text(size = 30, colour="white"))

g <- gA | (gC/gD)  # + plot_layout(guides = 'collect')

ggsave("./plots/plotSite.pdf", plot = g, width = 30, height = 18, device = "pdf", units = "cm")
