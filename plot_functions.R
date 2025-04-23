#### plot functions ####

## add stations

plot.station <- function(g, site = NULL, type = NULL, size = NULL, stroke = NULL, alpha = 0.7,
                  data = shp.stat, site.name = "station.project_id", type.name = "type") {
  
  require(tidyverse)
  require(sf)
  require(ggnewscale)
  
  if(is.null(site)) stop("Define a site (e.g. 'maisC'")
  if(is.null(size)) size <- 3
  if(is.null(stroke)) stroke <- 1.2
  
  data$site.name <- unlist(st_drop_geometry(data[, site.name]))
  data$type.name <- unlist(st_drop_geometry(data[, type.name]))
  
  g <- g +
    new_scale_color() +
    new_scale_fill() +
    
    geom_sf(data = data[data$site.name %in% site & data$type %in% type,], 
            aes(pch = type.name, fill = type.name), color = "black", stroke = stroke, size = size, alpha = alpha) +
    scale_shape_manual("station type", values = c("direct" = 23, "omni" = 21),
                       guide = guide_legend(order = 2)) +
    scale_fill_manual("station type", values = c("direct" = "white", "omni" = "grey80"),
                      guide = guide_legend(order = 2))
  
  return(g)
  
}


plot.raster <- function(g, leg = " ", lim = c(0, 200)){
  
  require(tidyverse)
  require(sf)
  require(ggnewscale)
  
  g <- g +
    scale_fill_viridis_c(leg, direction = -1, end = 0.9, option = "inferno", limits = lim, 
                         oob = scales::oob_squish, na.value = c(NA, "black")) +
    scale_color_viridis_c(leg, direction = -1, end = 0.9, option = "inferno", limits = lim, 
                          oob = scales::oob_squish, na.value = c(NA, "black")) +
    ggtitle(paste0(m, " - r ", d)) +
    facet_wrap(~Individual) +
    annotation_scale(width_hint = 0.2, bar_cols = c("grey60", "white"), 
                     location = "tl", height = unit(0.15, "cm"))
  return(g)
}

plot.theme <- function(x) {
  
  require(ggplot2)
  
  g <- g  +
    coord_sf(xlim = c(dim[1], dim[3]), ylim = c(dim[2], dim[4])) +
    theme_light() +
    theme(
      ## rotates x axis labels by 45° 
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}
