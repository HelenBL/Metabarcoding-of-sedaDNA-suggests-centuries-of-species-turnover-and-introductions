### Mapping using ggOceanmaps
######## Include Topography 

library(sp)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(sf)
library(tmap)
library(geodata)
library(ggnewscale)
library(maps)
library(cowplot)
library(ggspatial)
library(mapdata)
library(ggOceanMaps)
library(maptiles)
install.packages("marmap")
library(marmap)
library(leaflet)

# Download bathimetry of NOAA from England
bathy <- getNOAA.bathy(lon1 = -10.5, lon2 = 2, lat1 = 49, lat2 = 61, resolution = 0.1)
bathy_df <- fortify.bathy(bathy)

# Separate Ocean and Earth
bathy_r <- bathy_df
bathy_r$z[bathy_r$z >= 0] <- NA

relief_r <- bathy_df
relief_r$z[relief_r$z < 0] <- NA

# Convert to dataframes
bathy_df <- as.data.frame(bathy_r, xy = TRUE)
colnames(bathy_df) <- c("x", "y", "depth")
bathy_df <- na.omit(bathy_df)

relief_df <- as.data.frame(relief_r, xy = TRUE)
colnames(relief_df) <- c("x", "y", "elevation")
relief_df <- na.omit(relief_df)

bathy_palette <- colorRampPalette(c("#000080", "#0000CD", "#1E90FF", "#87CEFA"))
relief_palette <- colorRampPalette(c("#424242", "#6E8B3D", "#FFFFFF"))

bathy_df$color <- bathy_palette(100)[cut(bathy_df$depth, breaks = 100, labels = FALSE)]
relief_df$color <- relief_palette(100)[cut(relief_df$elevation, breaks = 100, labels = FALSE)]

ggplot() +
  # Bathimetry
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = I(color))) +
  
  geom_raster(data = relief_df, aes(x = x, y = y, fill = I(color)), alpha = 0.8) +
  
  coord_quickmap() +
  theme_minimal()



uk_map <- map_data("world", region = "UK")

map_main <- ggplot() +
  geom_raster(data = bathy_df, aes(x = x, y = y, fill = I(color))) +
  geom_raster(data = relief_df, aes(x = x, y = y, fill = I(color)), alpha = 0.8) +
 # geom_path(data = uk_map, aes(x = long, y = lat, group = group),
 #           color = "black", size = 0.5) +
  geom_rect(aes(xmin = -1.7, xmax = -1, ymin = 50.5, ymax = 51.1), color = "red", fill = NA, size = 1) +
  coord_sf(xlim = c(-10.5, 2), ylim = c(49, 61), crs = 4326)+
  theme_minimal() +
  annotation_scale(
    location = "bl",       # "bl" = bottom left
    width_hint = 0.2       
  )+
  annotation_north_arrow(
    location = "tl", which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1, "cm"), width = unit(1, "cm"),
    pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm")
  )
print(map_main)
europe_map <- map_data("world")

map_europe <- ggplot() +
  geom_polygon(data = europe_map, aes(x = long, y = lat, group = group), fill = "gray90", color = "black") +
  geom_rect(aes(xmin = -10.5, xmax = 2, ymin = 49, ymax = 61), color = "red", fill = NA, size = 1) +
  coord_quickmap(xlim = c(-25, 30), ylim = c(35, 70)) +
  theme_void()

print(map_europe)


leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addMarkers(lng = -1.4, lat = 50.9, popup = "Southampton")


# Sampling Sites
muestreo <- data.frame(
  site = c("Hythe", "Bunny Meadows"),
  lon = c(-1.4, -1.3),
  lat = c(50.9, 50.95),
  tipo = c("Hythe", "Bunny Meadows")  
)

muestreo_sf <- st_as_sf(muestreo, coords = c("lon", "lat"), crs = 4326)

# Interactive map
leaflet(muestreo_sf) %>%
  addProviderTiles("CartoDB.Positron") %>%  # Fondo claro y minimalista
  addCircleMarkers(
    radius = 8,
    color = ~ifelse(tipo == "Hythe", "#3FB8AF", "#F4D35E"),
    stroke = TRUE,
    fillOpacity = 0.9,
    label = ~tipo
  ) %>%
  addRectangles(
    lng1 = -1.7, lat1 = 50.5, lng2 = -1, lat2 = 51.1,
    color = "red", fillOpacity = 0, weight = 2
  )


bbox_soton <- st_as_sfc(st_bbox(c(
  xmin = -1.7, xmax = -1,
  ymin = 50.5, ymax = 51.1
), crs = 4326))

# Download CartoDB.Positron
tiles <- get_tiles(bbox_soton, provider = "CartoDB.Positron", crop = TRUE, zoom = 13)

# Convert raster to a dataframe for ggplot
tiles_df <- as.data.frame(tiles, xy = TRUE)
names(tiles_df)[3:5] <- c("R", "G", "B")

muestreo_sf <- st_as_sf(muestreo, coords = c("lon", "lat"), crs = 4326)

rectangulo <- st_as_sfc(st_bbox(c(
  xmin = -1.7, xmax = -1,
  ymin = 50.5, ymax = 51.1
), crs = 4326))



map_soton_svg <- ggplot() +
  geom_raster(data = tiles_df, aes(x = x, y = y, fill = rgb(R, G, B, maxColorValue = 255))) +
  scale_fill_identity() +
  
  # Sampling Point
  geom_sf(data = muestreo_sf, aes(color = tipo), size = 4) +
  scale_color_manual(values = c("Hythe" = "#3FB8AF", "Bunny Meadows" = "#F4D35E")) +
  
  # Red square
  geom_sf(data = rectangulo, fill = NA, color = "red", size = 1) +
  
  coord_sf(xlim = c(-1.7, -1), ylim = c(50.5, 51.1), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(color = "Sampling Site")+
  annotation_north_arrow(
    location = "tl",            # Top-left
    which_north = "true",       # Norte geográfico
    style = north_arrow_fancy_orienteering,
    height = unit(1.2, "cm"), width = unit(1.2, "cm"),
    pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm")
  ) +
  
  annotation_scale(
    location = "bl",            # Bottom-left
    width_hint = 0.2,           
    bar_cols = c("black", "white"),
    text_cex = 0.8,
    line_width = 1
  )


print(map_soton_svg)
