

#This script is part of the global PFAS modeling project
#Created February 2, 2026 by RV
#Last edited: February 2, 2026 by RV
#Edits will be uploaded to Github for easy access

#working with global caravan data. 
#this data comes from caravan-qual lite. Released December 11, 2025
#link to lite dataset: https://zenodo.org/records/17787066 

#"heavy" dataset can be found here: https://github.com/SustainableWaterSystems/Caravan-Qual 
library(readr) #bring in the spatial coordinates with coresponding PFAS data
library(sf) #plotting the spatial objects
library(ggplot2)
library(rnaturalearth) #for plotting country  outlines
library(rnaturalearthdata)

#pull in global caravan data from my github. 
PFOA_Raw <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Caravan_PFOA.csv")
PFOS_Raw = read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Caravan_PFOS.csv")
site_info = read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Caravan_wqms_site_info.csv")

#combine PFOA and PFOS by column names
PFAS_Raw = rbind(PFOA_Raw,PFOS_Raw)

#the raw datasets are in micrograms per liter. Convert to nanograms per liter
  PFAS <- PFAS_Raw %>% 
  mutate(obs = obs* 1000) #converts ug/L column to ng/L by simple multiplying

#combine site coordinate information with PFAS concentrations
PFAS <- left_join(PFAS, site_info, by = "wqms_id")


#convert to SF object for plotting
PFAS= st_as_sf(PFAS,
                  coords = c("wqms_lon", "wqms_lat"),  # x = Long, y = Lat
                  crs = 4326) 

#plot them on a world map
world <- ne_countries(scale = "medium", returnclass = "sf")
#update the coordinate system to match the map
PFAS = st_transform(PFAS, st_crs(world))

ggplot(data = world) +
  geom_sf(fill = "gray95", color = "gray20") +
  geom_sf(data = PFOA_sf, color = "red", size = 2, shape = 19) +
  geom_sf(data = PFOS_sf, color = "blue", size = 2, shape = 19) +
  ggtitle("Global Map of PFOA and PFAS") +
  theme_classic()

#everything above is personal code
#following is leaflet code from umass gen AI
library(sf)
library(leaflet)
library(rnaturalearth)
library(htmlwidgets)   # for saving
#initial world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# single coordinate system
world  <- st_transform(world, 4326)
PFOA_sf <- st_transform(PFOA_sf, 4326)
PFOS_sf <- st_transform(PFOS_sf, 4326)

# 3. optional: build HTML popup strings from attributes (safe even if different columns)
make_popup <- function(sf_obj) {
  df <- sf::st_drop_geometry(sf_obj)
  if (ncol(df) == 0) return(rep("", nrow(df)))
  apply(df, 1, function(r) {
    paste0("<b>", names(df), "</b>: ", ifelse(is.na(r), "", r), collapse = "<br/>")
  })
}
PFOA_sf$popup <- make_popup(PFOA_sf)
PFOS_sf$popup <- make_popup(PFOS_sf)

#build the leaflet map
Caravan_PFOA_PFAS_map <- leaflet(options = leafletOptions(minZoom = 2, maxZoom = 18)) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  # world polygon
  addPolygons(data = world,
              color = "#444444", weight = 1,
              fillColor = "lightgrey", fillOpacity = 0.5,
              group = "World",
              popup = ~name) %>%           # change to whatever world attribute you prefer
  # PFOA points
  addCircleMarkers(data = PFOA_sf,
                   color = "red", fillColor = "red",
                   radius = 5, stroke = FALSE, fillOpacity = 0.9,
                   popup = ~popup,
                   group = "PFOA",
                   clusterOptions = markerClusterOptions()) %>%
  # PFOS points
  addCircleMarkers(data = PFOS_sf,
                   color = "blue", fillColor = "blue",
                   radius = 5, stroke = FALSE, fillOpacity = 0.9,
                   popup = ~popup,
                   group = "PFOS",
                   clusterOptions = markerClusterOptions()) %>%
  addLayersControl(overlayGroups = c("World", "PFOA", "PFOS"),
                   options = layersControlOptions(collapsed = FALSE)) %>%
  addLegend(position = "topright",
            colors = c("red", "blue"),
            labels = c("PFOA", "PFOS"),
            title = "Compounds")

# show map
Caravan_PFOA_PFAS_map

# (optional) save to an HTML file
# saveWidget(Caravan_PFOA_PFAS_map, "Caravan_PFOA_PFAS_map.html", selfcontained = TRUE)
