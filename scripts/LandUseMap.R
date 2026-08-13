
# you don't need all these packages. but i don't remeber which ones are for what
#download all of these packages and then load
library(remotes)
library(StreamCatTools)
library(dplyr) #need dplyr loaded in order to pipe. Don't forget 
library(readr) 
library(sf)
library(nhdplusTools) #pulls in flowlines and associated COMIDS
library(writexl)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(viridis)
library(stringr)
library(nngeo)




start_point <- st_sfc(st_point(c(-72.5837008452757, 42.38344445112979)), crs = 4269)
start_comid <- discover_nhdplus_id(start_point)
start_comid

#now that i Have the starting comid, I collect all flowlines upstream
flowline <- navigate_nldi(
  list(featureSource = "comid", featureID = start_comid),
  mode = "upstreamTributaries",
  distance_km = 100) #goes 100km upstream from starting point
#this will be increased for larger watershed in the future

#save the flowlines 
subset_file <- "LWMR_nhd_subset.gpkg"

#this can take a while. 
#Occasionally the NHDplus website is down, so it won't work at a given moment.
#just give it an hour and it will be back up, usually...
subset <- subset_nhdplus(
  comids = as.integer(flowline$UT$nhdplus_comid),
  output_file = subset_file,
  nhdplus_data = "download",
  flowline_only = FALSE,
  return_data = TRUE,
  overwrite = TRUE
)
#make objects for each thing and plot
flowline  <- sf::read_sf(subset_file, "NHDFlowline_Network")
catchment <- sf::read_sf(subset_file, "CatchmentSP")
waterbody <- sf::read_sf(subset_file, "NHDWaterbody")
plot(st_geometry(catchment), col = "black")
plot(st_geometry(flowline), add = TRUE, col = "lightblue")
plot(start_point, add = TRUE, col = "red", pch = 19)
plot(st_geometry(waterbody), add = TRUE, col = "lightblue")
#make sure everything is in the same projection. It should be already.

#removed 1 problematic flowline at atkins reservoir
flowline = filter(flowline, !(objectid %in% c(47442)))


#replotting with the removed problematic flowline
plot(st_geometry(flowline), col = "blue")
plot(start_point, add = TRUE, col = "red", pch = 19)
#flowlines and prediction points are in!!

catchment_valid <- st_make_valid(catchment)
catchment <- st_union(catchment_valid) #combining my subcatchments



plot(st_geometry(flowline), col = "blue")
#plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)
plot(st_geometry(catchment), add = TRUE, border = "black", lwd = 4)

#renaming obs_clip into obs with different CRS (in meters)
flowline <- st_transform(flowline, crs =5070)
catchment <- st_transform(catchment_union, crs =5070)
st_area(catchment)

# Now plotting everything together
plot(st_geometry(flowline), col = "blue")
plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)


#  using NLCD raster for making a land cover map
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(tidyterra)  
library(FedData)

# Downloading NLCD for your watershed extent
nlcd <- get_nlcd(
  template = catchment,  # uses your catchment boundary as the extent
  label = "mill_river",
  year = 2021,           # most recent available
  dataset = "landcover"
)
# Check what land cover classes exist in your downloaded NLCD
unique(values(nlcd))
# Plot NLCD raster with cleaner rendering
# First convert nlcd values to factor
nlcd_factor <- as.factor(nlcd)

ggplot() +
  geom_spatraster(data = nlcd_factor, aes(fill = Class)) +
  scale_fill_manual(
    values = c(
      "11" = "#5475A8",
      "21" = "#DDC9C9",
      "22" = "#D89382",
      "23" = "#ED0000",
      "24" = "#AA0000",
      "31" = "#B2B2B2",
      "41" = "#38814E",
      "42" = "#1B6330",
      "43" = "#5CA253",
      "52" = "#CDB577",
      "71" = "#E8D44D",
      "81" = "#FBF65D",
      "82" = "#CA9146",
      "90" = "#7BA7BC",
      "95" = "#BAD8EA"
    ),
    labels = c(
      "11" = "Open water",
      "21" = "Developed open",
      "22" = "Developed low",
      "23" = "Developed medium",
      "24" = "Developed high",
      "31" = "Barren land",
      "41" = "Deciduous forest",
      "42" = "Evergreen forest",
      "43" = "Mixed forest",
      "52" = "Shrub/scrub",
      "71" = "Grassland",
      "81" = "Pasture",
      "82" = "Cultivated crops",
      "90" = "Woody wetlands",
      "95" = "Emergent wetlands"
    ),
    na.translate = FALSE,
    name = "Land cover"
  ) +
  geom_sf(data = flowline, color = "#5DADE2", linewidth = 0.4) +
  geom_sf(data = catchment, fill = NA, color = "black", linewidth = 0.8) +
  coord_sf(crs = st_crs(catchment)) +
  theme_void() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9))

# Reclassify raster first
nlcd_reclass <- classify(nlcd, 
                         rbind(
                           c(11, 1),   # open water
                           c(21, 2),   # developed
                           c(22, 2),   # developed
                           c(23, 2),   # developed
                           c(24, 2),   # developed
                           c(31, 3),   # barren
                           c(41, 4),   # forest
                           c(42, 4),   # forest
                           c(43, 4),   # forest
                           c(52, 5),   # shrub/scrub -> group with forest or barren
                           c(71, 5),   # grassland
                           c(81, 6),   # agricultural
                           c(82, 6),   # agricultural
                           c(90, 7),   # woody wetlands
                           c(95, 7)    # emergent wetlands
                         )
)

nlcd_reclass <- as.factor(nlcd_reclass)

