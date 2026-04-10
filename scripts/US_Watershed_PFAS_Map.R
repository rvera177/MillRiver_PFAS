# Plotting PFAS data across the united states.
library(readr) #bring in the spatial coordinates with coresponding PFAS data
library(sf) #plotting the spatial objects
library(ggplot2)
#install.packages("rnaturalearth")
library(rnaturalearth) #for plotting country  outlines
#install.packages("rnaturalearthdata")
library(rnaturalearthdata)
library(dplyr) #need this for pipes

#I want to bring in the papers our group has been extracting data from

Caravan_PFAS = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Caravan_PFAS_2026.csv")
Camacho_2024 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Camacho_et_al_2024_Florida.csv")
Sims_2025 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Sims_et_al_2025_%20Western_United_States.csv")
NH_DES_2026 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/NewHampshire_DES_PFAS_Data_Dump.csv")
Ahrens_2023 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Ahrens_et_al_2023_Arctic.csv")
Breitmeyer_2023 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Breitmeyer_et_al_2023_Pennsylvania.csv")
Sharma_2016 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Sharma_et_al_2016_Ganges_River.csv")
Zhang_2016 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Zhang_et_al_2016_RI_NY.csv")
Scott_2009 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Scott_et_al_2009_Canada.csv")
Goodrow_2020 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Goodrow_et_al_2020_New_Jersey.csv")
Bai_Son_2021 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Bai_and_Son_2021_Renoe_LasVegas.csv")
Teymoorian_2021 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Teymoorian_2025_Montreal.csv")
Maine_DEP_2026 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/MaineDEP_2026_Datadump_cleaned.csv")
WQP_USA_2026 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/WQP_USA_Data_complete.csv")
Hayworth_2022 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Hayworth_et_al_2022_Alabama_cleaned.csv")
Dunn_2023 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Dunn_et_al_2023_RhodeIsland_complete.csv")
AustraliaMap_2026 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Australia_Government_PFAS_CHEM_MAP_Clean.csv")
Forster_2024 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Forster_et_al_2024_SouthCarolina_cleaned.csv")
Penland_2020 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Penland_2020_SC_NC_cleaned.csv")
Labad_2025 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Labad_et_al_2025_Georgia.csv")
Webb_2026 = read_csv("https://raw.githubusercontent.com/rvera177/GlobalPFAS/refs/heads/main/data/complete/Webb_et_al_2026_Savannah.csv")

All_PFAS = bind_rows(Camacho_2024, Sims_2025, NH_DES_2026, 
                     Breitmeyer_2023, Ahrens_2023, Sharma_2016, Caravan_PFAS,
                     Zhang_2016, Scott_2009, Goodrow_2020, Bai_Son_2021, 
                     Teymoorian_2021, Maine_DEP_2026, WQP_USA_2026, 
                     Hayworth_2022, Dunn_2023, AustraliaMap_2026, Forster_2024, 
                     Penland_2020, Labad_2025, Webb_2026)

All_PFAS_SW <- subset(All_PFAS, `Sample Type` == "Surface Water")

#this is the total number of sites
All_PFAS_SW %>%
  distinct(Latitude, Longitude) %>%
  nrow()

#ggplot of an sf object won't work if their are NA's in the lat and long
All_PFAS_sf <- All_PFAS_SW %>% #remove NA's and make numeric. if not numeric, turns into an NA
  mutate(across(-all_of(c("Article (Author et al YYYY)", "Sample Name", "Sample Type",
                          "Sample Date (MM/DD/YYY)", "Sample Time", "Analysis Method")), ~ as.numeric(.))) %>%
  filter(!is.na(`Latitude`), !is.na(`Longitude`))
#convert to SF object for plotting
All_PFAS_sf= st_as_sf(All_PFAS_sf,
                      coords = c("Longitude", "Latitude"),  # x = Long, y = Lat
                      crs = 4326) 


#everything above is my own code
#the following is leaflet code made with help from UMass Gen AI because I never made this before
#install.packages("leaflet")
#install.packages("htmlwidgets")

library(leaflet)
library(htmlwidgets)

# Build popup only from data columns
make_popup <- function(df) {
  if (inherits(df, "sf")) df <- sf::st_drop_geometry(df)
  
  df_chr <- as.data.frame(lapply(df, function(x) {
    ifelse(is.na(x), "", as.character(x))
  }), check.names = FALSE)
  
  apply(df_chr, 1, function(r) {
    paste0("<b>", names(df_chr), "</b>: ", r, collapse = "<br/>")
  })
}
All_PFAS_SW <- All_PFAS_SW %>% select(-any_of("popup"))
All_PFAS_SW$popup <- make_popup(All_PFAS_SW)

# Build the leaflet map
RV_Teja_map <- leaflet(options = leafletOptions(minZoom = 2, maxZoom = 18)) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(data = All_PFAS_SW %>% filter(!is.na(PFOA)),
                   lng = ~Longitude,   # ← added
                   lat = ~Latitude,    # ← added
                   color = "red", fillColor = "red",
                   radius = 5, stroke = FALSE, fillOpacity = 0.9,
                   popup = ~popup,
                   group = "PFOA") %>%
  addCircleMarkers(data = All_PFAS_SW %>% filter(!is.na(PFOS)),
                   lng = ~Longitude,   # ← added
                   lat = ~Latitude,    # ← added
                   color = "blue", fillColor = "blue",
                   radius = 3, stroke = FALSE, fillOpacity = 0.9,
                   popup = ~popup,
                   group = "PFOS") %>%
  addLayersControl(overlayGroups = c("PFOA", "PFOS"),
                   options = layersControlOptions(collapsed = FALSE)) %>%
  addLegend(position = "topright",
            colors = c("red", "blue"),
            labels = c("PFOA", "PFOS"),
            title = "Compounds")

RV_Teja_map


#okay so now i want to start experimenting with an ML model
#starting off with just the united states, so i'm going to bound
#my data to only points in the US
remotes::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
# Load U.S. state boundaries

USA_states <- ne_states(country = "united states of america", returnclass = "sf")

# Clipping data to whats inside of the usa
US_PFAS <- st_join(All_PFAS_sf, USA_states, join = st_intersects, left = FALSE)
#okay so i have 5515 data points in the united states
ggplot(data = USA_states) +
  geom_sf(fill = "gray95", color = "gray20") +
  geom_sf(data = US_PFAS %>% filter(!is.na(PFOS)), color = "red", size = 2, shape = 19, na.rm=TRUE) +
  ggtitle("United States Map of PFOA") +
  theme_classic()

# crop in lon/lat
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, xmax = -66, ymin = 22, ymax = 51), crs = st_crs(4326)))
USA_states_crop <- st_crop(USA_states, conus_bbox)
US_PFAS_crop <- st_crop(US_PFAS, conus_bbox)

ggplot() +
  geom_sf(data = USA_states_crop, fill = "white", color = "gray2") +
  geom_sf(data = US_PFAS_crop %>% filter(!is.na(PFOA)), color = "black",fill  = "red",  size = 2, stroke = 0.5,  shape = 21, na.rm = TRUE) +
  theme_classic()

#only looking at New england huc2
library(nhdplusTools)
library(dataRetrieval)

# Get all HUC2 boundaries directly
huc2 <- st_read("C:/Users/Marston User/Documents/Global Papers/HUC_Boundaries/Watershed_Boundary_Dataset_HUC_2s_Source_Resilience_Climate_GovDataset/HU02.shp")
#huc2 <- st_read("Watershed_Boundary_Dataset_HUC_2s_-6033951573301036398/HU02.shp")

huc2 %>% st_drop_geometry() %>% select(huc2, name)


# Five HUC regions at the moment have the most data. Ploting them here. 
# Huc = "01", "02", "03", "04", "15"
##NEW ENGLAND

# Filter HUC2 to New England only (huc2 == "01")
#ideally, it would be nice to have a map of the whole use, with different HUC2 boundaries different colours
ne_huc2 <- huc2 %>%
  filter(huc2 == "01")

# Match CRS to your PFAS sf object
ne_huc2 <- st_transform(ne_huc2, st_crs(All_PFAS_sf))

# Clip PFAS observations to New England
sf_use_s2(FALSE)
NE_PFAS <- st_intersection(All_PFAS_sf, ne_huc2)
sf_use_s2(TRUE)

# Quick check
cat("New England observations:", nrow(NE_PFAS), "\n")
cat("Unique sites:", NE_PFAS %>% distinct(geometry) %>% nrow(), "\n")

# Plot to verify
ggplot() +
  geom_sf(data = ne_huc2, fill = "gray95", color = "gray20") +
  geom_sf(data = NE_PFAS %>% filter(!is.na(PFOS)), 
          color = "black", fill = "red",
          size = 2, stroke = 0.5, shape = 21, na.rm = TRUE) +
  ggtitle("New England PFAS observations (HUC2 = 01)") +
  theme_classic()

#now for mid atlantic
midatlantic_huc2 <- huc2 %>%
  filter(huc2 == "02")

# Match CRS to your PFAS sf object
midatlantic_huc2 <- st_transform(midatlantic_huc2, st_crs(All_PFAS_sf))

# Clip PFAS observations to midatlantic
sf_use_s2(FALSE)
midatlantic_PFAS <- st_intersection(All_PFAS_sf, midatlantic_huc2)
sf_use_s2(TRUE)

cat("Mid Atlantic observations:", nrow(midatlantic_PFAS), "\n")
cat("Unique sites:", midatlantic_PFAS %>% distinct(geometry) %>% nrow(), "\n")

ggplot() +
  geom_sf(data = midatlantic_huc2, fill = "gray95", color = "gray20") +
  geom_sf(data = midatlantic_PFAS %>% filter(!is.na(PFOS)), 
          color = "black", fill = "red",
          size = 2, stroke = 0.5, shape = 21, na.rm = TRUE) +
  ggtitle("mid atlantic PFAS observations (HUC2 = 02)") +
  theme_classic()


##SOUTHEAST
#(huc2 == "03")
SE_huc3 <- huc2 %>%
  filter(huc2 == "03")

# Match CRS to your PFAS sf object
SE_huc3 <- st_transform(SE_huc3, st_crs(All_PFAS_sf))

# Clip PFAS observations to New England
sf_use_s2(FALSE)
SE_PFAS <- st_intersection(All_PFAS_sf, SE_huc3)
sf_use_s2(TRUE)

# Quick check
cat("Southeast observations:", nrow(SE_PFAS), "\n")
cat("Unique sites:", SE_PFAS %>% distinct(geometry) %>% nrow(), "\n")

# Plot to verify
ggplot() +
  geom_sf(data = SE_huc3, fill = "gray95", color = "gray20") +
  geom_sf(data = SE_PFAS %>% filter(!is.na(PFOS)), 
          color = "black", fill = "red",
          size = 2, stroke = 0.5, shape = 21, na.rm = TRUE) +
  ggtitle("Southeast PFAS observations (HUC2 = 03)") +
  theme_classic()

#Fireeee 
#A little concerned about the lack of data in the upper right hand corner but
#I think we just have less data there... we will see




#NORTHWEST
#(huc2 == "17")
#going to start naming these huc(watershed#) because I think I will get more 
#confused if I just do them numerically

LowerColorado_huc15 <- huc2 %>%
  filter(huc2 == "15")

# Match CRS to your PFAS sf object
LowerColorado_huc15 <- st_transform(LowerColorado_huc15, st_crs(All_PFAS_sf))

# Clip PFAS observations to New England
sf_use_s2(FALSE)
LowerColorado_PFAS <- st_intersection(All_PFAS_sf, LowerColorado_huc15)
sf_use_s2(TRUE)

# Quick check
cat("LowerColorado observations:", nrow(LowerColorado_PFAS), "\n")
cat("Unique sites:", LowerColorado_PFAS %>% distinct(geometry) %>% nrow(), "\n")

# Plot to verify
ggplot() +
  geom_sf(data = LowerColorado_huc15, fill = "gray95", color = "gray20") +
  geom_sf(data = LowerColorado_PFAS %>% filter(!is.na(PFOS)), 
          color = "black", fill = "red",
          size = 2, stroke = 0.5, shape = 21, na.rm = TRUE) +
  ggtitle("LowerColorado PFAS observations (HUC2 = 15)") +
  theme_classic()


#I am wondering if there is an easier way to do this where we
#can highlight the watersheds without making them their own individual plots
#just kidding my r keeps crashing when I try to run this

#install.packages("rmapshaper")
library(rmapshaper)
#CONUS = continuous united states map 

##Map of all watersheds colored in 
#------------------
# Crop watersheds to CONUS using approximate long/lat just to remove Alaska/Hawaii
huc2_simple <- ms_simplify(huc2, keep = 0.05, keep_shapes = TRUE)
huc2_simple <- st_transform(huc2_simple, 4326)

# Automatically fix self-intersections / duplicate vertices
huc2_simple <- st_make_valid(huc2_simple)

# Optional: force geometry type to MULTIPOLYGON to avoid other errors
huc2_simple <- st_cast(huc2_simple, "MULTIPOLYGON")

conus <- st_crop(huc2_simple,
                 st_bbox(c(xmin = -145, xmax = -66, ymin = 18, ymax = 50),
                         crs = st_crs(4326)))



# Filter PFAS points to the cropped watersheds
pfas_conus <- st_filter(All_PFAS_sf, conus)



# Plot
ggplot() +
  geom_sf(data = conus, aes(fill = factor(huc2)), color = "gray40", alpha = 0.5) +
  geom_sf(data = pfas_conus %>% filter(!is.na(PFOA)),
          color = "black", size = 0.4, alpha = 0.7) +
  scale_fill_viridis_d(name = "HUC2 Region") +
  theme_classic() +
  ggtitle("PFOA Observations Across US Watersheds") +
  coord_sf(expand = FALSE)  # <-- removes padding but keeps automatic limits

#I feel like the northeast looks kind of weird and the map looks distorted but I am not sure why
#I am guessing it's just the map projection? 
#Will come back to try and fix this later






#Now I want to try making it so that only certain regions are highlighted 

# Define which HUC2s to highlight (just randomly choosing these as key watersheds because I think we have the most data here)
highlight_hucs <- c("01", "02", "03", "04", "15")

# Separate highlighted and background watersheds
highlight_watersheds <- conus %>% dplyr::filter(huc2 %in% highlight_hucs)
background_watersheds <- conus %>% dplyr::filter(!huc2 %in% highlight_hucs)

# Plot
ggplot() +
  # Background watersheds in light gray
  geom_sf(data = background_watersheds, fill = "tan", color = "gray40", alpha = 0.5) +
  # Highlighted watersheds with distinct fill
  geom_sf(data = highlight_watersheds, aes(fill = factor(huc2)), color = "gray40", alpha = 0.7) +
  # PFAS points
  geom_sf(data = pfas_conus %>% filter(!is.na(PFOA)),
          color = "black", size = 0.2, alpha = 0.7) +
  # Fill scale for highlighted HUC2s
  scale_fill_viridis_d(name = "Highlighted HUC2 Region") +
  theme_classic() +
  ggtitle("PFOA Observations Across Selected US Watersheds") +
  coord_sf(expand = FALSE)


#The distortion seen in this map looks the same as the poster map Raul sent so I think it's ok









#UNFINISHED CODE BELOW, amalgamated from 3 or so chatgpt prompts but my R studio keeps crashing :(

#Okay now I want to add the 4s lines 
huc4 <- st_read("Watershed_Boundary_Dataset_HUC_4s_-978931873943038265/HU04.shp")

huc4 %>% st_drop_geometry() %>% select(huc4, name)

# Reproject all layers to a common CRS (Albers Equal Area for CONUS)
target_crs <- 5070
conus_proj <- st_transform(conus, crs = target_crs)
pfas_conus_proj <- st_transform(pfas_conus, crs = target_crs)
huc4_proj <- st_transform(huc4, crs = target_crs)

# Plot
ggplot() +
  geom_sf(data = conus_proj, aes(fill = factor(huc2)), color = "gray40", alpha = 0.5) +
  geom_sf(data = huc4_proj, color = "white", fill = NA, size = 0.3) +  # HUC4 boundaries on top
  geom_sf(data = pfas_conus_proj %>% filter(!is.na(PFOA)),
          color = "black", size = 0.4, alpha = 0.7) +
  scale_fill_viridis_d(name = "HUC2 Region") +
  theme_classic() +
  ggtitle("PFOA Observations Across US Watersheds") +
  coord_sf(expand = FALSE)









library(sf)
library(ggplot2)
library(dplyr)
library(viridis)

# --- 1. Read HUC4 shapefile ---

#these huc boundaries were downloaded from;
#https://resilience.climate.gov/datasets/esri::watershed-boundary-dataset-huc-4s/explore?location=28.414693%2C0.313494%2C0
huc4 <- st_read("C:/Users/Marston User/Documents/Global Papers/HUC_Boundaries/Watershed_Boundary_Dataset_HUC_4s_-978931873943038265/HU04.shp")
#I also want to add major lakes
# Download lakes (scale = 'medium' or 'large' for more detail)
lakes <- ne_download(scale = "large", type = "lakes", category = "physical", returnclass = "sf")

# Project to match your map CRS
lakes_proj <- st_transform(lakes, crs = target_crs)
# Clip lakes to only where they intersect with your CONUS HUC layer
lakes_proj <- st_intersection(lakes_proj, st_union(conus_proj))
# --- 2. Crop HUC2 to CONUS extent ---
conus <- st_crop(
  huc2_simple,
  st_bbox(c(xmin = -125, xmax = -66, ymin = 22, ymax = 60),
          crs = st_crs(4326))
)

# --- 3. Filter PFAS points to cropped watersheds ---
pfas_conus <- st_filter(All_PFAS_sf, conus)

# --- 4. Define target CRS (Albers Equal Area for CONUS) ---
target_crs <- 5070
conus_proj <- st_transform(conus, crs = target_crs)
pfas_proj <- st_transform(pfas_conus, crs = target_crs)
huc4_proj <- st_transform(huc4, crs = target_crs)

# Base plot: CONUS + HUC4 +lakes 
ggplot() +
  geom_sf(data = conus_proj, aes(fill = factor(huc2)), color = "gray40", alpha = 0.5) +
  geom_sf(data = huc4_proj, fill = NA, color = "white", size = 0.3) +
  geom_sf(data = lakes_proj, fill = "steelblue", color = "steelblue", size = 0.3) +
  scale_fill_viridis_d(name = "HUC2 Region") +
  geom_sf(data = pfas_proj, color = "black", size = 0.4, alpha = 0.7) +
  theme_classic() +
  ggtitle("PFOA Observations Across US Watersheds") +
  coord_sf(
    xlim = c(-125, -66),
    ylim = c(24, 55),
    crs = 4326,# interpret xlim/ylim as lon/lat
    expand = FALSE
  )


