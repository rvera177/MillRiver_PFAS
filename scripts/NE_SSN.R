
#New England SSN model for PFAS compounds
#created in January 2026 - RV
#updated on 1/21/26 - RV


library(remotes)
library(StreamCatTools)
library(dplyr)
library(readr) #bring in the spatial coordinates with coresponding PFAS data
library(sf)
library(nhdplusTools) #pulling in flowlines and associated COMIDS
library(SSN2)
library(SSNbler)
library(writexl)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(viridis)
library(stringr)

setwd("C:/Users/Marston User/Documents/LWMR Isoscape")
#you can put whatever folder makes sense for you. 
#the location doesn't matter, it's just where final plots will be added to.
#you don't need to download any data before hand. 

PFAS_Spatial_Oct_2025 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Spatial_1_PFASResults.csv")
S1 <- PFAS_Spatial_Oct_2025 #need dplyr loaded in order to pipe. Don't forget lol

S1 <- S1 %>%#remove NA's for coordinates
  filter(!is.na(Lat), !is.na(Long))

clean_names <- function(x) {
  x |>
    gsub(" ", "_", x = _) |>
    gsub("[-:]", ".", x = _) |>
    make.names() |>
    gsub("^X(\\d)", "\\1", x = _)
}

colnames(S1) <- clean_names(colnames(S1))

#Naming family groups that I can add together
pfas_groups <- list(
  PFCA = c("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA",
           "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTrDA", "PFTeDA"),
  PFSA = c("PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS",
           "PFNS", "PFDS", "PFDoS"),
  FTSA = c("4.2FTS", "6.2FTS", "8.2FTS"),
  PFOSA = c("PFOSA", "NMeFOSA", "NEtFOSA"),
  FOSAA = c("NMeFOSAA", "NEtFOSAA"),
  PFECA = c("HFPO.DA", "ADONA", "PFMPA", "PFMBA", "NFDHA"),
  PFESA = c("9Cl.PF3ONS", "11Cl.PF3OUdS", "PFEESA"),
  FTCA = c("3.3_FTCA", "5.3_FTCA", "7.3_FTCA"))

#create new columns in dataset with Sum of corresponding species in the family
# Create new columns in dataset with sum of corresponding species in the family
S1 <- S1 %>%
  bind_cols(
    lapply(names(pfas_groups), function(fam_name) {
      # Get the column names in S1 for this family
      cols <- paste0(pfas_groups[[fam_name]], "_Results")
      cols_existing <- intersect(cols, names(S1))  # only existing columns
      
      tibble(
        !!fam_name := rowSums(S1[cols_existing], na.rm = TRUE)
      )
    })
  )

S1 <- S1 %>%
  mutate(OBSPRED_ID = row_number())


#snapping sites to nearest flowline so it gets the nearest comid
sites <- st_as_sf(S1,
                  coords = c("Long", "Lat"),
                  crs = 4326)

#the following chuck is for snapping prediction sites to a flowline.
bb <- st_bbox(sites) #bb=bounding box around each Observation sites
bb_poly <- st_as_sfc(st_bbox(sites))  # convert bbox to sfc
bb_sf <- st_sf(geometry = bb_poly)# convert to sf object
flines <- get_nhdplus(AOI = bb_sf, realization = "flowline") #Download flowlines each bounding box
idx <- get_flowline_index(flines, sites, max_matches = 1) #Snap points to nearest flowline & get COMID
S1$COMID <- as.numeric(idx$COMID) #comid ID for each point! Don't need to use idx anymore.

###------Exploring the existing data --------------
# (use StreamCatTools documentation for metric names) 
#https://www.epa.gov/national-aquatic-resource-surveys/streamcat-metrics-and-definitions
streamcat_data_cat <- sc_get_data( comid = S1$COMID, 
                                   metric = c("conn", "npdesdens", "pctimp2019",
                                              "pcturbhi2019", "pcturblo2019",
                                              "pcturbmd2019", "pcturbop2019",
                                              "huden2010", "rdcrs"),
                                   aoi = "Cat" ) #area of interst, whole watershed
streamcat_data_ws <- sc_get_data( comid = S1$COMID, 
                                  metric = c("conn", "npdesdens", "pctimp2019",
                                             "pcturbhi2019", "pcturblo2019",
                                             "pcturbmd2019", "pcturbop2019",
                                             "huden2010", "rdcrs"),
                                  aoi = "ws" ) #area of interst, whole watershed

#renaming since there was a name conversion using sc_get_data
streamcat_data_cat <- streamcat_data_cat %>% rename(COMID = comid)
S1_Cat <- left_join(S1, streamcat_data_cat, by ="COMID") 

#Make sure that your working model uses S1!
streamcat_data_ws <- streamcat_data_ws %>% rename(COMID = comid)
S1 <- left_join(S1, streamcat_data_ws, by ="COMID") 

cor_results_ws <- S1 %>%
  select(PFAS40, PFCA, PFSA, PFOA_Results, connws, npdesdensws, pctimp2019ws, pcturbhi2019ws, pcturblo2019ws, pcturbmd2019ws, pcturbop2019ws) %>% 
  cor(method = "spearman", use = "complete.obs")
cor_long_ws <- melt(cor_results_ws)

cor_results_Cat <- S1_Cat %>%
  select(PFAS40, PFCA, PFSA, PFOA_Results, conncat, npdesdenscat, pctimp2019cat, pcturbhi2019cat, pcturblo2019cat, pcturbmd2019cat, pcturbop2019cat) %>% 
  cor(method = "spearman", use = "complete.obs")
cor_long_cat <- melt(cor_results_Cat)

#plotting the correlation results. 

pheatmap(cor_results_Cat,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = TRUE,
         breaks = seq(-1, 1, length.out = 51),
         main = "Spearman Correlation Heatmap COMID scale")


pheatmap(cor_results_ws,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = TRUE,
         breaks = seq(-1, 1, length.out = 51),
         main = "Spearman Correlation Heatmap watershed scale")


cor_results_ws <- S1 %>%
  select(PFAS40, PFCA, PFSA, PFOA_Results, connws, npdesdensws, pctimp2019ws, 
         pcturbhi2019ws, pcturblo2019ws, pcturbmd2019ws, pcturbop2019ws) %>% 
  cor(method = "spearman", use = "complete.obs")

# Mask lower triangle + diagonal and then melt it
cor_results_ws[lower.tri(cor_results_ws, diag = TRUE)] <- NA
cor_long_ws <- melt(cor_results_ws, varnames = c("Variable1", "Variable2"), value.name = "Spearman_r", na.rm = TRUE)

# Reordering Variable2 so labels are aligned with color with tiles
cor_long_ws$Variable2 <- factor(cor_long_ws$Variable2, levels = rev(unique(cor_long_ws$Variable2)))

# Plot
ggplot(cor_long_ws, aes(x = Variable1, y = Variable2, fill = Spearman_r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Spearman_r, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1)
  )

#--------Setting up flowpaths for SSN-----------------------
#Using NHDPlus data set from R package
#changed to all of NEW England. All flowlines were downloaded from online.
# https://www.epa.gov/waterdata/nhdplus-northeast-data-vector-processing-unit-01
# Download Folder: NHDSnapshot_04.7Z 
# Load subfolder: NHDFlowline.shp

library(sf)
flowline <- read_sf("C:/Users/Marston User/Documents/NE PFAS Data/NHDPlusV21_NE_01_NHDSnapshot_04/NHDPlusNE/NHDPlus01/NHDSnapshot/Hydrography/NHDFlowline.shp")
catchment <- read_sf("C:/Users/Marston User/Documents/NE PFAS Data/NHDPlusV21_NE_01_NHDPlusCatchment_01/NHDPlusNE/NHDPlus01/NHDPlusCatchment/Catchment.shp")
sf::sf_use_s2(FALSE) #Sets the S2 spherical geometry to False. This forces sf to use planar geometry... aparently idk
catchment <- st_union(catchment) #this merges all catchments into one big one. 
catchment <- st_sf(catchment) #making sure the merged results is an sf object
plot(st_geometry(catchment), col = "green4")
plot(st_geometry(flowline), add=TRUE, col = "blue")

#removed problematic flowlines. These are mostly downstream divergence issues,
# caused by islands, ponds, and/or swamp area.
# flowline = filter(flowline, !(objectid %in% c(50681,46590,52094,49521,
#                                               52242, 52243,51998,47024,46976, 47722,43199,44620,44653, 43148
#                                               , 46282,47090,47133,43611,43116, 43115,44449, 45041,40981, 41421
#                                               ,43277, 43278, 44955,51168,51163,51359,51326,51231,52448,51342
#                                               ,50437,50692,50785, 51759,52145,52822,52760,52256,52809,52724,51333
#                                               ,50226,51609, 51881,51885,52110, 52108,52105,52454,52521,46664
#                                               ,47536,51784,51753,51874,51878,51877,51749,51748,52541,52540,52554
#                                               ,47442,47613,48856,49891,52177,52814,52562,47574,47424,50016,52090
#                                               ,47396,46738,51922,48768,52073,52072,52069,48128,49249,47897, 51909
#                                               ,47583,48952,52152,49702,48675, 46795,49320,49311,49950,46052,48168
#                                               ,46081,45878,48021,48020,48214,46083,47282,47283,47284,47290,46952
#                                               ,46527,48011,46085, 48000,45435,45436,43205,43208,44697,45733,44675
#                                               ,44694,45774,45423, 43162,44599,43072,42771,43150,46282,47089,44073
#                                               ,44582,44561,43757,44539,43247,43463,44061, 44537,44533,44521, 44520
#                                               ,44519,44518,44538,42710,41721,43576,44427, 41958,44262,45079,42993
#                                               ,43568, 43280,41278,41279,41204, 40933,42963,43517,43502,43501,40909
#                                               ,41263,41269,41270,41209, 41256,41262,42420,40545,40588, 41106
#                                               ,49728, 49729, 49730, 49731, 49732
# )))

#create NSI (National Stream Internet) prediction points at the center of each flowline
#need to drop the Measurement = M dimension in the XYZM geometry of flowline
flowline = st_zm(flowline, drop=TRUE, what="ZM") #drop Elevation (Z) and Measure (M, usually Time) dimension from flowline

nsi_PredPoints <- flowline %>%
  st_centroid() %>%            # centroid works for LINESTRING
  st_cast("POINT") %>%         # ensures geometry is POINT
  mutate(comid = flowline$comid)  # carry over COMID into the flowline attribute table

plot(st_geometry(catchment), col = "green4")
plot(st_geometry(flowline), add=TRUE, col = "blue")
plot(st_geometry(nsi_PredPoints), add = TRUE, col = "red", pch = 19, cex = 1)
#flowlines and prediction points are in!!

# Convert the data from S1 to an sf object using Lat/Long
S1_sf <- st_as_sf(S1, 
                  coords = c("Long", "Lat"),  # x = Long, y = Lat
                  crs = 4326)                # WGS84

# 2. Transform to same CRS as flowline
S1_sf <- st_transform(S1_sf, st_crs(flowline))

# check the geometry
head(S1_sf)
plot(st_geometry(S1_sf), add = TRUE)

# Save as geopackage file so i can work on it in Arcgis if i want. simplicity too
st_write(S1_sf, "obs.gpkg", delete_layer = TRUE)

#error is due to time and name of a column. 
#it's still was created, so don't worry about it 

obs <- st_read("obs.gpkg")  
obs <- st_transform(obs, st_crs(flowline))
#obs_clip <- obs[st_within(obs, catchment_union, sparse = FALSE), ]

plot(st_geometry(catchment), col = "green4")
plot(st_geometry(flowline), add=TRUE, col = "blue")
plot(st_geometry(nsi_PredPoints), add = TRUE, col = "red", pch = 19)
plot(st_geometry(obs), add = TRUE, col = "black", pch = 19)

#giving Prediction points obspred ID numbers
nsi_PredPoints <- nsi_PredPoints %>% 
  mutate(
    OBSPRED = row_number() + 100000)

# get StreamCat data for the predicion point COMIDs
streamcat_data <- sc_get_data(
  comid = nsi_PredPoints$COMID,
  metric = c("npdesdens", "pctimp2019", "pcturblo2019"),  # add any other metrics you want
  aoi = "ws")  # watershed scale area of interest

nsi_PredPoints <- nsi_PredPoints %>% rename(comid = COMID) #new name = old name
streamcat_data <- streamcat_data %>% rename(comid = COMID) #new name = old name
head(streamcat_data)
#add streamcat data to prediction points
nsi_PredPoints <- left_join(nsi_PredPoints, streamcat_data, by = "comid")
#shorthen prediction points variables down to only the ones i need.
nsi_PredPoints <- nsi_PredPoints %>%
  select(OBSPRED, comid, totdasqkm, npdesdensws, pctimp2019ws, geom)

#renaming obs_clip into obs, and putting everything in coordinates that are meters based
flowline <- st_transform(flowline, crs =5070)
catchment <- st_transform(catchment, crs =5070)
obs <- st_transform(obs, crs =5070)
pred <- st_transform(nsi_PredPoints_clipped, crs =5070)


# Now plotting everything together with new coordinate system
plot(st_geometry(catchment), col = "green4")
plot(st_geometry(flowline), add=TRUE, col = "blue")
plot(st_geometry(obs), add = TRUE, col = "black", pch = 19)
plot(st_geometry(pred), add = TRUE, col = "red", pch = 19)


temp_dir <- "C:/Users/Marston User/Documents/LWMR Isoscape"
#change this to somewhere on your computer that makes sense.
dir.create(temp_dir, showWarnings = FALSE)
library(sf)

ssn_path <- file.path(temp_dir, "NE_model.ssn")
#won't properly run if you have edges up on ArcGIS
flowlines_2 = lines_to_lsn(flowline, 
                           lsn_path = temp_dir,
                           overwrite = TRUE)

#if you get any topology errors, follow next step.
#called Node Correction
#Otherwise, skip and move on to SSN Assemble stage

#------------Node Error Correction-----------------
node_errors = st_read("node_errors.gpkg")  
plot(st_geometry(node_errors), add = TRUE, col = "black", pch = 19)

#found a few topological errors. Going to download the flowpath to correct complex confluences
gpkg_path <- file.path("C:/Users/Marston User/Documents/ArcGIS/Projects/NE_SSN/NE_SSN.gdb")
st_write(flowline,
         dsn = gpkg_path,
         layer = "flowline_corrected",
         driver = "GPKG",
         quiet = FALSE)
#fixed in ARCGIS, now pulling in the fixed file into R
flowpathfixed<- "C:/Users/Marston User/Documents/LWMR Isoscape/flowline_corrected.gpkg"
sf::st_layers(flowpathfixed)
layers_info <- sf::st_layers(flowpathfixed)
print(layers_info)
layer_names <- layers_info$name
print(layer_names)

#read the flowlines (use the layer name listed by st_layers; "flowline_corrected")
flowline <- sf::st_read(dsn =flowpathfixed, layer = "flowline_corrected", quiet = FALSE)

#---------Assemble the SSN-----------------
#Run again!!

flowlines_2 = lines_to_lsn(flowline, 
                           lsn_path = temp_dir,
                           overwrite = TRUE)

#For LWMR catchment, the furthest observation from flowline is 138 meters
#This is at the Leverett Pond. I took the water sample at the opposite end of the pond 
#i put snap tolerence to 150m
obs <- sites_to_lsn(
  sites = obs,
  edges = flowlines_2,
  lsn_path = temp_dir,
  file_name = "obs",
  snap_tolerance = 150,
  save_local = TRUE,
  overwrite = TRUE)
#For lwmk catchment, the furthest prediction from flowline is 199 meters
#i put snap tolerence to 350m.
#this takes a while because there are so many prediction points
# one pred for each flowline.
preds <- sites_to_lsn(
  sites = pred,
  edges = flowlines_2,
  save_local = TRUE,
  lsn_path = temp_dir,
  file_name = "nsi_PredPoints",
  snap_tolerance = 350,
  overwrite = TRUE)
#12305 out of 12447 snapped at 350m distance. This is crazy far so IDK if I want to push it further
#yes
edges <- updist_edges(
  edges = flowlines_2,
  save_local = TRUE,
  lsn_path = temp_dir,
  calc_length = TRUE)
#yes
site.list <- updist_sites(
  sites = list(
    obs = obs,
    preds = preds),
  edges = edges,
  length_col = "Length",
  save_local = TRUE,
  lsn_path = temp_dir)

names(site.list) ## View output site.list names

edges <- afv_edges(
  edges = edges,
  infl_col = "totdasqkm",
  segpi_col = "areaPI",
  afv_col = "afvArea",
  lsn_path = temp_dir)

site.list <- afv_sites(
  sites = site.list,
  edges = edges,
  afv_col = "afvArea",
  save_local = TRUE,
  lsn_path = temp_dir)

names(site.list$preds) ## View column names in pred1km
names(edges) ## Look at edges column names

ggplot() +
  geom_sf(data = edges, aes(color = upDist)) +
  geom_sf(data = site.list$obs, aes(color = upDist)) +
  coord_sf(datum = st_crs(obs)) +
  scale_color_viridis_c()
#everything looks preped for assembly
NE_ssn <- ssn_assemble(
  edges = edges,
  lsn_path = temp_dir,
  obs_sites = site.list$obs,
  preds_list = site.list[c("preds")],
  ssn_path = paste0(temp_dir, "/NE.ssn"),
  import = TRUE,
  check = TRUE,
  afv_col = "afvArea",
  overwrite = TRUE
)
#SSN assembled. Great!

#-----------SSN Prediction------------
library(SSN2)

#plotting base SSN
ggplot() +
  geom_sf(
    data = NE_ssn$edges,
    color = "medium blue",
    aes(linewidth = totdasqkm)) +
  scale_linewidth(range = c(0.5, 2.0)) +
  geom_sf(
    data = NE_ssn$preds$preds,
    size = 0.5,
    shape = 21,
    fill = "white",
    color = "dark grey") +
  geom_sf(
    data = NE_ssn$obs,
    size = 1,
    aes(color = PFAS40)) +
  coord_sf(datum = st_crs(obs)) +
  scale_color_viridis_c() +
  labs(color = "PFAS (ppt)", linewidth = "WS Area") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

#plotting SSN without the prediction points for visibility
ggplot() +
  geom_sf(
    data = NE_ssn$edges,
    color = "medium blue",
    aes(linewidth = totdasqkm)) +
  scale_linewidth(range = c(0.5, 2.5)) +
  geom_sf(
    data = NE_ssn$obs,
    size = 1,
    aes(color = PFAS40)) +
  coord_sf(datum = st_crs(obs)) +
  scale_color_viridis_c() +
  labs(color = "PFAS (ppt)", linewidth = "WS Area") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

## Generate hydrologic distance matrices
ssn_create_distmat(NE_ssn)

#making a copy to a temporary directory so I'm not editing original ssn
path <- system.file("temp_dir/NE.ssn", package = "SSN2")

NE_pred <- ssn_import(
  path = NE_ssn$path,
  predpts = c("preds"),
  overwrite = TRUE
)
summary(NE_pred)

names(NE_pred$preds)

ggplot() +
  geom_sf(data = NE_pred$edges) +
  geom_sf(data = NE_pred$preds$preds, pch = 17, color = "red", size=0.3) +
  geom_sf(data = NE_pred$obs, color = "blue", size = 2) +
  theme_bw()

#hydrologic distance matrices that preserve directionality, 
#which are required for statistical modeling
#why did i do this twice? First time was for CR_SSN. Second for CR_pred. 
#shouldn't they be replicates? Why does the second take longer?
ssn_create_distmat(
  ssn.object = NE_pred,
  predpts = c("preds"),
  among_predpts = TRUE,
  overwrite = TRUE)

ggplot() +
  geom_sf(data = NE_pred$edges) +
  geom_sf(data = NE_pred$obs, aes(color = PFAS40), size = 2) +
  scale_color_viridis_c(limits = c(0, 60), option = "H") +
  theme_bw()

Togregram <- Torgegram(
  formula = PFAS40 ~ npdesdensws + pctimp2019ws,
  ssn.object = NE_pred,
  type = c("flowcon", "flowuncon", "euclid")
)
plot(Togregram)

#name the columns in the S1 file you want models for!
# generate all column names
all_cols <- names(NE_ssn$obs)
#view(all_cols)
# find the positions of the start and end columns that you want to make models for
start <- match("PFAS40", all_cols)
end   <- match("FTCA", all_cols)
#might have to change this "end" if having issues
# subset the column names
pfas_cols <- all_cols[start:end]
pfas_cols
#now i'm looking at all PFAS compounds
#so fun

#blank list for modeled compounds
models <- list()
skipped_compounds <- c()  # keep track of skipped compounds

#this takes my computer around 40 seconds to run a spatial model for ALL compounds. Sick.
for (compound in pfas_cols) {
  response <- NE_ssn$obs[[compound]]  # extract the response values
  # Check variability: at least 2 unique, non-NA values
  if(length(unique(na.omit(response))) < 2) {
    message(paste("Skipping", compound, "- not enough variability"))
    skipped_compounds <- c(skipped_compounds, compound)
    next  # skip this compound
  }
  
  # Fit model for each individual compound. Takes a few seconds
  form <- as.formula(paste0(compound, " ~ npdesdensws + pctimp2019ws"))
  models[[compound]] <- ssn_lm(
    formula = form,
    ssn.object = CR_ssn,
    tailup_type = "exponential",
    euclid_type = "gaussian",
    additive = "afvArea")
}

#comparing multiple mods and the original

#glances(PFAS40_mod, PFCA_mod, PFSA_mod, PFOA_mod)
#tidy(PFCA_mod, conf.int = TRUE)
#glance(PFCA_mod)
#logLik is the log likelihood
#plot(PFCA_mod, which = 1)

#PFOA is the best! 
#need to use the lowest AIC and AICc,
#which in this case comes from ss_mod.

#predict the compound concentration at each edge aka stream segment

preds <- list()
for (compound in pfas_cols) {
  preds[[compound]] <- augment(
    models[[compound]],
    newdata = "preds",  # <-- use all prediction sites in CR_ssn
    pred.type = "preds"
  )
}

for (compound in pfas_cols) {
  df <- preds[[compound]] %>%
    as.data.frame() %>%
    select(-comid) %>%  # remove any existing comid
    left_join(
      NE_ssn$preds$preds %>% st_drop_geometry() %>% select(pid, comid),
      by = "pid"
    ) %>%
    select(comid, .fitted) %>%
    rename(pred = .fitted)
  
  NE_pred$edges <- NE_pred$edges %>%
    left_join(df, by = "comid") %>%
    rename(!!paste0(compound, "_pred") := pred)
}
#need to fix this. Download the predicted edges to edit in ArcGIS
ggsave(NE_pred$edges,file=.gpkg)
#okay apparently ggsave is only for standard images, not spatial data
#save spatial data which is an sf object using st_write
st_write(obj=NE_pred$edges, dsn="Pred_edges.gpkg", layer="edges", driver="GPKG")

#SSN prediction pfas concentrations for all streams. Now to plot it
#------------------Plotting the SSN----------------------
#save plots to the working director, in a new folder called out_dir
out_dir <- "NE_PFAS_maps" #folder where the pfas maps are going to go
if (!dir.exists(out_dir)) dir.create(out_dir)

library(ggplot2)
library(sf)
library(dplyr)
library(wesanderson)
library(rlang)
library(scales)
library(broom)

# fixed scale limits
#check your maximum PFAS concentration in the prediction area before setting a scale_max
#view(CR_pred$edges)
scale_min <- 0
scale_max <- 200

# ensure obs_sf is an sf and matches edges CRS (as before)
if (!inherits(NE_ssn$obs, "sf")) {
  obs_sf <- st_as_sf(NE_ssn$obs)
} else {
  obs_sf <- NE_ssn$obs
}
edges_crs <- st_crs(NE_pred$edges)
if (is.na(st_crs(obs_sf)) || st_crs(obs_sf) != edges_crs) {
  obs_sf <- st_transform(obs_sf, edges_crs)
}
if (is.na(st_crs(catchment)) || st_crs(catchment) != edges_crs) {
  catchment <- st_transform(catchment, edges_crs)
}

# --- basic overview map of the SSN set up ---
edges_overview <- NE_pred$edges %>% mutate(feature = "Predicted ng/L")
obs_overview   <- obs_sf %>% mutate(feature = "Observed ng/L")

p_overview <- ggplot() +
  geom_sf(data = catchment, fill = "grey95", color = "black", size = 0.6) +
  geom_sf(data = edges_overview, aes(color = feature), linewidth = 0.6, inherit.aes = FALSE) +
  geom_sf(data = obs_overview, aes(shape = feature, fill = feature), size = 3, stroke = 0.6, inherit.aes = FALSE) +
  scale_color_manual(name = "Map key",
                     values = c("Predicted ng/L" = "black"),
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(name = "Map key",
                     values = c("Observed ng/L" = 21),
                     guide = guide_legend(order = 1)) +
  scale_fill_manual(name = "Map key",
                    values = c("Observed ng/L" = "black"),
                    guide = guide_legend(order = 1)) +
  labs(title = "Sample sites and stream segments") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

ggsave(filename = file.path(out_dir, "overview_map.png"), plot = p_overview,
       width = 8, height = 6, dpi = 300)

#--- forloop plots for all compounds------------------

# make sure obs_sf is an sf and matches edges CRS
if (!inherits(NE_ssn$obs, "sf")) obs_sf <- st_as_sf(NE_ssn$obs) else obs_sf <- NE_ssn$obs
edges_crs <- st_crs(NE_pred$edges)
if (is.na(st_crs(obs_sf)) || st_crs(obs_sf) != edges_crs) obs_sf <- st_transform(obs_sf, edges_crs)
if (is.na(st_crs(catchment)) || st_crs(catchment) != edges_crs) catchment <- st_transform(catchment, edges_crs)

# precompute glance table for R^2 lookup. This is for later
glance_df <- map_df(models, glance, .id = "compound")

# palette (dark blue into wes colors)
pal <- colorRampPalette(c("#010A4A", "#012A4A", wes_palette("Zissou1", type = "continuous"), "#8B0000"))(256)
#sqrt(256) = 16. changed to sqrt(400) with no noticible change
for (compound in pfas_cols) {
  if (!compound %in% names(models)) next
  
  pred_col <- paste0(compound, "_pred")
  obs_col  <- compound
  
  # get R^2 or pseudo-R^2
  r2_val <- glance_df %>%
    filter(compound == !!compound) %>%
    { if (nrow(.) == 0) NA_real_ else
      if ("pseudo.r.squared" %in% names(.)) .$pseudo.r.squared else
        if ("r.squared" %in% names(.)) .$r.squared else NA_real_ }
  r2_label <- if (is.na(r2_val)) "R² = NA" else paste0("R\u00B2 = ", format(r2_val, digits = 3))
  
  p <- ggplot() +
    # predicted stream segments (continuous color)
    geom_sf(data = NE_pred$edges,
            aes(color = !!sym(pred_col), linewidth = !!sym("totdasqkm")),
            show.legend = TRUE)+
    scale_linewidth(range = c(0.2, 2.0))+
    
    # observed sample points colored by observed values (no extra legend)
    geom_sf(data = obs_sf,
            aes(color = !!sym(obs_col)),
            shape = 21, fill = "white",
            size = 0.5, stroke = 0.5,
            inherit.aes = FALSE,
            show.legend = FALSE) +
    
    # shared continuous colorbar
    scale_color_gradientn(
      colors = pal, #the colour palette I made earlier
      limits = c(scale_min, scale_max),
      oob = scales::squish,
      na.value = "grey80",
      name = "PFAS (ng/L)",
      breaks = c(0, 1, 5,10, 15,25, 35),
      labels = scales::number_format(accuracy = 1.0),
      trans = "sqrt" #what is this?
      #trans means "Deprecated in favour of transform".. What?
    ) +
    
    guides(color = guide_colorbar(
      barwidth  = grid::unit(0.6, "cm"),
      barheight = grid::unit(6, "cm"),
      label.theme = element_text(size = 12),
      title.theme = element_text(size = 13, face = "bold"),
      title.position = "top"
    )) +
    
    labs(title = paste("Predicted", compound), subtitle = r2_label) +
    
    theme_classic() +
    theme(
      plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title    = element_text(size = 12),
      axis.text     = element_text(size = 11),
      legend.title  = element_text(size = 13, face = "bold"),
      legend.text   = element_text(size = 12),
      legend.position = "right",
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(compound, "_map.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
}
# Plots are created for all of your compounds. Check out your folder!!


#---------Summaries and Post Analysis------------------

#this next step reports all your models back to you

coeffs_df <- map_df(models, ~tidy(.x), .id = "compound")
glance_df <- map_df(models, glance, .id = "compound")
full_model_summary <- coeffs_df %>%
  left_join(glance_df, by = "compound") %>%
  select(compound, term, estimate, std.error, p.value,
         n, p, npar, value, AIC, AICc, logLik, deviance, pseudo.r.squared)
full_model_summary


#not standardized, so coefficient values are all over the place.
standardize_ssn_coefs <- function(model, data, response) {
  # raw coefficients
  raw_coefs <- coef(model)
  
  # get SD of response
  sd_y <- sd(data[[response]], na.rm = TRUE)
  
  # predictors (exclude intercept)
  preds <- names(raw_coefs)[names(raw_coefs) != "(Intercept)"]
  
  # compute standardized coefs for predictors
  std_coefs <- sapply(preds, function(p) {
    sd_x <- sd(data[[p]], na.rm = TRUE)
    raw_coefs[p] * (sd_x / sd_y)
  })
  
  # return intercept + standardized coefs
  c("(Intercept)" = raw_coefs["(Intercept)"], std_coefs)
}
std_coef_list <- list()

for (compound in pfas_cols) {
  std_coef_list[[compound]] <- standardize_ssn_coefs(
    model = models[[compound]],
    data  = NE_ssn$obs,
    response = compound
  )
}
std_coef_df <- bind_rows(lapply(names(std_coef_list), function(comp) {
  tibble(
    compound = comp,
    predictor = names(std_coef_list[[comp]])[-1],   # drop intercept
    std_coef = as.numeric(std_coef_list[[comp]][-1])
  )
}))

ggplot(std_coef_df, aes(x = predictor, y = compound, fill = std_coef)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    midpoint = 0,
    low = "blue",
    high = "red",
    mid = "white"
  ) +
  labs(
    title = "Standardized Coefficients for PFAS Models",
    x = "Predictor",
    y = "Compound",
    fill = "Std. Coef"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8)
  )

