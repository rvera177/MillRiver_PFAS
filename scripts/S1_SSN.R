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

PFAS_Spatial_Oct_2025 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Spatial_1_PFASResults.csv")
S1 <- PFAS_Spatial_Oct_2025 #need dplyr loaded in order to pipe. Don't forget lol

S1 <- S1 %>%#remove NA's for coordinates
  filter(!is.na(Lat), !is.na(Long))

#Naming family groups that I can add together
pfas_groups <- list(
  PFCA = c("PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA",
    "PFNA", "PFDA", "PFUnA", "PFDoA", "PFTrDA", "PFTeDA"),
  PFSA = c("PFBS", "PFPeS", "PFHxS", "PFHpS", "PFOS",
    "PFNS", "PFDS", "PFDoS"),
  FTSA = c("4:2FTS", "6:2FTS", "8:2FTS"),
  PFOSA = c("PFOSA", "NMeFOSA", "NEtFOSA"),
  FOSAA = c("NMeFOSAA", "NEtFOSAA"),
  PFECA = c("HFPO-DA", "ADONA", "PFMPA", "PFMBA", "NFDHA"),
  PFESA = c("9Cl-PF3ONS", "11Cl-PF3OUdS", "PFEESA"),
  FTCA = c("3-3 FTCA", "5-3 FTCA", "7-3 FTCA"))

#create new columns in dataset with Sum of corresponding species in the family
S1 <- S1 %>%
  bind_cols(
    lapply(names(pfas_groups), function(fam_name) {
      
      # PFBA → "PFBA Results"
      cols <- paste0(pfas_groups[[fam_name]], " Results")
      cols_existing <- intersect(cols, names(S1))
      
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
  select(PFAS40, PFCA, PFSA, `PFOA Results`, connws, npdesdensws, pctimp2019ws, pcturbhi2019ws, pcturblo2019ws, pcturbmd2019ws, pcturbop2019ws) %>% 
  cor(method = "spearman", use = "complete.obs")
cor_long_ws <- melt(cor_results_ws)

cor_results_Cat <- S1_Cat %>%
  select(PFAS40, PFCA, PFSA, `PFOA Results`, conncat, npdesdenscat, pctimp2019cat, pcturbhi2019cat, pcturblo2019cat, pcturbmd2019cat, pcturbop2019cat) %>% 
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


#keeping pctimp2019ws and npdesdensws
#percent impervious 2019 and NPDES Density (potential point sources)
#https://www.epa.gov/npdes
#maybeee including connws after the model is built

#I'm going to go with PFOA
#okay, now get the same two covariates for my prediction points. 

#Prediction_Points <- read_csv("~/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape/Prediction_Points.csv")
#head(Prediction_Points)

#pred_streamcat <- sc_get_data(
#  comid = Prediction_Points$COMID,
#  metric = c("conn", "npdesdens", "pctimp2019", "pcturblo2019"),
#  aoi = "ws")
#pred_streamcat <- pred_streamcat %>%
#  rename(COMID = comid)

#Prediction_Points <- left_join(Prediction_Points, pred_streamcat, by = "COMID")

#write_xlsx(Prediction_Points, "Prediction_CovariatesAdded.xlsx")
#write_xlsx(S1, "Observation_CovariatesAdded.xlsx")

#went to arcgis and created the following set up file in order to be 
#able to run the model and create the SSN object
#S1_ready <- S1 %>% mutate(Type = "Observation")
#Prediction_ready <- Prediction_Points %>% mutate(Type = "Prediction")

# Combine them
#ObservationsPredictionsMerged <- bind_rows(S1_ready, Prediction_ready)

#write_xlsx(ObservationsPredictionsMerged,
#  "~/Soil&Water lab/Spatial Stream Networks/ObservationsPredictionsMerged.xlsx")
#View(ObservationsPredictionsMerged)

#ObservationsPredictionsMerged <- read_csv("~/Soil&Water lab/Spatial Stream Networks/ObservationsPredictionsMerged.csv")
#View(ObservationsPredictionsMerged)
# this is the data I'm working with just to see. The data is actually in arcgis as a shapefile

setwd("~/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape")
#you can put whatever folder makes sense for you. 
#the location doesn't matter, it's just where final plots will be added to.
#you don't need to download any data before hand. 

#-------Using NHDPlus data set from R package ------------------
#Starting over using !
#42.38344445112979, -72.5837008452757
#coordinates at the lake warner outlet.

start_point <- st_sfc(st_point(c(-72.5837008452757, 42.38344445112979)), crs = 4269)
start_comid <- discover_nhdplus_id(start_point)
start_comid

#now that i Have the starting comid, I collect all flowlines upstream
flowline <- navigate_nldi(
  list(featureSource = "comid", featureID = start_comid),
  mode = "upstreamTributaries",
  distance_km = 200)

#save the flowlines
subset_file <- "LWMR_nhd_subset.gpkg"

#this can take a while. 
#Occasionally the NHDplus website is down, so it won't work at a given moment. 
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
plot(st_geometry(flowline), col = "blue")
plot(start_point, add = TRUE, col = "red", pch = 19)
#make sure everything is in the same projection. It should be already.

#removed 2 problematic flowlines at atkins reservoir
flowline = filter(flowline, !(objectid %in% c(47439,47740)))

#these are NSI prediction points at the center of each flowline
nsi_PredPoints_clipped <- flowline %>%
  st_centroid() %>%            # centroid works for LINESTRING
  st_cast("POINT") %>%         # ensures geometry is POINT
  mutate(comid = flowline$comid)  # carry over COMID

plot(st_geometry(flowline), col = "blue")
plot(start_point, add = TRUE, col = "red", pch = 19)
plot(st_geometry(nsi_PredPoints_clipped), add = TRUE, col = "red", pch = 19, cex = 1)
#flowlines are in!!

# Convert S1 to an sf object using Lat/Long
S1_sf <- st_as_sf(S1, 
                  coords = c("Long", "Lat"),  # x = Long, y = Lat
                  crs = 4326)                # WGS84

# 2. Transform to same CRS as flowline
S1_sf <- st_transform(S1_sf, st_crs(flowline))

# check the geometry
head(S1_sf)
plot(st_geometry(S1_sf), add = TRUE)

# Save as shapefile
st_write(S1_sf, "obs.gpkg", delete_layer = TRUE)

#error is due to time and name of a column. 
#it's still was created, so don't worry about it 



library(nngeo)  # for snapping
# I bring in my prediction points,
# and snap prediction points inside of the watershed to a flowline
#nsi_PredPoints <- st_read("C:/Users/Ruli's computer/OneDrive/Documents/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape/NE_PredictionPoints_NSI.shp")

# Transform NSI points to match flowlines CRS
#nsi_PredPoints are all the prediction points from new england
#nsi_PredPoints <- st_transform(nsi_PredPoints, st_crs(flowline))
catchment_valid <- st_make_valid(catchment)
catchment_union <- st_union(catchment_valid) #combining by subcatchments

#clip prediction points to points inside the catchment
#nsi_PredPoints_clipped <- nsi_PredPoints[st_within(nsi_PredPoints, catchment_union, sparse = FALSE), ]

obs <- st_read("obs.gpkg")  
obs <- st_transform(obs, st_crs(flowline))
obs_clip <- obs[st_within(obs, catchment_union, sparse = FALSE), ]
#obs_clip <- obs_clip %>%
#  rename(TotDASqKM = TtDASKM)%>%
#  rename(npdesdensws = npdsdns)%>%
#  rename(pctimp2019ws = pct2019)

plot(st_geometry(flowline), col = "blue")
#plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)
plot(st_geometry(catchment_union), add = TRUE, border = "black", lwd = 4)
plot(st_geometry(nsi_PredPoints_clipped), add = TRUE, col = "red", pch = 19)
plot(st_geometry(obs), add = TRUE, col = "blue", pch = 19)
#plot(st_geometry(obs_clip), add = TRUE, col = "blue", pch = 19)

#giving Prediction points obspred ID numbers and a PFAS40 row
nsi_PredPoints_clipped <- nsi_PredPoints_clipped %>% 
  mutate(
    OBSPRED = row_number() + 100000)

# get StreamCat data for the predicion point COMIDs
streamcat_data <- sc_get_data(
  comid = nsi_PredPoints_clipped$comid,
  metric = c("npdesdens", "pctimp2019", "pcturblo2019"),  # add any other metrics you want
  aoi = "ws")  # watershed scale area of interest
#add streamcat data to prediction points
nsi_PredPoints_clipped <- left_join(nsi_PredPoints_clipped, streamcat_data, by = "comid")
#shorthen prediction points variables down to only the ones i need.
nsi_PredPoints_clipped <- nsi_PredPoints_clipped %>%
  select(OBSPRED, comid, totdasqkm, npdesdensws, pctimp2019ws, geom)

#renaming obs_clip into obs
flowline <- st_transform(flowline, crs =5070)
obs <- st_transform(obs_clip, crs =5070)
pred <- st_transform(nsi_PredPoints_clipped, crs =5070)
catchment <- st_transform(catchment_union, crs =5070)


# Now plotting everything together
plot(st_geometry(flowline), col = "blue")
plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)
plot(st_geometry(obs), add = TRUE, col = "blue", pch = 19)
plot(st_geometry(pred), add = TRUE, col = "red", pch = 19)


temp_dir <- "C:/Users/Ruli's computer/OneDrive/Documents/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape"
#change this to somewhere on your computer that makes sense.
dir.create(temp_dir, showWarnings = FALSE)
library(sf)

ssn_path <- file.path(temp_dir, "PFAS_model.ssn")

flowlines_2 = lines_to_lsn(flowline, 
             lsn_path = temp_dir,
             overwrite = TRUE)
#For this catchment, the furthest observation from flowline is 138 meters
#i put snap tolerence to 150m
obs <- sites_to_lsn(
  sites = obs,
  edges = flowlines_2,
  lsn_path = temp_dir,
  file_name = "obs",
  snap_tolerance = 150,
  save_local = TRUE,
  overwrite = TRUE)
#For this catchment, the furthest prediction from flowline is 199 meters
#i put snap tolerence to 200m
preds <- sites_to_lsn(
  sites = pred,
  edges = flowlines_2,
  save_local = TRUE,
  lsn_path = temp_dir,
  file_name = "nsi_PredPoints_clipped",
  snap_tolerance = 200,
  overwrite = TRUE)
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


PFAS_ssn <- ssn_assemble(
  edges = edges,
  lsn_path = temp_dir,
  obs_sites = site.list$obs,
  preds_list = site.list[c("preds")],
  ssn_path = paste0(temp_dir, "/PFAS.ssn"),
  import = TRUE,
  check = TRUE,
  afv_col = "afvArea",
  overwrite = TRUE
)

#SSN created. Great!
# Plotting nodes for reference if you want
#nodes <- st_read(file.path(temp_dir, "nodes.gpkg"))
#nodes_proj <- st_transform(nodes, st_crs(flowline))
#plot(st_geometry(nodes_proj), add = TRUE, col = "black", pch = 19)
#text(st_coordinates(nodes_proj), labels = nodes_proj$pointid, pos = 3, cex = 0.7)
# In case you get a node_erros notification, plot them. 
#plot(st_geometry(node_errors), add=TRUE, col = "blue")


#plotting SSN
ggplot() +
  geom_sf(
    data = PFAS_ssn$edges,
    color = "medium blue",
    aes(linewidth = totdasqkm)) +
  scale_linewidth(range = c(0.5, 2.5)) +
  geom_sf(
    data = PFAS_ssn$preds$preds,
    size = 2.5,
    shape = 21,
    fill = "white",
    color = "dark grey") +
  geom_sf(
    data = PFAS_ssn$obs,
    size = 2.5,
    aes(color = PFAS40)) +
  coord_sf(datum = st_crs(obs)) +
  scale_color_viridis_c() +
  labs(color = "PFAS (ppt)", linewidth = "WS Area") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

library(SSN2)

## Generate hydrologic distance matrices
ssn_create_distmat(PFAS_ssn)


#don't fit for now!
## Fit the model
#PFAS40_mod <- ssn_lm(
#  formula = PFAS40 ~ npdesdensws + pctimp2019ws,
#  ssn.object = PFAS_ssn,
#  tailup_type = "exponential",
#  euclid_type = "gaussian",
#  additive = "afvArea")

#fitted ssn model statistics
#summary(PFAS40_mod)
#varcomp(PFAS40_mod)  #nugget is what is not being explained by the covariates

#making a copy to a temporary directory so I'm not editing original ssn
path <- system.file("temp_dir/PFAS.ssn", package = "SSN2")

PFAS_pred <- ssn_import(
  path = PFAS_ssn$path,
  predpts = c("preds"),
  overwrite = TRUE
)
summary(PFAS_pred)

names(PFAS_pred$preds)

ggplot() +
  geom_sf(data = PFAS_pred$edges) +
  geom_sf(data = PFAS_pred$preds$preds, pch = 17, color = "red") +
  geom_sf(data = PFAS_pred$obs, color = "blue", size = 2) +
  theme_bw()

#hydrologic distance matrices that preserve directionality, 
#which are required for statistical modeling
ssn_create_distmat(
  ssn.object = PFAS_pred,
  predpts = c("preds"),
  among_predpts = TRUE,
  overwrite = TRUE)

ggplot() +
  geom_sf(data = catchment, color= "black", lwd = 1.5) +
  geom_sf(data = PFAS_pred$edges) +
  geom_sf(data = PFAS_pred$obs, aes(color = PFAS40), size = 5) +
  scale_color_viridis_c(limits = c(0, 60), option = "H") +
  theme_bw()

Togregram <- Torgegram(
  formula = PFAS40 ~ npdesdensws + pctimp2019ws,
  ssn.object = PFAS_pred,
  type = c("flowcon", "flowuncon", "euclid")
)
plot(Togregram)

#quick modification of an ssn including tail down
#name the columns you want models for!
# get all column names
all_cols <- names(PFAS_ssn$obs)
#view(all_cols)
# find the positions of the start and end columns 
# that you want to make models for
start <- match("PFAS40", all_cols)
end   <- match("FTCA", all_cols)
# subset the column names
pfas_cols <- all_cols[start:end]
pfas_cols
#now i'm looking at all PFAS compounds
#so fun

models <- list()
skipped_compounds <- c()  # keep track of skipped compounds

#this takes my computer around 40 seconds to run a spatial model for ALL compounds. Sick.
for (compound in pfas_cols) {
  response <- PFAS_ssn$obs[[compound]]  # extract the response values
  # Check variability: at least 2 unique, non-NA values
  if(length(unique(na.omit(response))) < 2) {
    message(paste("Skipping", compound, "- not enough variability"))
    skipped_compounds <- c(skipped_compounds, compound)
    next  # skip this compound
  }
  
  # Fit model for each individual compound
  form <- as.formula(paste0(compound, " ~ npdesdensws + pctimp2019ws"))
  models[[compound]] <- ssn_lm(
    formula = form,
    ssn.object = PFAS_ssn,
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

#predict the compound concentration at each edge 

preds <- list()
for (compound in pfas_cols) {
  preds[[compound]] <- augment(
    models[[compound]],
    newdata = "preds",  # <-- tells augment to use all prediction sites in PFAS_ssn
    pred.type = "preds"
  )
}

for (compound in pfas_cols) {
  
  df <- preds[[compound]] %>%
    as.data.frame() %>%
    select(-comid) %>%  # remove any existing comid
    left_join(
      PFAS_ssn$preds$preds %>% st_drop_geometry() %>% select(pid, comid),
      by = "pid"
    ) %>%
    select(comid, .fitted) %>%
    rename(pred = .fitted)
  
  PFAS_pred$edges <- PFAS_pred$edges %>%
    left_join(df, by = "comid") %>%
    rename(!!paste0(compound, "_pred") := pred)
}


#should be something like this
#"geometry" ".fitted" ".resid" ".hat" ".std.resid"
out_dir <- "PFAS_maps"
if (!dir.exists(out_dir)) dir.create(out_dir)

for (compound in pfas_cols) {
  
  pred_col <- paste0(compound, "_pred")
  
  p <- ggplot() +
    geom_sf(data = catchment, fill = NA, color = "darkgreen", size = 1)+
    geom_sf(data = PFAS_pred$edges,
            aes(color = !!sym(pred_col)),
            linewidth = 1) +
    scale_color_viridis_c(option = "plasma") +
    labs(
      title = paste("Predicted", compound),
      color = "Predicted PFAS"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(out_dir, paste0(compound, "_map.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
}
