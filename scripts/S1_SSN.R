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


PFAS_Spatial_Oct_2025 <- read_csv("~/Soil&Water lab/WetCenterPFASResults.csv")
S1 <- PFAS_Spatial_Oct_2025 #need dplyr loaded in order to pipe. Don't forget lol

S1 <- S1 %>%#remove NA's for coordinates
  filter(!is.na(Lat), !is.na(Long))

#snapping sites to nearest flowline so it gets the nearest comid
sites <- st_as_sf(
  S1,
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
streamcat_data <- sc_get_data( comid = S1$COMID, 
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
streamcat_data <- streamcat_data %>% rename(COMID = comid)
S1 <- left_join(S1, streamcat_data, by ="COMID") 

streamcat_data_ws <- streamcat_data_ws %>% rename(COMID = comid)
S1_ws <- left_join(S1, streamcat_data_ws, by ="COMID") 

cor_results_ws <- S1_ws %>%
  select(PFAS40, `PFOA Results`, connws, npdesdensws, pctimp2019ws, pcturbhi2019ws, pcturblo2019ws, pcturbmd2019ws, pcturbop2019ws) %>% 
  cor(method = "spearman", use = "complete.obs")
cor_long_ws <- melt(cor_results_ws)

cor_results_Cat <- S1 %>%
  select(PFAS40, `PFOA Results`, conncat, npdesdenscat, pctimp2019cat, pcturbhi2019cat, pcturblo2019cat, pcturbmd2019cat, pcturbop2019cat) %>% 
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

Prediction_Points <- read_csv("~/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape/Prediction_Points.csv")
head(Prediction_Points)

pred_streamcat <- sc_get_data(
  comid = Prediction_Points$COMID,
  metric = c("conn", "npdesdens", "pctimp2019"),
  aoi = "ws")
pred_streamcat <- pred_streamcat %>%
  rename(COMID = comid)

Prediction_Points <- left_join(Prediction_Points, pred_streamcat, by = "COMID")

write_xlsx(Prediction_Points, "Prediction_CovariatesAdded.xlsx")
write_xlsx(S1, "Observation_CovariatesAdded.xlsx")

#went to arcgis and created the following set up file in order to be 
#able to run the model and create the SSN object

ObservationsPredictionsMerged <- read_csv("~/Soil&Water lab/Spatial Stream Networks/ObservationsPredictionsMerged.csv")
View(ObservationsPredictionsMerged)
# this is the data I'm working with just to see. The data is actually in arcgis as a shapefile

setwd("~/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape")
gdb_path <- "C:/Users/Ruli's computer/OneDrive/Documents/Soil&Water lab/Spatial Stream Networks/rmrs-preparing-input-data-for-generating-the-ssnobject_exampledata/MyProject/MyProject.gdb"
st_layers("C:/Users/Ruli's computer/OneDrive/Documents/Soil&Water lab/Spatial Stream Networks/rmrs-preparing-input-data-for-generating-the-ssnobject_exampledata/MyProject/MyProject.gdb")
ssn_dir <- "C:/Users/Ruli's computer/OneDrive/Documents/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape"

OBSPred_sf <- st_read(dsn = gdb_path, layer = "ObservationsPred_Merge")
# Separate into obs (1:17) and preds (18:47)
obs <- OBSPred_sf %>% filter(OBSPRED_ID <= 17)
preds <- OBSPred_sf %>% filter(OBSPRED_ID > 17)

# Write them as new shapefiles
st_write(obs, "obs.shp", delete_layer = TRUE)
st_write(preds, "preds.shp", delete_layer = TRUE)

# Read edges and points
edges <- st_read("LWMR_StreamNetwork_NSI.shp")  # projected
obs_sites <- st_read("obs.shp")  # only observations
pred_sites <- st_read("preds.shp")

# Transform points to match edges CRS
obs_sites <- st_transform(obs_sites, st_crs(edges))
pred_sites <- st_transform(pred_sites, st_crs(edges))

#Convert flowlines edges to LSN (network)
lines_to_lsn(
  streams = edges,
  lsn_path = ssn_dir,
  overwrite = TRUE)


ssn.obj <- ssn_assemble(
  edges = edges,
  lsn_path = ssn_dir,
  obs_sites = obs_sites,
  preds_list = list(preds = pred_sites),
  ssn_path = file.path(ssn_dir, "PFAS_model.ssn"),
  import = TRUE,
  overwrite = TRUE)

#somethings wrong. 

#-------Using NHDPlus data set from R ------------------
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
#make objects for eacth think and plot
flowline  <- sf::read_sf(subset_file, "NHDFlowline_Network")
catchment <- sf::read_sf(subset_file, "CatchmentSP")
waterbody <- sf::read_sf(subset_file, "NHDWaterbody")
plot(st_geometry(flowline), col = "blue")
plot(start_point, add = TRUE, col = "red", pch = 19)

library(nngeo)  # for snapping
# I bring in my prediction points,
# and snap prediction points inside of the watershed to a flowline
nsi_PredPoints <- st_read("C:/Users/Ruli's computer/OneDrive/Documents/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape/NE_PredictionPoints_NSI.shp")

# 2. Transform NSI points to match flowlines CRS
nsi_PredPoints <- st_transform(nsi_PredPoints, st_crs(flowline))
catchment_valid <- st_make_valid(catchment)
catchment_union <- st_union(catchment_valid) #combining by subcatchments

#clip prediction points to points inside the catchment
nsi_PredPoints_clipped <- nsi_PredPoints[st_within(nsi_PredPoints, catchment_union, sparse = FALSE), ]

obs <- st_read("obs.shp")  
obs <- st_transform(obs, st_crs(flowline))
obs_clip <- obs[st_within(obs, catchment_union, sparse = FALSE), ]
obs_clip <- obs_clip %>%
  rename(TotDASqKM = TtDASKM)%>%
  rename(npdesdensws = npdsdns)%>%
  rename(pctimp2019ws = pct2019)

plot(st_geometry(flowline), col = "blue")
#plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)
plot(st_geometry(catchment_union), add = TRUE, border = "black", lwd = 4)
plot(st_geometry(nsi_PredPoints_clipped), add = TRUE, col = "red", pch = 19)
plot(st_geometry(obs_clip), add = TRUE, col = "blue", pch = 19)

#giving Prediction points obspred ID numbers and a PFAS40 row
nsi_PredPoints_clipped <- nsi_PredPoints_clipped %>% 
  mutate(
    OBSPRED = row_number() + 100000,
    PFAS40 = 0)

# get StreamCat data for the predicion point COMIDs
streamcat_data <- sc_get_data(
  comid = nsi_PredPoints_clipped$COMID,
  metric = c("npdesdens", "pctimp2019"),  # add any other metrics you want
  aoi = "ws"  # watershed
) %>%
  rename(COMID = comid)
#add streamcat data to prediction points
nsi_PredPoints_clipped <- left_join(nsi_PredPoints_clipped, streamcat_data, by = "COMID")
#shorthen prediction points variables down to only the ones i need.
nsi_PredPoints_clipped <- nsi_PredPoints_clipped %>%
  select(OBSPRED, PFAS40, COMID, TotDASqKM, npdesdensws, pctimp2019ws, geometry)


#make sure everything is in the same projection. It should be already.
flowline <- st_transform(flowline, crs = 5070)
flowline = filter(flowline, !(objectid %in% c(47439,47740)))


obs <- st_transform(obs_clip, crs = 5070)
pred <- st_transform(nsi_PredPoints_clipped, crs = 5070)
catchment <- st_transform(catchment_union, crs = 5070)
# Snap within 100 m 
obs_snapped_geom <- st_snap(st_geometry(obs), st_geometry(flowline), tolerance = 100)
pred_snapped_geom <- st_snap(st_geometry(pred), st_geometry(flowline), tolerance = 100)

obs_snapped  <- obs
pred_snapped <- pred

st_geometry(obs_snapped)  <- obs_snapped_geom
st_geometry(pred_snapped) <- pred_snapped_geom

plot(st_geometry(flowline), col = "blue")
plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)
plot(st_geometry(obs_snapped), add = TRUE, col = "blue", pch = 19)
plot(st_geometry(pred_snapped), add = TRUE, col = "red", pch = 19)


temp_dir <- "C:/Users/Ruli's computer/OneDrive/Documents/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape"
dir.create(temp_dir, showWarnings = FALSE)
library(sf)

ssn_path <- file.path(temp_dir, "PFAS_model.ssn")

flowlines_2 = lines_to_lsn(flowline, 
             lsn_path = temp_dir,
             overwrite = TRUE)

obs <- sites_to_lsn(
  sites = obs_snapped,
  edges = flowlines_2,
  lsn_path = temp_dir,
  file_name = "obs",
  snap_tolerance = 100,
  save_local = TRUE,
  overwrite = TRUE
)

preds <- sites_to_lsn(
  sites = pred_snapped,
  edges = flowlines_2,
  save_local = TRUE,
  lsn_path = temp_dir,
  file_name = "nsi_PredPoints_clipped",
  snap_tolerance = 100,
  overwrite = TRUE
)

edges <- updist_edges(
  edges = flowlines_2,
  save_local = TRUE,
  lsn_path = temp_dir,
  calc_length = TRUE
)

site.list <- updist_sites(
  sites = list(
    obs = obs,
    preds = preds
  ),
  edges = edges,
  length_col = "Length",
  save_local = TRUE,
  lsn_path = temp_dir
)

names(site.list) ## View output site.list names

edges <- afv_edges(
  edges = edges,
  infl_col = "totdasqkm",
  segpi_col = "areaPI",
  afv_col = "afvArea",
  lsn_path = temp_dir
)


site.list <- afv_sites(
  sites = site.list,
  edges = edges,
  afv_col = "afvArea",
  save_local = TRUE,
  lsn_path = temp_dir
)

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

#SSN created. 
# Plotting nodes for reference if you want
#nodes <- st_read(file.path(temp_dir, "nodes.gpkg"))
#nodes_proj <- st_transform(nodes, st_crs(flowline))
#plot(st_geometry(nodes_proj), add = TRUE, col = "black", pch = 19)
#text(st_coordinates(nodes_proj), labels = nodes_proj$pointid, pos = 3, cex = 0.7)
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

## Fit the model
ssn_mod <- ssn_lm(
  formula = PFAS40 ~ npdesdensws + pctimp2019ws,
  ssn.object = PFAS_ssn,
  tailup_type = "exponential",
  euclid_type = "gaussian",
  additive = "afvArea")

#fitted ssn model statistics
summary(ssn_mod)
varcomp(ssn_mod)  #nugget is what is not being explained by the covariates

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
  geom_sf(data = PFAS_pred$preds$preds, pch = 17, color = "blue") +
  geom_sf(data = PFAS_pred$obs, color = "brown", size = 2) +
  theme_bw()

#hydrologic distance matrices that preserve directionality, 
#which are required for statistical modeling
ssn_create_distmat(
  ssn.object = PFAS_pred,
  predpts = c("preds"),
  among_predpts = TRUE,
  overwrite = TRUE
)

ggplot() +
  geom_sf(data = catchment, color= "black", lwd = 1.5) +
  geom_sf(data = PFAS_pred$edges) +
  geom_sf(data = PFAS_pred$obs, aes(color = PFAS40), size = 5) +
  scale_color_viridis_c(limits = c(0, 50), option = "H") +
  theme_bw()

Togregram <- Torgegram(
  formula = PFAS40 ~ npdesdensws + pctimp2019ws,
  ssn.object = PFAS_pred,
  type = c("flowcon", "flowuncon", "euclid")
)
plot(Togregram)

#quick modification of an ssn including tail down
PFAS_mod <- ssn_lm(
  formula = PFAS40 ~ npdesdensws + pctimp2019ws,
  ssn.object = PFAS_pred,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)
tidy(PFAS_mod, conf.int = TRUE)
glance(PFAS_mod)
#logLik is the log likelihood
plot(PFAS_mod, which = 1)
tidy(PFAS_pred, conf.int = TRUE)
#comparing multiple mods and the original
ssn_mod2 <- ssn_lm(
  formula = PFAS40 ~ npdesdensws + pctimp2019ws,
  ssn.object = PFAS_pred,
  taildown_type = "spherical"
)
glances(ssn_mod, PFAS_mod, ssn_mod2)
#ssn_mod might be the best one to go with 
#need to use the lowest AIC and AICc,
#which in this case comes from ss_mod.

predict(ssn_mod, newdata = "preds")
aug_preds <- augment(ssn_mod, newdata = "preds")
aug_preds[, ".fitted"]

#this was an attempt to make edge/stream reaches same
#colour as corresponding PFAS predict point. 
edge_preds <- augment(ssn_mod, pred.type = "edge")

#saving predicted values (.fitted) to a new object
pred_vals <- aug_preds %>%
  st_drop_geometry() %>%
  select(COMID, .fitted) %>%
  distinct()
#joining to edges so it has the predicted values
PFAS_pred$edges <- PFAS_pred$edges %>%
  left_join(pred_vals,
            by = c("comid" = "COMID"))
st_write(PFAS_pred$edges, "PFAS_predicted_edges.gpkg", delete_dsn = TRUE)


ggplot() +
  geom_sf(data = catchment, color= "black", lwd = 1.5, fill="black") +
  geom_sf(data = PFAS_pred$edges, aes(color = .fitted, linewidth=2)) +
  scale_linewidth(range = c(0.7, 2))+
  geom_sf(data = PFAS_pred$obs, aes(fill = PFAS40), color = "white", size = 5, stroke = 0.5, pch =19)+
  #geom_sf(data = aug_preds, aes(color = .fitted), size = 4) +
  scale_color_viridis_c(limits = c(0, 61), option = "H") +
  labs(
    title = "PFAS Spatial Stream Network",
    subtitle = "Mill River Watershed: Amherst, MA",
    color = "Predicted PFAS (ng/L)" )+
  theme_classic()

summary(PFAS_pred)


st_write(aug_preds, paste0(tempdir(), "/aug_preds.gpkg"))
predict(ssn_mod)
predict(ssn_mod, newdata = "all")

predict(ssn_mod, newdata = "preds", block = TRUE, interval = "prediction")

