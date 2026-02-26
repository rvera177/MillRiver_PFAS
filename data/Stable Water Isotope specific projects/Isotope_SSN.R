#Raul Vera
#Created Feb 25 2026
#Last updated: 

#Hello and welcome to my PFAS SSN model
#this script models 17 PFAS stream observations from Amherst
#Spatial 1 corresponds to October 12th sampling 

#first, set your working directory. This is where your files will be saved
#getwd = current working directory/folder. 
#getwd()

#if you like this location, then leave as is.
#If you want to change the folder, use setwd
setwd("~/Soil&Water lab/Spatial Stream Networks/LWMR Isoscape")
#you can put whatever folder makes sense for you. 
#the location doesn't matter, it's just where final plots will be added to.
#you don't need to download any data before hand. 


#download all of these packages and then load
library(remotes)
library(StreamCatTools)
library(dplyr)
library(readr) 
library(sf)
library(nhdplusTools) #pulls in flowlines and associated COMIDS
library(SSN2)
library(SSNbler)
library(writexl)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(viridis)
library(stringr)
library(nngeo)

Isotope <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Stable%20Water%20Isotope%20specific%20projects/Isotope_SSN.csv")
S1 <- Isotope #need dplyr loaded in order to pipe. Don't forget lol


S1 <- S1 %>%
  mutate(OBSPRED_ID = row_number()) #creates a new column with row number


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
#this next step runs your sites through the StreamCAt database, and returns information for your metrics of interest
#there are two AOI's you can choose. 
#Choose wisely...
streamcat_data_cat <- sc_get_data( comid = S1$COMID, 
                                   metric = c("conn", "npdesdens", "pctimp2019",
                                              "pcturbhi2019", "pcturblo2019",
                                              "pcturbmd2019", "pcturbop2019",
                                              "huden2010", "rdcrs"),
                                   aoi = "Cat" ) #area of interst, imediate stream segment 
streamcat_data_ws <- sc_get_data( comid = S1$COMID, 
                                  metric = c("conn", "npdesdens", "pctimp2019",
                                             "pcturbhi2019", "pcturblo2019",
                                             "pcturbmd2019", "pcturbop2019",
                                             "huden2010", "rdcrs"),
                                  aoi = "ws" ) #area of interst, imediate stream segmetn and everything upstream of it

#renaming since there was a name conversion using sc_get_data
streamcat_data_cat <- streamcat_data_cat %>% rename(COMID = comid)
S1_Cat <- left_join(S1, streamcat_data_cat, by ="COMID") 

#Make sure that your working model uses S1!
streamcat_data_ws <- streamcat_data_ws %>% rename(COMID = comid)
S1 <- left_join(S1, streamcat_data_ws, by ="COMID") 

cor_results_ws <- S1 %>%
  select(TOC, DOC, Calcium, Sulfate, connws, npdesdensws, pctimp2019ws, pcturbhi2019ws, pcturblo2019ws, pcturbmd2019ws, pcturbop2019ws) %>% 
  cor(method = "spearman", use = "complete.obs")
cor_long_ws <- melt(cor_results_ws)

cor_results_Cat <- S1_Cat %>%
  select(TOC, DOC, Calcium, Sulfate, conncat, npdesdenscat, pctimp2019cat, pcturbhi2019cat, pcturblo2019cat, pcturbmd2019cat, pcturbop2019cat) %>% 
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
  select(TOC, DOC, Calcium, Sulfate, connws, npdesdensws, pctimp2019ws, 
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

# I'm keeping pctimp2019ws and npdesdensws
#percent impervious 2019 (non-point sources) and NPDES Density (potential point sources)
#https://www.epa.gov/npdes
#this was done based on decently ok corelation coeffictions, and my own intuition 
#based on the literature and common sense


#-------Using the NHDPlus data set from R package ------------------
#Starting over using following coordinates!
#42.38344445112979, -72.5837008452757
#coordinates at the lake warner outlet.

start_point <- st_sfc(st_point(c(-72.5837008452757, 42.38344445112979)), crs = 4269)
start_comid <- discover_nhdplus_id(start_point)
start_comid

#now that i Have the starting comid, I collect all flowlines upstream
flowline <- navigate_nldi(
  list(featureSource = "comid", featureID = start_comid),
  mode = "upstreamTributaries",
  distance_km = 200) #goes 200km upstream from starting point
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
plot(st_geometry(flowline), col = "blue")
plot(start_point, add = TRUE, col = "red", pch = 19)
plot(waterbody, add = TRUE, col = "blue")
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
#flowlines and prediction points are in!!

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

#warning message is due to time and name of a column. 
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
#st_write(obs_clip, "obs_clip.gpkg", delete_dsn = TRUE)


plot(st_geometry(flowline), add=TRUE, col = "blue")
#plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)
plot(st_geometry(catchment_union), add = TRUE, border = "black", lwd = 4)
plot(st_geometry(nsi_PredPoints_clipped), add = TRUE, col = "red", pch = 19)
plot(st_geometry(obs), add = TRUE, col = "blue", pch = 19)
#plot(st_geometry(obs_clip), add = TRUE, col = "blue", pch = 19)

#giving Prediction points obspred ID numbers and a PFAS40 row
#predictions have obspred_ID starting at 100,000. That's just the standard
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

#renaming obs_clip into obs with different CRS (in meters)
flowline <- st_transform(flowline, crs =5070)
obs <- st_transform(obs_clip, crs =5070)
pred <- st_transform(nsi_PredPoints_clipped, crs =5070)
catchment <- st_transform(catchment_union, crs =5070)


# Now plotting everything together
plot(st_geometry(flowline), col = "blue")
plot(st_geometry(catchment), add = TRUE, border = "darkgreen", lwd = 2)
plot(st_geometry(obs), add = TRUE, col = "blue", pch = 19)
plot(st_geometry(pred), add = TRUE, col = "red", pch = 19)


temp_dir <- "/Users/katarzynawisnauckas"
#change this to somewhere on your computer that makes sense.


#-----DON't Change ANyThing HERE!!--------------
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
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
  snap_tolerance = 300,
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

#names(site.list) ## View output site.list names

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

####-------Okay now you are allowed to change things----------------
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
    aes(color = TOC)) +
  coord_sf(datum = st_crs(obs)) +
  scale_color_viridis_c() +
  labs(color = "TOC (mg/L)", linewidth = "WS Area") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10))

library(SSN2)

## Generate hydrologic distance matrices
ssn_create_distmat(PFAS_ssn)

#nugget is what is not being explained by the covariates

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
  geom_sf(data = PFAS_pred$obs, aes(color = TOC), size = 5) +
  scale_color_viridis_c(limits = c(0, 6), option = "H") +
  theme_bw()

Togregram <- Torgegram(
  formula = TOC ~ npdesdensws + pctimp2019ws,
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
start <- match("TOC", all_cols)
end   <- match("Sulfate", all_cols)
# subset the column names
pfas_cols <- c("TOC","DOC", "UV254")
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


#predict the compound concentration at each edge aka stream segment

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

#save plots to the working director, in a new foler called out_dir
out_dir <- "PFAS_maps"
if (!dir.exists(out_dir)) dir.create(out_dir)

library(ggplot2)
library(sf)
library(dplyr)
library(wesanderson)
library(rlang)
library(scales)
library(broom)


# ensure obs_sf is an sf and matches edges CRS (as before)
if (!inherits(PFAS_ssn$obs, "sf")) {
  obs_sf <- st_as_sf(PFAS_ssn$obs)
} else {
  obs_sf <- PFAS_ssn$obs
}
edges_crs <- st_crs(PFAS_pred$edges)
if (is.na(st_crs(obs_sf)) || st_crs(obs_sf) != edges_crs) {
  obs_sf <- st_transform(obs_sf, edges_crs)
}
if (is.na(st_crs(catchment)) || st_crs(catchment) != edges_crs) {
  catchment <- st_transform(catchment, edges_crs)
}

# --- 1) Overview map (blacked out) with combined map key ---
edges_overview <- PFAS_pred$edges %>% mutate(feature = "Predicted ng/L")
obs_overview   <- obs_sf %>% mutate(feature = "Observed ng/L")

p_overview <- ggplot() +
  geom_sf(data = catchment, fill = "grey95", color = "black", size = 0.6) +
  geom_sf(data = edges_overview, aes(color = feature), linewidth = 0.6, inherit.aes = FALSE) +
  geom_sf(data = obs_overview, aes(shape = feature, fill = feature), size = 3, stroke = 0.6, inherit.aes = FALSE) +
  scale_color_manual(name = "Map key",
                     values = c("Predicted mg/L" = "black"),
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(name = "Map key",
                     values = c("Observed mg/L" = 21),
                     guide = guide_legend(order = 1)) +
  scale_fill_manual(name = "Map key",
                    values = c("Observed mg/L" = "black"),
                    guide = guide_legend(order = 1)) +
  labs(title = "Sample sites and stream segments") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )
p_overview 

ggsave(filename = file.path(out_dir, "overview_map.png"), plot = p_overview,
       width = 8, height = 6, dpi = 300)

#--- forloop plots for all compounds------------------

# fixed scale limits
scale_min <- 1.5
scale_max <- 6 #this is the maximum concentration predicted and observed. 
#edit the max if needed in future model creations

# make sure obs_sf is an sf and matches edges CRS
if (!inherits(PFAS_ssn$obs, "sf")) obs_sf <- st_as_sf(PFAS_ssn$obs) else obs_sf <- PFAS_ssn$obs
edges_crs <- st_crs(PFAS_pred$edges)
if (is.na(st_crs(obs_sf)) || st_crs(obs_sf) != edges_crs) obs_sf <- st_transform(obs_sf, edges_crs)
if (is.na(st_crs(catchment)) || st_crs(catchment) != edges_crs) catchment <- st_transform(catchment, edges_crs)

glance_df <- map_df(models, glance, .id = "compound") # precompute glance table for R^2 lookup

# this is the colour palette for legend (dark blue and wes anderson colors)
pal <- colorRampPalette(c("#012A4A", wes_palette("Zissou1", type = "continuous")))(256)

for (compound in pfas_cols) {
  if (!compound %in% names(models)) next
  
  pred_col <- paste0(compound, "_pred")
  obs_col  <- compound
  
  r2_val <- glance_df %>% # get R^2 or pseudo-R^2
    filter(compound == !!compound) %>%
    { if (nrow(.) == 0) NA_real_ else
      if ("pseudo.r.squared" %in% names(.)) .$pseudo.r.squared else
        if ("r.squared" %in% names(.)) .$r.squared else NA_real_ }
  r2_label <- if (is.na(r2_val)) "R² = NA" else paste0("R\u00B2 = ", format(r2_val, digits = 3))
  
  p <- ggplot() +
    geom_sf(data = catchment, fill = NA, color = "black", size = 0.6) +
    # predicted stream segments (continuous color)
    geom_sf(data = PFAS_pred$edges,
            aes(color = !!sym(pred_col)),
            linewidth = 1.5, show.legend = TRUE) +
    # sample locations colored by observed values
    geom_sf(data = obs_sf,
            aes(color = !!sym(obs_col)),
            shape = 21, fill = "white",
            size = 2, stroke = 2.5,
            inherit.aes = FALSE,
            show.legend = FALSE) +
    
    #continuous colorbar for all plots
    scale_color_gradientn(
      colors = pal,
      limits = c(scale_min, scale_max),
      oob = scales::squish,
      na.value = "grey80",
      name = "PFAS (ng/L)",
      breaks = c(0, 1,2,3,4,5),
      labels = scales::number_format(accuracy = 1.0),
      trans = "sqrt") +
    
    guides(color = guide_colorbar(
      barwidth  = grid::unit(0.6, "cm"),
      barheight = grid::unit(6, "cm"),
      label.theme = element_text(size = 12),
      title.theme = element_text(size = 13, face = "bold"),
      title.position = "top" )) +
    
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
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6))
  
  ggsave(
    filename = file.path(out_dir, paste0(compound, "_map.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
}

p
#this next step reports all your model summaries back to you
coeffs_df <- map_df(models, ~tidy(.x), .id = "compound")
glance_df <- map_df(models, glance, .id = "compound")
full_model_summary <- coeffs_df %>%
  left_join(glance_df, by = "compound") %>%
  select(compound, term, estimate, std.error, p.value,
         n, p, npar, value, AIC, AICc, logLik, deviance, pseudo.r.squared)
full_model_summary
#for comparing multiple mods. You don't need to worry about this, these are notes to self
#glances(PFAS40_mod, PFCA_mod, PFSA_mod, PFOA_mod)
#tidy(PFCA_mod, conf.int = TRUE)
#glance(PFCA_mod)
#logLik is the log likelihood
#plot(PFCA_mod, which = 1)
#need to use the lowest AIC and AICc,

#not standardized, so coefficient values are all over the place.
#this next step will standardize everything.
standardize_ssn_coefs <- function(model, data, response) {
  raw_coefs <- coef(model) # raw coefficients numbers
  sd_y <- sd(data[[response]], na.rm = TRUE)  # get SD of response
  preds <- names(raw_coefs)[names(raw_coefs) != "(Intercept)"]  # predictors (exclude intercept)
  
  #Now you standardize the coefs for predictors
  std_coefs <- sapply(preds, function(p) {
    sd_x <- sd(data[[p]], na.rm = TRUE)
    raw_coefs[p] * (sd_x / sd_y)
  })
  # return the intercept + standardized coefs
  c("(Intercept)" = raw_coefs["(Intercept)"], std_coefs)
}
std_coef_list <- list()

for (compound in pfas_cols) {
  std_coef_list[[compound]] <- standardize_ssn_coefs(
    model = models[[compound]],
    data  = PFAS_ssn$obs,
    response = compound
  )
}

std_coef_df <- bind_rows(lapply(names(std_coef_list), function(comp) {
  tibble(
    compound = comp,
    predictor = names(std_coef_list[[comp]])[-1],   # drop intercept
    std_coef = as.numeric(std_coef_list[[comp]][-1]))
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
    x = "Predictor Metric",
    y = "Modeled Compound",
    fill = "Std. Coef"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8)
  )



#That's alot of models, so i'll make a heat map for only a few that would be interesting in talking about
#want_comps are the compounds or families to put on the figure. 
want_comps <- c("PFCA", "PFSA")

# renaming the prediction labels for easier readibility
pred_labels <- c(
  npdesdensws   = "NDPES Density",
  pctimp2019ws  = "Percent Impervious")

# preparing plotting set up
plot_df <- std_coef_df %>%
  filter(compound %in% want_comps) %>%
  mutate(
    predictor = as.character(predictor),
    predictor = sub("^([^\\.]+)\\..*$", "\\1", predictor),
    compound = factor(compound, levels = want_comps)
  ) %>%
  mutate(
    predictor_label = ifelse(predictor %in% names(pred_labels),
                             pred_labels[predictor],
                             predictor),
    predictor_label = factor(predictor_label, levels = unique(predictor_label)),
    label = sprintf("%.2f", std_coef),
    text_color = ifelse(abs(std_coef) > 0.5, "white", "black")
    #text_color will update based on color of the background tile
  )

# plotting coefficient tiles with labels on top
ggplot(plot_df, aes(x = predictor_label, y = compound, fill = std_coef)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label, color = text_color), size = 10) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Std coef") +
  labs(title = "Standardized
model coefficients",
       x = "Predictor metric", y = "Modeled Compound(s)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 25, face = "bold", hjust=0.5)
  )