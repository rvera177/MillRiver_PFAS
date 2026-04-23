
#Raul Vera
#Created April 18th
#Last updated: December 17th

#Hello and welcome to my PFAS SSN model for Spatial 3
#this script models 32 PFAS stream observatios from Amherst
#Spatial 3 corresponds to March 28th Amherst Mill River sampling 

#first, set your working directory. This is where your files will be saved
#getwd = current working directory/folder. 
getwd()
#if you like this location, then leave as is.
#If you want to change the folder, use setwd
setwd("C:/Users/Marston User/Documents/LWMR Isoscape")
#you can put whatever folder makes sense for you. 
#the location doesn't matter, it's just where final plots will be added to.

#download all of these packages and then load
library(remotes)
library(StreamCatTools)
library(dplyr) #need dplyr loaded in order to pipe. Don't forget 
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

PFAS_Spatial_March_2026 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Spatial_3_PFASresults.csv")
S1 <- PFAS_Spatial_March_2026 

S1 <- S1 %>%#remove NA's for coordinates
  filter(!is.na(Latitude), !is.na(Longitude))

#these changes the name of some of the headers for easier manipulation later
clean_names <- function(x) {
  x |>
    gsub(" ", "_", x = _) |>
    gsub("[-:]", ".", x = _) |>
    make.names() |>
    gsub("^X(\\d)", "\\1", x = _)
}

colnames(S1) <- clean_names(colnames(S1))

#Naming functional head groups that I can add together later
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

#create new columns in S1 with the sum of corresponding species in the family
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


# Snap sites to nearest flowline to get COMID
sites <- st_as_sf(S1,
                  coords = c("Longitude", "Latitude"),
                  crs = 4326)

bb_poly <- st_as_sfc(st_bbox(sites))
bb_sf   <- st_sf(geometry = bb_poly)
flines  <- get_nhdplus(AOI = bb_sf, realization = "flowline")

idx <- get_flowline_index(flines, sites, max_matches = 1)

# --- FIX: safe join instead of direct assignment ---
# idx contains a column 'point_id' that corresponds to the row number of `sites`
S1$COMID <- NA_real_  # initialize with NA so unmatched sites get NA

S1$COMID[idx$id] <- as.numeric(idx$COMID)

# Check which sites didn't get a COMID
missing <- which(is.na(S1$COMID))
cat("Sites without a COMID:", length(missing), "\n")
print(S1[missing, c("Sample_Name_on_LC.MS", "Name")])  # inspect the problem sites


#this bb box method works 100% of the time, 20% of the time. 


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
#renaming since there was a name conversion using sc_get_data
streamcat_data_cat <- streamcat_data_cat %>% rename(COMID = comid)
S1_Cat <- left_join(S1, streamcat_data_cat, by ="COMID") 

cor_results_Cat <- S1_Cat %>%
  select(PFAS40, PFCA, PFSA, PFOA_Results, conncat, npdesdenscat, pctimp2019cat, pcturbhi2019cat, pcturblo2019cat, pcturbmd2019cat, pcturbop2019cat) %>% 
  cor(method = "spearman", use = "complete.obs")
cor_long_cat <- reshape2::melt(cor_results_Cat)



#trying even more predictors at a watershed scale. 
# Pull a comprehensive set of StreamCat metrics all at once
streamcat_big <- sc_get_data(
  comid = S1$COMID,
  metric = c(
    # urban/developed
    "pctimp2019", "pcturbhi2019", "pcturblo2019", "pcturbmd2019", "pcturbop2019",
    "huden2010", "rdcrs",
    # point sources
    "npdesdens",
    "superfunddens",
    "septic",
    # agricultural
    "fert", "pctagdrainage", "p_ff_2017", "manure",
    # connectivity
    "conn",
    # nitrogen
    "n_hw_2012", "n_usgsww_2012", "n_dep_2017",
    # phosphorus  
    "p_ags_2017",
    # forest/natural
    "pctmxfst2019", "pctdecid2019", "pctconif2019",
    # wetlands
    "pctwdwet2019",
    # geology/soils
    "bfi",        # baseflow index
    "wtdep"      # water table depth
  ),
  aoi = "ws"
)

streamcat_big <- streamcat_big %>% rename(COMID = comid)
# Automatically identify all StreamCat predictor columns
predictor_cols <- names(streamcat_big %>% select(-COMID))
cat("Total predictors pulled:", length(predictor_cols), "\n")

S1 <- left_join(S1, streamcat_big, by = "COMID")


# Correlate all predictors against PFAS responses automatically
pfas_response_cols <- c("PFAS40", "PFCA", "PFSA", "PFOA_Results", "PFOS_Results")

cor_big <- S1 %>%
  select(all_of(c(pfas_response_cols, predictor_cols))) %>%
  cor(method = "spearman", use = "complete.obs")

# Extract only PFAS vs predictor block
cor_block <- cor_big[predictor_cols, pfas_response_cols, drop = FALSE]

# Convert to long and rank by absolute correlation with PFAS40
cor_long_big <- reshape2::melt(cor_block,
                     varnames = c("Predictor", "PFAS"),
                     value.name = "Spearman_r") %>%
  group_by(Predictor) %>%
  mutate(max_abs_r = max(abs(Spearman_r))) %>%
  ungroup() %>%
  arrange(desc(max_abs_r)) %>%
  mutate(text_color = ifelse(abs(Spearman_r) > 0.7, "white", "black"))

ggplot(cor_long_big, aes(x = PFAS, y = Predictor, fill = Spearman_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Spearman_r, 2), color = text_color), size = 3.5) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman's \u03c1"   # unicode rho symbol
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y  = element_text(hjust = 1, size = 10),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    panel.grid   = element_blank()
  )

# Print ranked table - useful for predictor selection
cor_long_big %>%
  filter(PFAS == "PFAS40") %>%
  select(Predictor, Spearman_r) %>%
  arrange(desc(abs(Spearman_r))) %>%
  print(n = Inf)

#okay so now i run a PCA on Spearman Ranked 

library(FactoMineR)
library(factoextra)

# Step 1 - select your variables
pca_vars <- S1 %>%
  st_drop_geometry() %>%
  select(
    PFAS40, PFCA, PFSA, PFOA_Results,
    npdesdensws, pctimp2019ws, fertws,
    huden2010ws, rdcrsws, connws,
    bfiws, n_hw_2012ws, septicws,
    pctagdrainagews, pctmxfst2019ws
  ) %>%
  na.omit()

# Step 2 - rank transform everything
pca_ranked <- pca_vars %>%
  mutate(across(everything(), rank))

# Step 3 - run PCA on ranked data
pca_result <- PCA(
  pca_ranked,
  scale.unit = TRUE,   # standardize to unit variance
  ncp = 5,             # keep 5 components
  graph = FALSE
)

# Step 4 - check variance explained
fviz_eig(pca_result, addlabels = TRUE) +
  labs(title = "Variance explained by each PCA component")

# Step 5 - variable loadings plot
fviz_pca_var(
  pca_result,
  col.var = "contrib",    # color by contribution
  gradient.cols = c("#012A4A", "#5DADE2", "#8B0000"),
  repel = TRUE,           # avoid label overlap
  title = "PCA variable loadings — Spearman rank transformed"
)

# Step 6 - biplot showing sites and variables
fviz_pca_biplot(
  pca_result,
  repel = TRUE,
  col.var = "#185FA5",
  col.ind = "#D85A30",
  title = "PCA biplot"
)

# Print loadings for first 3 components
print(pca_result$var$coord[, 1:3])

# Print variance explained
print(pca_result$eig)

# Visualize the separation clearly
fviz_pca_var(
  pca_result,
  axes = c(1, 2),           # PC1 vs PC2
  col.var = "contrib",
  gradient.cols = c("#012A4A", "#5DADE2", "#8B0000"),
  repel = TRUE,
  title = "PCA variable loadings (PC1 vs PC2)"
)

fviz_pca_var(
  pca_result,
  axes = c(1, 3),           # PC1 vs PC3
  col.var = "contrib",
  gradient.cols = c("#012A4A", "#5DADE2", "#8B0000"),
  repel = TRUE,
  title = "PCA variable loadings (PC1 vs PC3)"
)

library(patchwork)
library(factoextra)

# PCA biplot
p_biplot <- fviz_pca_biplot(
  pca_result,
  repel = TRUE,
  col.var = "#185FA5",      # predictors in blue
  col.ind = "#D85A30",      # sites in orange
  alpha.ind = 0.6,
  title = "PCA biplot — Spearman rank transformed",
  ggtheme = theme_classic()
)

# Or a cleaner variables-only plot
p_vars <- fviz_pca_var(
  pca_result,
  axes = c(1, 2),
  col.var = "contrib",
  gradient.cols = c("#012A4A", "#5DADE2", "#8B0000"),
  repel = TRUE,
  title = "Variable loadings — PC1 vs PC2",
  ggtheme = theme_classic()
) 

print(p_vars)
# Check loadings of your three chosen predictors across PC1 and PC2
chosen <- c("npdesdensws", "pctimp2019ws", "pctagdrainagews")
pca_result$var$coord[chosen, 1:3]

cor(S1$pctimp2019ws, S1$pctagdrainagews, method = "spearman", use = "complete.obs")

# I'm keeping pctimp2019ws and npdesdensws and agriculture!
#percent impervious 2019 (non-point sources) and NPDES Density (potential point sources)
#https://www.epa.gov/npdes
#this was done based on XYZ


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

#removed 1 problematic flowline at atkins reservoir
flowline = filter(flowline, !(objectid %in% c(47442)))

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
                  coords = c("Longitude", "Latitude"),  # x = Long, y = Lat
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


plot(st_geometry(flowline), col = "blue")
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
  metric = c("npdesdens", "pctimp2019", "pctagdrainage"),  # add any other metrics you want
  aoi = "ws")  # watershed scale area of interest
#add streamcat data to prediction points
nsi_PredPoints_clipped <- left_join(nsi_PredPoints_clipped, streamcat_data, by = "comid")
#shorthen prediction points variables down to only the ones i need.
nsi_PredPoints_clipped <- nsi_PredPoints_clipped %>%
  select(OBSPRED, comid, totdasqkm, npdesdensws, pctimp2019ws, pctagdrainagews, geom)

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


temp_dir <- "C:/Users/Marston User/Documents/LWMR Isoscape/S3"
#change this to somewhere on your computer that makes sense.


#-----DON't Change ANyThing HERE!!--------------
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
  snap_tolerance = 250,
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

####-------Okay now you are allowed to change things----------------
names(site.list$preds) ## View column names in pred
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
  formula = PFAS40 ~ npdesdensws + pctimp2019ws + pctagdrainagews,
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


#first, exploring different SSN model structure for various compounds

# Run covariance comparison for all key compounds
key_compounds <- c("PFAS40", "PFSA", "PFCA", "PFOS_Results", "PFOA_Results")

covar_results <- list()

for (compound in key_compounds) {
  
  form <- as.formula(paste0(compound, " ~ npdesdensws + pctimp2019ws + pctagdrainagews"))
  
  m_tu <- ssn_lm(form, PFAS_ssn, tailup_type = "exponential", 
                 additive = "afvArea", estmethod = "ml")
  m_td <- ssn_lm(form, PFAS_ssn, taildown_type = "exponential",
                 additive = "afvArea", estmethod = "ml")
  m_tu_td <- ssn_lm(form, PFAS_ssn, tailup_type = "exponential",
                    taildown_type = "exponential",
                    additive = "afvArea", estmethod = "ml")
  
  aic_vals <- c(TU = AIC(m_tu), TD = AIC(m_td), TU_TD = AIC(m_tu_td))
  r2_vals  <- c(
    TU    = as.numeric(summary(m_tu)$pseudoR2),
    TD    = as.numeric(summary(m_td)$pseudoR2),
    TU_TD = as.numeric(summary(m_tu_td)$pseudoR2)
  )
  
  covar_results[[compound]] <- data.frame(
    compound  = compound,
    model     = names(aic_vals),
    AIC       = round(aic_vals, 2),
    delta_AIC = round(aic_vals - min(aic_vals), 2),
    pseudo_R2 = round(r2_vals, 3)
  )
}

# Combine and print the comparison summary.
#the goal is to have a model with a small AIC and large pseudo_R2
covar_summary <- bind_rows(covar_results)
print(covar_summary)

# results show that diffent compounds show different AIC depending on Tail up, Tail Down 

models <- list()
skipped_compounds <- c()

for (compound in pfas_cols) {
  response <- PFAS_ssn$obs[[compound]]
  if(length(unique(na.omit(response))) < 2) {
    message(paste("Skipping", compound, "- not enough variability"))
    skipped_compounds <- c(skipped_compounds, compound)
    next
  }
  
  form <- as.formula(paste0(compound, " ~ npdesdensws + pctimp2019ws + pctagdrainagews"))
  models[[compound]] <- ssn_lm(
    formula = form,
    ssn.object = PFAS_ssn,
    tailup_type = "exponential",
    additive = "afvArea",
    estmethod = "ml")
} 

# Predict compound concentration at each edge
preds <- list()
for (compound in pfas_cols) {
  if (!compound %in% names(models)) {
    message("Skipping ", compound, " - no model fit")
    next
  }
  preds[[compound]] <- augment(
    models[[compound]],
    newdata = "preds",
    pred.type = "preds"
  )
}

# Remove any existing _pred columns to avoid duplicates
cols_to_remove <- grep("_pred", names(PFAS_pred$edges), value = TRUE)
cat("Removing these columns:\n")
print(cols_to_remove)
PFAS_pred$edges <- PFAS_pred$edges %>%
  select(-all_of(cols_to_remove))

# Join predictions to edges
for (compound in pfas_cols) {
  
  if (compound %in% skipped_compounds || is.null(preds[[compound]])) {
    message("Skipping ", compound)
    next
  }
  
  if (!"pid" %in% names(st_drop_geometry(preds[[compound]]))) {
    message("Skipping ", compound, " - no pid column")
    next
  }
  
  df <- preds[[compound]] %>%
    st_drop_geometry() %>%
    select(-comid) %>%
    left_join(
      PFAS_ssn$preds$preds %>% st_drop_geometry() %>% select(pid, comid),
      by = "pid"
    ) %>%
    select(comid, .fitted) %>%
    rename(!!paste0(compound, "_pred") := .fitted)
  
  PFAS_pred$edges <- PFAS_pred$edges %>%
    left_join(df, by = "comid")
}

# Verify
pred_cols_present <- grep("_pred$", names(PFAS_pred$edges), value = TRUE)
cat("Pred columns in edges:", length(pred_cols_present), "\n")




#save plots to the working director, in a new folder called out_dir
out_dir <- "PFAS_maps"
if (!dir.exists(out_dir)) dir.create(out_dir)


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
p_overview 

ggsave(filename = file.path(out_dir, "overview_map.png"), plot = p_overview,
       width = 8, height = 6, dpi = 300)

#--- forloop plots for all compounds------------------

#left off here. Title map?
# fixed scale limits
scale_min <- 0
scale_max <- 60 #this is the maximum concentration predicted and observed. 
#edit the max if needed in future model creations

# make sure obs_sf is an sf and matches edges CRS
if (!inherits(PFAS_ssn$obs, "sf")) obs_sf <- st_as_sf(PFAS_ssn$obs) else obs_sf <- PFAS_ssn$obs
edges_crs <- st_crs(PFAS_pred$edges)
if (is.na(st_crs(obs_sf)) || st_crs(obs_sf) != edges_crs) obs_sf <- st_transform(obs_sf, edges_crs)
if (is.na(st_crs(catchment)) || st_crs(catchment) != edges_crs) catchment <- st_transform(catchment, edges_crs)

glance_df <- map_df(models, glance, .id = "compound") # precompute glance table for R^2 lookup
glance_df
# this is the colour palette for legend (dark blue and wes anderson colors)
pal <- colorRampPalette(c("#012A4A", "#014F86",wes_palette("Zissou1", type = "continuous"),
                          "#8B0000"))(256)

for (compound in pfas_cols) {
  if (!compound %in% names(models)) next
  
  pred_col <- paste0(compound, "_pred")
  obs_col  <- compound
  
  r2_val <- glance_df %>%
    filter(compound == !!compound) %>%
    { if (nrow(.) == 0) NA_real_ else
      if ("pseudo.r.squared" %in% names(.)) .$pseudo.r.squared else
        if ("r.squared" %in% names(.)) .$r.squared else NA_real_ }
  r2_label <- if (is.na(r2_val)) "R² = NA" else paste0("R\u00B2 = ", format(r2_val, digits = 3))
  
  p <- ggplot() +
    geom_sf(data = catchment, fill = NA, color = "black", size = 0.6) +
    geom_sf(data = PFAS_pred$edges,
            aes(color = .data[[pred_col]]),        # fixed
            linewidth = 1.5, show.legend = TRUE) +
    geom_sf(data = obs_sf,
            aes(color = .data[[obs_col]]),          # fixed
            shape = 21, fill = "white",
            size = 2, stroke = 2.5,
            inherit.aes = FALSE,
            show.legend = FALSE) +
    scale_color_gradientn(
      colors = pal,
      limits = c(scale_min, scale_max),
      oob = scales::squish,
      na.value = "grey80",
      name = "PFAS (ng/L)") +
    guides(color = guide_colorbar(
      barwidth  = grid::unit(0.6, "cm"),
      barheight = grid::unit(6, "cm"),
      label.theme = element_text(size = 12),
      title.theme = element_text(size = 13, face = "bold"),
      title.position = "top")) +
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
  
  message("Saved: ", compound)
}

p


#--- multipanel plot -----------------------------------

library(patchwork)
library(ggspatial)

wanted <- c("PFAS40", "PFSA", "PFOS_Results", "PFCA", "PFOA_Results")

title_map <- c(
  PFAS40        = "PFAS40",
  PFCA          = "PFCA",
  PFSA          = "PFSA",
  PFOS_Results  = "PFOS",
  PFOA_Results  = "PFOA",
  PFBS_Results  = "PFBS",
  PFHxS_Results = "PFHxS",
  PFNA_Results  = "PFNA"
)

plot_list <- list()
scale_max_multiplot <- 60

for (compound in wanted) {
  
  pred_col    <- paste0(compound, "_pred")
  obs_col     <- compound
  clean_title <- ifelse(compound %in% names(title_map),
                        title_map[[compound]],
                        compound)
  
  r2_val <- glance_df %>%
    filter(compound == !!compound) %>%
    { if (nrow(.) == 0) NA_real_ else
      if ("pseudo.r.squared" %in% names(.)) .$pseudo.r.squared else
        if ("r.squared" %in% names(.)) .$r.squared else NA_real_ }
  r2_label <- if (is.na(r2_val)) "R² = NA" else
    paste0("R\u00B2 = ", format(r2_val, digits = 3))
  
  plot_list[[compound]] <- ggplot() +
    geom_sf(data = catchment, fill = NA, color = "black", linewidth = 0.6) +
    geom_sf(data = PFAS_pred$edges,
            aes(color = .data[[pred_col]]),
            linewidth = 1.5, show.legend = TRUE) +
    geom_sf(data = obs_sf,
            aes(color = .data[[obs_col]]),
            shape = 21, fill = "white",
            size = 2, stroke = 2.5,
            inherit.aes = FALSE,
            show.legend = FALSE) +
    scale_color_gradientn(
      colors   = pal,
      limits   = c(scale_min, scale_max_multiplot),
      oob      = scales::squish,
      na.value = "grey80",
      name     = "PFAS (ng/L)",
      trans    = "log1p",
      breaks   = c(0, 5, 10, 20, 40, 60),
      labels   = c("0", "5", "10", "20", "40", "60")
    ) +
    guides(color = guide_colorbar(
      barwidth       = grid::unit(0.4, "cm"),
      barheight      = grid::unit(4, "cm"),
      label.theme    = element_text(size = 9),
      title.theme    = element_text(size = 10, face = "bold"),
      title.position = "top")) +
    labs(title = clean_title, subtitle = r2_label) +
    theme_classic() +
    theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,                 hjust = 0.5),
      axis.title    = element_blank(),
      axis.text     = element_text(size = 7),
      legend.title  = element_text(size = 10, face = "bold"),
      legend.text   = element_text(size = 9),
      legend.position = "right",
      plot.margin   = margin(t = 4, r = 4, b = 4, l = 4))
}

# Add north arrow and scale bar to PFAS40, remove axis text from all
for (compound in names(plot_list)) {
  
  if (compound == "PFAS40") {
    plot_list[[compound]] <- plot_list[[compound]] +
      annotation_north_arrow(
        location = "tl",
        style    = north_arrow_minimal(text_size = 8),
        height   = unit(0.8, "cm"),
        width    = unit(0.8, "cm")) +
      annotation_scale(
        location   = "bl",
        width_hint = 0.3,
        text_size  = 8)
  }
  
  plot_list[[compound]] <- plot_list[[compound]] +
    guides(color = guide_colorbar(
      barwidth       = grid::unit(0.4, "cm"),
      barheight      = grid::unit(8, "cm"),
      label.theme    = element_text(size = 9),
      title.theme    = element_text(size = 10, face = "bold"),
      title.position = "top")) +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
}

# Assemble with patchwork
final_plot <- plot_list[["PFAS40"]] +
  (plot_list[["PFSA"]]  / plot_list[["PFOS_Results"]]) +
  (plot_list[["PFCA"]]  / plot_list[["PFOA_Results"]]) +
  plot_layout(widths = c(2, 1, 1), guides = "collect") +
  plot_annotation(
    title    = "Predicted PFAS concentrations — Mill River Watershed",
    subtitle = "Spatial 3 (March 2026, n = 27) · Tailup exponential SSN model",
    caption  = "Circles = observed values · Stream color = predicted concentration",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      plot.caption  = element_text(size = 9,  hjust = 0.5)
    )
  )

final_plot

ggsave(
  filename = file.path(out_dir, "S3_multipanel_SSN_map.png"),
  plot     = final_plot,
  width    = 14, height = 8, dpi = 400
)


# 1. Check your residuals - which sites are poorly fit?
residuals_df <- augment(models[["PFAS40"]]) %>%
  st_drop_geometry() %>%
  select(pid, .fitted, .resid) %>%
  mutate(pid = as.integer(pid)) %>%           # convert to integer to match
  left_join(
    PFAS_ssn$obs %>% st_drop_geometry() %>% 
      select(pid, Name, PFAS40),
    by = "pid"
  ) %>%
  arrange(desc(abs(.resid)))

print(residuals_df, n = Inf)

# Visualize observed vs fitted
ggplot(residuals_df, aes(x = .fitted, y = PFAS40)) +
  geom_point(size = 3, color = "#185FA5") +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray50") +
  geom_text(aes(label = Name), 
            size = 2.8, hjust = -0.1, vjust = 0.5,
            check_overlap = TRUE) +
  labs(x = "Fitted (ng/L)", y = "Observed (ng/L)",
       title = "Observed vs fitted — PFAS40") +
  theme_classic()







#thats it for the SSN model!! now, the following are for supplementary information

#  using NLCD raster for land cover map
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(tidyterra)  
library(FedData)

# Download NLCD for your watershed extent
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

# Get catchments for your watershed COMIDs
comids <- PFAS_pred$edges %>% 
  st_drop_geometry() %>% 
  pull(comid)

catchments_sub <- get_nhdplus(
  comid = comids,
  realization = "catchment"
)

# Check it downloaded correctly
plot(st_geometry(catchments_sub))
ggplot() +
  geom_spatraster(data = nlcd_reclass, aes(fill = Class)) +
  scale_fill_manual(
    values = c(
      "1" = "#A8D8EA",
      "2" = "#AA0000",
      "3" = "#B2B2B2",
      "4" = "#38814E",
      "5" = "#CDB577",
      "6" = "#CA9146",
      "7" = "#3B9C8C"
    ),
    labels = c(
      "1" = "Open water",
      "2" = "Developed",
      "3" = "Barren",
      "4" = "Forest",
      "5" = "Shrub/grassland",
      "6" = "Agricultural",
      "7" = "Wetlands"
    ),
    na.translate = FALSE,
    name = "Land cover"
  ) +
  geom_sf(data = catchments_sub,        # subcatchment borders
          fill = NA, 
          color = "gray20", 
          linewidth = 0.2,
          linetype = "dashed") +
  geom_sf(data = flowline, color = "#A8D8EA", linewidth = 0.5) +
  geom_sf(data = catchment, fill = NA, color = "black", linewidth = 0.8) +
  coord_sf(crs = st_crs(catchment)) +
  theme_void() +
  theme(legend.position = "left",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 9))


#this makes maps of the streamcat metrics on the study area
# Get all COMIDs from your flowlines
all_comids <- PFAS_pred$edges %>%
  st_drop_geometry() %>%
  pull(comid) %>%
  unique() %>%
  na.omit()

# Download StreamCat for ALL flowline COMIDs
streamcat_edges <- sc_get_data(
  comid = all_comids,  # all flowline COMIDs from earlier
  metric = c("npdesdens", "pctimp2019", "fert",
             "huden2010", "rdcrs", "conn",
            "pctagdrainage", "septic"),
  aoi = "ws"
)

cat("Rows returned:", nrow(streamcat_edges), "\n")
cat("COMIDs requested:", length(all_comids), "\n")

# Rejoin to catchments
catchments_joined <- catchments_sub %>%
  mutate(featureid = as.integer(featureid)) %>%
  left_join(
    streamcat_edges %>%
      mutate(comid = as.integer(comid)) %>%
      rename(featureid = comid),
    by = "featureid"
  )

# Check coverage
missing <- sum(is.na(catchments_joined$fertws))
cat("Catchments still missing fertws:", missing, "\n")


# Define predictors and clean labels for supplementary
supp_predictors <- c("pctimp2019ws", "npdesdensws", "fertws", 
                     "huden2010ws", "rdcrsws", "pctagdrainagews", 
                     "connws", "septicws")

supp_labels <- c(
  pctimp2019ws = "% impervious",
  npdesdensws  = "NPDES density",
  fertws       = "Fertilizer application",
  huden2010ws  = "Human population density",
  rdcrsws      = "Road-stream crossings",
  pctagdrainagews = "% Agricultural drainage",
  connws       = "Watershed connectivity",
  septicws     = "Septic System Density"
)

# Just plot to screen first to check they look right
for (pred in supp_predictors) {
  
  p <- ggplot() +
    geom_sf(data = catchments_joined,
            aes(fill = .data[[pred]]),
            color = "gray40", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = supp_labels[[pred]],
      option = "magma",
      na.value = "gray90"
    ) +
    geom_sf(data = flowline, color = "#5DADE2", linewidth = 0.5) +
    geom_sf(data = catchment, fill = NA, color = "black", linewidth = 0.8) +
    labs(title = supp_labels[[pred]]) +
    coord_sf(crs = st_crs(catchment)) +
    theme_void() +
    theme(
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9)
    )
  
  print(p)  # explicitly print so it shows in loop
  Sys.sleep(1)  # pause 1 second between plots so you can see each one
}






#check the influence of each predictor on the 
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
want_comps <- c("PFCA", "PFSA", "PFAS40")

# renaming the prediction labels for easier readibility
pred_labels <- c(
  npdesdensws  = "NPDES density",
  pctimp2019ws = "% impervious",
  pctagdrainagews = "% Agricultural"
)

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
    text_color = ifelse(abs(std_coef) > 0.35, "white", "black")
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

plot_df

#this is trying to understand relevance of a few new areas, that were deemed hotspots in the first model

# Pull S3_14 and S3_15 compound profiles
obs_sf %>%
  st_drop_geometry() %>%
  filter(Name %in% c("UMass Farm", "Brandywine Pond")) %>%
  select(Name, PFAS40, PFOA_Results, PFOS_Results, PFBS_Results, 
         PFBA_Results, PFHxS_Results, PFCA, PFSA, FTCA) %>%
  print()

# Check FTCA and FTS compounds specifically - these are biosolid related
obs_sf %>%
  st_drop_geometry() %>%
  filter(Name %in% c("UMass Farm", "Brandywine Pond")) %>%
  select(Name, FTCA, FTSA, X6.2FTS_Results, X8.2FTS_Results, 
         X4.2FTS_Results, PFECA, PFESA) %>%
  print()


# Compare Knightly PFCA:PFSA ratio to all other sites
S1 %>%
  mutate(ratio_CA_SA = PFCA / PFSA) %>%
  select(Name, PFAS40, PFCA, PFSA, ratio_CA_SA) %>%
  arrange(desc(ratio_CA_SA)) %>%
  print(n = Inf)
