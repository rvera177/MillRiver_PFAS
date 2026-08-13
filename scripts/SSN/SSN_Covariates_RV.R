
#  title: "Mill River PFAS: Covariate development (generalized pipeline)"
#author: "adapted from Matt Fuller's original script"
#date: "8/13/26"
# output: html_document

  
  # Purpose
  # This code generalizes Matt's basin/RCA covariate workflow so new predictor rasters
  # (and re-derived layers can be added by editing
  # ONE registry list below, instead of hand-writing a new code block per raster.
  #
  # Two computation paths are preserved on purpose:
  #   - BASIN path: basins are nested/overlapping -> loop over basin_stack layers,
  #     mask() each covariate raster to each layer. Required because overlapping
  #     polygons cannot be represented as a single categorical zone raster.
  #   - RCA path: RCAs are mutually exclusive -> single categorical zone raster,
  #     use zonal() once across all reaches. Much faster; don't loop here.
  
library(whitebox); wbt_init()
library(sf)
library(terra)
library(tidyverse)
getwd()
list.files("../outputs/basins/", pattern = "basin_", full.names = TRUE)
setwd("C:/Users/Marston User/Documents/LWMR Isoscape/PFAS_MillRiver_SSN-main/PFAS_MillRiver_SSN-main/data")

# ---- 1. Read in basin stack, RCA zone raster, and pour points --------------

reach_basins <-
  list.files(path = "../outputs/basins/", pattern = "basin_", full.names = TRUE) |>
  lapply(terra::rast)

basin_stack <- rast(reach_basins)
template <- basin_stack[[1]]

reach_rca_r <- rast("../data/raster/final_pour_point_RCAs.tif")
cell_area_m2 <- prod(res(reach_rca_r))

final_pour_points <-
  st_read("../data/shp/final_pour_points_manual_corrections_20260518.shp") |>
  st_transform(crs(basin_stack))


# ---- 2. Alignment helper (crops in native CRS BEFORE reprojecting) ---------
# NOTE: this fixes the bug in the original AgTile block, which skipped
# project(crs(template)) while NLCD and NPDES did not.
#
# For CONUS-wide (or otherwise huge) source rasters, reprojecting the FULL
# extent first is very slow and unnecessary. Instead: reproject just the
# template's bounding box into the source raster's native CRS, crop to that
# (cheap - no resampling), THEN reproject the small clipped piece. A buffer
# is added before the native-CRS crop because reprojecting a rectangle can
# distort it into a non-rectangular shape - the buffer guards against
# clipping off real edge pixels before the final resample.

get_native_crop_extent <- function(template, target_crs, buffer_m = 1000) {
  template_poly <- as.polygons(ext(template), crs = crs(template))
  template_poly_proj <- project(template_poly, target_crs)
  ext(template_poly_proj) + buffer_m
}

align_to_template <- function(path, template, method = c("near", "bilinear"),
                              pre_clip_buffer_m = 1000, cache_dir = NULL) {
  method <- match.arg(method)
  
  # optional cache: skip re-cropping/reprojecting huge CONUS rasters on reruns
  if (!is.null(cache_dir)) {
    cache_path <- file.path(cache_dir, paste0(tools::file_path_sans_ext(basename(path)), "_aligned.tif"))
    if (file.exists(cache_path)) return(rast(cache_path))
  }
  
  r <- rast(path)
  
  if (!same.crs(r, template)) {
    native_ext <- get_native_crop_extent(template, crs(r), buffer_m = pre_clip_buffer_m)
    r <- crop(r, native_ext)          # cheap crop, still in native CRS/resolution
    r <- project(r, crs(template))    # reproject only the small clipped piece
  } else {
    r <- crop(r, ext(template) + pre_clip_buffer_m)
  }
  
  r_aligned <- r |>
    crop(template) |>
    extend(template) |>
    resample(template, method = method)
  
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    writeRaster(r_aligned, cache_path, overwrite = TRUE)
  }
  
  r_aligned
}



# ---- 3. Covariate registry --------------------------------------------------
# This is the only thing you edit to add a new predictor.
#   type:
#     "percent_cover" - binary/categorical raster (e.g. AgTile, wetlands mask);
#                        reports % of cells that are non-NA / "present"
#     "continuous"     - continuous raster (e.g. NLCD imperviousness, slope,
#                        precip); reports the mean value
#     "point_density"  - point vector layer (e.g. NPDES, dams); reports count
#                        and density per km2
#   resample_method: "near" for categorical/binary, "bilinear" for continuous
#
# For point layers, `path` should point to the vector file, not a raster.
# For raster layers, `path` should point to the GeoTIFF.
# reclassify_class: for binary land cover masks, specify the NLCD class(es) to keep as 1

covariate_registry <- list(
  agtile = list(
    path = "../data/raster/AgTile_clip_resample10m.tif",
    type = "percent_cover",
    resample_method = "near",
    out_name = "agtile"),
  #2025 impervious data source. https://www.mrlc.gov/data/type/fractional-impervious-surface
  nlcd_imp2025 = list( 
    path = "../data/raster/Annual_NLCD_FractionalImperviousSurface_2025_CU_C1V2/Annual_NLCD_FctImp_2025_CU_C1V2.tif",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "nlcd_imp2025"),
  nlcd_imp2019 = list(
    path = "../data/raster/NLCD2019_imp_clip_resample10m.tif",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "nlcd_imp2019"),
  # 2025 NLCD Cultivated Crops (class 81): Areas used for production of annual crops 
  # Data source: https://www.usgs.gov/centers/eros/science/annual-nlcd-land-cover-classification 
  # (corn, soybeans) and perennial woody crops (orchards). Crop vegetation > 20%.
  nlcd_cultivated_crops_2025 = list(
    path = "../data/raster/Annual_NLCD_LndCov_2025_CU_C1V2/Annual_NLCD_LndCov_2025_CU_C1V2.tif",
    type = "percent_cover",
    resample_method = "near",
    out_name = "nlcd_cultivated_crops_2025",
    reclassify_class = 81),
  # 2025 NLCD Pasture/Hay (class 82): Grasses, legumes, or grass/legume mixtures
  # planted for livestock grazing or hay/seed production on perennial cycle. 
  # Pasture/Hay vegetation > 20%.
  nlcd_pasture_hay_2025 = list(
    path = "../data/raster/Annual_NLCD_LndCov_2025_CU_C1V2/Annual_NLCD_LndCov_2025_CU_C1V2.tif",
    type = "percent_cover",
    resample_method = "near",
    out_name = "nlcd_pasture_hay_2025",
    reclassify_class = 82),
  # 2025 NLCD Planted/Cultivated (classes 81 + 82): Combined cultivated crops and 
  # pasture/hay. Both fall under Planted/Cultivated classification and are sources 
  # of non-point contamination from pesticides and biosolid applications.
  nlcd_planted_cultivated_2025 = list(
    path = "../data/raster/Annual_NLCD_LndCov_2025_CU_C1V2/Annual_NLCD_LndCov_2025_CU_C1V2.tif",
    type = "percent_cover",
    resample_method = "near",
    out_name = "nlcd_planted_cultivated_2025",
    reclassify_class = c(81, 82)),
  npdes = list(
    path = "../data/shp/USEPA_NPDES_pts.shp",
    type = "point_density",
    out_name = "npdes"),
  # Baseflow Index (BFI): 1km resolution, interpolated from USGS streamgages
  # Represents base flow component of streamflow from ground-water discharge
  # link to source: https://www.sciencebase.gov/catalog/item/631405c5d34e36012efa3192
  # actual raster data is in the pfi48grd.zip file
  bfi = list(
    path = "../data/raster/BaseFlowIndex_USGS_48grd/bfi48grd/w001001.adf",
    type = "continuous",
    resample_method = "bilinear",
    out_name = "bfi")
)
# ---- 4. Pre-align all registered rasters / load all point layers -----------

aligned_rasters <- list()
point_layers <- list()

for (nm in names(covariate_registry)) {
  cov <- covariate_registry[[nm]]
  if (cov$type == "point_density") {
    point_layers[[nm]] <- st_read(cov$path, quiet = TRUE) |>
      vect() |> project(crs(template))
  } else {
    r_aligned <- align_to_template(cov$path, template, cov$resample_method)
    
    # Apply reclassification if specified (for binary land cover masks)
    if (!is.null(cov$reclassify_class)) {
      r_aligned <- ifel(r_aligned %in% cov$reclassify_class, 1, NA)
    }
    
    aligned_rasters[[nm]] <- r_aligned
  }
}

# sanity check alignment
compareGeom(template, unlist(aligned_rasters)[[1]])


# ---- 5. Basin-side summary functions ----------------------------------------

summarize_percent_cover_basin <- function(basin_mask, data_raster) {
  data_masked <- mask(data_raster, basin_mask)
  n_cells <- sum(!is.na(values(basin_mask)))
  100 * sum(!is.na(values(data_masked))) / n_cells
}

summarize_continuous_mean_basin <- function(basin_mask, data_raster) {
  data_masked <- mask(data_raster, basin_mask)
  vals <- values(data_masked)
  mean(vals[!is.na(vals)])
}

summarize_point_density_basin <- function(basin_mask, pts) {
  n_cells <- sum(!is.na(values(basin_mask)))
  basin_area_km2 <- (n_cells * cell_area_m2) / 1e6
  
  pts_subset <- pts[ext(basin_mask), ]
  if (nrow(pts_subset) == 0) {
    n_pts <- 0
  } else {
    vals <- terra::extract(basin_mask, pts_subset)
    n_pts <- sum(!is.na(vals[, 2]))
  }
  list(n_points = n_pts, density_km2 = n_pts / basin_area_km2)
}


# ---- 6. Basin-side loop (generic over the registry) -------------------------

results_list <- lapply(seq_len(nlyr(basin_stack)), function(i) {
  
  basin_mask <- basin_stack[[i]]
  basin_id <- names(basin_mask)
  
  n_cells <- sum(!is.na(values(basin_mask)))
  basin_area_km2 <- (n_cells * cell_area_m2) / 1e6
  
  row <- data.frame(reach_id = basin_id, basin_area_km2 = basin_area_km2)
  
  for (nm in names(covariate_registry)) {
    cov <- covariate_registry[[nm]]
    
    if (cov$type == "percent_cover") {
      row[[paste0("pct_", cov$out_name, "_bas")]] <-
        summarize_percent_cover_basin(basin_mask, aligned_rasters[[nm]])
      
    } else if (cov$type == "continuous") {
      row[[paste0(cov$out_name, "_bas")]] <-
        summarize_continuous_mean_basin(basin_mask, aligned_rasters[[nm]])
      
    } else if (cov$type == "point_density") {
      pd <- summarize_point_density_basin(basin_mask, point_layers[[nm]])
      row[[paste0(cov$out_name, "_count_bas")]] <- pd$n_points
      row[[paste0(cov$out_name, "_density_km2_bas")]] <- pd$density_km2
    }
  }
  row
})

results <- do.call(rbind, results_list) |>
  mutate(reach_id = as.numeric(gsub("basin_", "", reach_id))) |>
  filter(!is.na(reach_id))
summary(results)

write_csv(results, "../data/RVsummary_basin_covariates.csv")
saveRDS(results, "../data/RVsummary_basin_covariates.rds")


# ---- 7. RCA-side: shared cell count -----------------------------------------

cell_count <- zonal(!is.na(reach_rca_r), reach_rca_r, fun = "sum") |>
  rename(reach_id = 1, cell_count = 2) |>
  mutate(rca_area_km2 = (cell_count * cell_area_m2) / 1e6)


# ---- 8. RCA-side loop (generic over the registry, uses fast zonal()) -------
# Point layers still need rasterizing to a count-per-cell raster before zonal()
# can sum them, matching the original npdes_r_aligned approach.

rca_results <- cell_count

for (nm in names(covariate_registry)) {
  cov <- covariate_registry[[nm]]
  
  if (cov$type == "continuous") {
    stat <- zonal(aligned_rasters[[nm]], reach_rca_r, fun = "mean", na.rm = TRUE) |>
      rename(reach_id = 1, value = 2) |>
      rename(!!paste0(cov$out_name, "_rca") := value)
    rca_results <- left_join(rca_results, stat, by = "reach_id")
    
  } else if (cov$type == "percent_cover") {
    covered_count <- zonal(aligned_rasters[[nm]], reach_rca_r, fun = "sum", na.rm = TRUE) |>
      rename(reach_id = 1, covered_count = 2) |>
      mutate(covered_count = ifelse(is.na(covered_count), 0, covered_count))
    stat <- left_join(covered_count, cell_count, by = "reach_id") |>
      mutate(!!paste0("pct_", cov$out_name, "_rca") := (covered_count / cell_count) * 100) |>
      select(reach_id, !!paste0("pct_", cov$out_name, "_rca"))
    rca_results <- left_join(rca_results, stat, by = "reach_id")
    
  } else if (cov$type == "point_density") {
    # rasterize points to a countable raster on the template grid first
    pts_r <- rasterize(point_layers[[nm]], template, fun = "length", background = 0)
    pt_count <- zonal(pts_r, reach_rca_r, fun = "sum", na.rm = TRUE) |>
      rename(reach_id = 1, n_pts = 2) |>
      mutate(n_pts = ifelse(is.na(n_pts), 0, n_pts))
    stat <- left_join(pt_count, cell_count, by = "reach_id") |>
      mutate(
        !!paste0(cov$out_name, "_count_rca") := n_pts,
        !!paste0(cov$out_name, "_density_km2_rca") := (n_pts / (cell_count * cell_area_m2)) * 1e6
      ) |>
      select(reach_id, contains(cov$out_name))
    rca_results <- left_join(rca_results, stat, by = "reach_id")
  }
}

summary(rca_results)

write_csv(rca_results, "../data/RVsummary_RCA_covariates.csv")
saveRDS(rca_results, "../data/RVsummary_RCA_covariates.rds")


# ---- 9. Visualize -----------------------------------------------------------

results |>
  pivot_longer(cols = -reach_id, names_to = "covariate", values_to = "value") |>
  ggplot(aes(x = covariate, y = value)) +
  geom_violin() +
  scale_y_log10() +
  facet_wrap(~covariate, scales = "free", nrow = 1)

rca_results |>
  pivot_longer(cols = -reach_id, names_to = "covariate", values_to = "value") |>
  ggplot(aes(x = covariate, y = value)) +
  geom_violin() +
  scale_y_log10() +
  facet_wrap(~covariate, scales = "free", nrow = 1)


# ---- 10. Visualization: Spatial distribution of predictors ----------------
library(tmap)
library(RColorBrewer)

# Convert reach/basin boundaries to spatial for mapping
final_pour_points_sf <- final_pour_points

reach_rca_sf <- as.polygons(reach_rca_r) |> st_as_sf()

basin_boundaries <- lapply(reach_basins, as.polygons) |> 
  lapply(st_as_sf) |> 
  bind_rows()

# Plot 1: Raw BFI raster with reach boundaries
tm_shape(aligned_rasters[["bfi"]]) + 
  tm_raster(title = "Base Flow Index", palette = "Blues") +
  tm_shape(reach_rca_sf) + 
  tm_borders(col = "red", lwd = 2) +
  tm_layout(title = "BFI Distribution across Reaches", frame = TRUE)

# Plot 2: Raw NLCD Planted/Cultivated with reach boundaries
tm_shape(aligned_rasters[["nlcd_planted_cultivated_2025"]]) + 
  tm_raster(title = "Planted/Cultivated (1=Yes, NA=No)", palette = "Greens") +
  tm_shape(reach_rca_sf) + 
  tm_borders(col = "red", lwd = 2) +
  tm_layout(title = "Planted/Cultivated Land across Reaches", frame = TRUE)

# Rename the raster value column to reach_id for joining
reach_rca_sf <- reach_rca_sf |>
  rename(reach_id = final_pour_point_RCAs)


rca_results_sf <- left_join(reach_rca_sf, rca_results, by = "reach_id")

# Plot 3: Choropleth of summarized BFI by reach (blue)
tm_shape(rca_results_sf) + 
  tm_fill("bfi_rca", col.legend = tm_legend(title = "Mean BFI"), col.scale = tm_scale_continuous(values = "brewer.blues")) +
  tm_borders() +
  tm_title("BFI: Reach-scale Summary")

# Plot 4: Choropleth of % Planted/Cultivated by reach (orange for farmland)
tm_shape(rca_results_sf) + 
  tm_fill("pct_nlcd_planted_cultivated_2025_rca", 
          palette = "Oranges",
          title = "% Planted/Cultivated") +
  tm_borders() +
  tm_title("Planted/Cultivated: Reach-scale Summary (%)")

# Plot 5: Choropleth of % NLCD Impervious 2025 by reach
tm_shape(rca_results_sf) + 
  tm_fill("nlcd_imp2025_rca", 
          palette = "Purples",
          title = "% Impervious") +
  tm_borders() +
  tm_title("NLCD Impervious 2025: Reach-scale Summary (%)")


basin_summary_sf <- lapply(seq_len(nlyr(basin_stack)), function(i) {
  basin_sf <- as.polygons(basin_stack[[i]]) |> st_as_sf()
  basin_sf <- basin_sf |> 
    select(geometry) |>  # keep only geometry, drop inconsistent raster value columns
    mutate(reach_id = results$reach_id[i])
  basin_sf
}) |> 
  bind_rows()

# Now join with results
basin_summary_sf <- left_join(basin_summary_sf, results, by = "reach_id")

# Plot: Choropleth for basin-scale planted/cultivated
tm_shape(basin_summary_sf) + 
  tm_fill("pct_nlcd_planted_cultivated_2025_bas", 
          palette = "Oranges",
          title = "% Planted/Cultivated") +
  tm_borders() +
  tm_title("Planted/Cultivated: Basin-scale Summary (%)")
