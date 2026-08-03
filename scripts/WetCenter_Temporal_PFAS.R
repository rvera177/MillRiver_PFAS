

# MILL RIVER PFAS ANALYSIS - COMPLETE WORKFLOW v3
# Reorganized: Flow Records First, Then PFAS Processing

# --- Load All Libraries Once ---
library(tidyverse)
library(lubridate)
library(readr)
library(ggplot2)
library(data.table)
library(highcharter)
library(readxl)

# --- CONFIGURATION (single source of truth) ---
CONFIG <- list(
  TIMEZONE = "America/New_York",
  DATE_FORMAT = "%m/%d/%y %H:%M",
  DATE_FORMAT_FULL = "%m/%d/%Y %H:%M",
  PAD_HOURS = 24,
  OVERLAP_FILTER_VALUE = 14.9,
  COR_MISSING_THRESH = 0.85,
  STAGE_SEQ_RANGE = c(0.20, 0.65),
  STAGE_SEQ_BY = 0.01
)

# --- DATA URLs (centralized for easy updates) ---
DATA_URLS <- list(
  water_pressure = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation4_WaterPressure.csv",
  air_pressure = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation3_AirPressure.csv",
  gauge_height = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenter%20GaugeHeight%20-%20Sheet1-2.csv",
  pfas_storm16 = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/RV_PFAS_Storm16.csv",
  precipitation = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PrecipitationData_Storms16.csv",
  pfas_inventory = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PFAS%20Inventory%20Downloaded%20July%2021%202026.csv",
  mdl_file = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/EPA1633_MDL_UMass_Amherst_Engineering.csv",
  wet_center_pfas = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenterPFASResults.csv",
  wet_center_discharge = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/2025_1YEAR_WET%20Center%20Umass%20AmherstWETCENTERUMASSMILL%20RIVER%20LEVEL%20SENSOR1_yearchannel1.csv",
  storm_event_15 = "https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Storm_Event_15_all_results.csv"
)

# --- Define PFAS compounds to analyze ---
PFAS_COMPOUNDS <- c(
  "11Cl-PF3OUdS", "3-3 FTCA", "4:2FTS", "5-3 FTCA", "6:2FTS", "7-3 FTCA", "8:2FTS",
  "9Cl-PF3ONS", "ADONA", "HFPO-DA", "NEtFOSA", "NEtFOSAA", "NEtFOSE", "NFDHA",
  "NMeFOSA", "NMeFOSAA", "NMeFOSE", "PFBA", "PFBS", "PFDA", "PFDoA", "PFDoS",
  "PFDS", "PFEESA", "PFHpA", "PFHpS", "PFHxA", "PFHxS", "PFMBA", "PFMPA",
  "PFNA", "PFNS", "PFOA", "PFOS", "PFOSA", "PFPeA", "PFPeS", "PFTeDA", "PFTrDA", "PFUnA"
)

# --- HELPER FUNCTIONS ---

# Safe CSV loading with error handling
load_csv_safely <- function(url, name = NULL) {
  tryCatch(
    read_csv(url, show_col_types = FALSE),
    error = function(e) {
      warning(paste("Failed to load:", name %||% url))
      NULL
    }
  )
}

# Load Excel file safely
load_excel_safely <- function(url, sheet = 1, name = NULL) {
  tryCatch(
    read_excel(url, sheet = sheet),
    error = function(e) {
      warning(paste("Failed to load Excel:", name %||% url))
      NULL
    }
  )
}

# Parse DateTime with flexible formats
parse_mill_datetime <- function(date_col, time_col = NULL, tz = CONFIG$TIMEZONE) {
  if (!is.null(time_col)) {
    paste_str <- paste(date_col, format(time_col, "%H:%M"))
    fmt <- CONFIG$DATE_FORMAT
  } else {
    paste_str <- date_col
    fmt <- CONFIG$DATE_FORMAT_FULL
  }
  as.POSIXct(paste_str, format = fmt, tz = tz)
}

# Compute flow from gauge height using power law
compute_flow <- function(height_m, a, b) {
  a * (height_m ^ b)
}

# Rolling nearest join wrapper
roll_join_nearest <- function(dt1, dt2, on_col = "DateTime") {
  setDT(dt1); setDT(dt2)
  dt2[dt1, on = on_col, roll = "nearest"]
}

# Standardized ggplot theme
theme_mill <- function() {
  theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 12)
    )
}

# Create event dataframe with parsed dates
create_events <- function(event_specs) {
  event_specs %>%
    mutate(
      start = as.POSIXct(start, tz = CONFIG$TIMEZONE),
      end = as.POSIXct(end, tz = CONFIG$TIMEZONE)
    )
}

# Define event color palette
event_colors <- c(
  "Storm" = "steelblue",
  "Baseflow" = "forestgreen",
  "Snowmelt" = "mediumpurple"
)

# ============================================================
# SECTION 1: PRESSURE SENSOR CALIBRATION
# ============================================================

cat("\n========== SECTION 1: PRESSURE SENSOR CALIBRATION ==========\n")

# --- Load raw sensor data ---
raw_data <- map(DATA_URLS[c("water_pressure", "air_pressure", "gauge_height")],
                ~load_csv_safely(.x))

# --- Parse sensor timestamps and standardize columns ---
water_ts <- raw_data$water_pressure %>%
  mutate(DateTime = parse_date_time(`Date-Time (EDT)`, 
                                    orders = c("mdy HM", "mdy HMS"),
                                    tz = CONFIG$TIMEZONE)) %>%
  rename(Water_Abs_kPa = `Absolute Pressure kPa`) %>%
  select(DateTime, Water_Abs_kPa)

air_ts <- raw_data$air_pressure %>%
  mutate(DateTime = parse_date_time(`Date-Time (EDT)`, 
                                    orders = c("mdy HM", "mdy HMS"),
                                    tz = CONFIG$TIMEZONE)) %>%
  rename(Air_Abs_kPa = `Absolute Pressure kPa`) %>%
  select(DateTime, Air_Abs_kPa)

# --- Barometric compensation via rolling nearest join ---
sensor_ts <- roll_join_nearest(water_ts, air_ts, "DateTime") %>%
  as_tibble() %>%
  mutate(Water_Pressure_Only = Water_Abs_kPa - Air_Abs_kPa)

# --- Parse manual gauge readings ---
gauge_clean <- raw_data$gauge_height %>%
  filter(!is.na(`GaugeHeight (ft)`)) %>%
  mutate(
    DateTime = parse_mill_datetime(Date, Time),
    GaugeHeight_ft = `GaugeHeight (ft)`
  ) %>%
  select(DateTime, GaugeHeight_ft)

# --- Match gauge readings to nearest sensor record ---
calib_points <- roll_join_nearest(gauge_clean, sensor_ts, "DateTime") %>%
  as_tibble() %>%
  filter(!is.na(Water_Pressure_Only), !is.na(GaugeHeight_ft))

cat("Calibration points matched:", nrow(calib_points), "\n")

# --- Fit pressure-to-height calibration model ---
calib_model <- lm(GaugeHeight_ft ~ Water_Pressure_Only, data = calib_points)
cat("\nCalibration Model Summary:\n")
print(summary(calib_model))

# --- Apply calibration to full time series ---
sensor_calibrated <- sensor_ts %>%
  mutate(
    Calibrated_Height_ft = predict(calib_model,
                                   newdata = data.frame(Water_Pressure_Only = Water_Pressure_Only))
  )

# --- Plot: Calibration verification ---
print(ggplot() +
        geom_line(data = sensor_calibrated,
                  aes(x = DateTime, y = Calibrated_Height_ft, color = "Sensor"),
                  linewidth = 0.5, alpha = 0.8) +
        geom_point(data = calib_points,
                   aes(x = DateTime, y = GaugeHeight_ft, color = "Manual Gauge"),
                   size = 3, shape = 21, fill = "red") +
        scale_color_manual(values = c("Sensor" = "steelblue", "Manual Gauge" = "red")) +
        labs(title = "Sensor Calibration: Pressure-to-Height",
             x = "Date/Time (EDT)", y = "Water Level (ft)") +
        theme_mill())

# ============================================================
# SECTION 2: RATING CURVE FITTING
# ============================================================

cat("\n========== SECTION 2: RATING CURVE FITTING ==========\n")

# --- Create rating curve dataset ---
rating_data <- tribble(
  ~Stage_m, ~Flow_m3s,
  0.23, 0.34,
  0.27, 0.42,
  0.58, 1.125,
  0.58, 1.146,
  0.26, 0.179,
  0.26, 0.176,
  0.37, 0.522,
  0.36, 0.505,
  0.29, 0.302,
  0.29, 0.299,
  0.49, 0.755,
  0.48, 0.724
) %>%
  group_by(Stage_m) %>%
  summarise(Flow_m3s = mean(Flow_m3s), .groups = "drop")

# --- Fit power law via NLS ---
rating_model <- nls(Flow_m3s ~ a * Stage_m^b,
                    data = rating_data,
                    start = list(a = 2.74, b = 1.73))

# Store rating curve parameters globally
RATING_PARAMS <- list(
  a = coef(rating_model)["a"],
  b = coef(rating_model)["b"]
)

cat("\n\nRating Curve Model:\n")
cat("Q =", round(RATING_PARAMS$a, 3), "* H^", round(RATING_PARAMS$b, 3), "\n")
print(summary(rating_model))

# --- Apply rating curve to calibrated sensor data ---
WetStation4Stream <- sensor_calibrated %>%
  rename(Time = DateTime) %>%
  mutate(
    GaugeHeight_m = Calibrated_Height_ft * 0.3048,  # feet to meters
    Flow = compute_flow(GaugeHeight_m, RATING_PARAMS$a, RATING_PARAMS$b)
  ) %>%
  select(Time, GaugeHeight_m, Calibrated_Height_ft, Flow) %>%
  arrange(Time)

# --- Deduplicate by averaging (if any exact duplicates) ---
WetStation4Stream <- WetStation4Stream %>%
  group_by(Time) %>%
  summarise(
    GaugeHeight_m = mean(GaugeHeight_m, na.rm = TRUE),
    Calibrated_Height_ft = mean(Calibrated_Height_ft, na.rm = TRUE),
    Flow = mean(Flow, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Time)

cat("\n\nWetStation4 Time Series:\n")
cat("Date range:", format(min(WetStation4Stream$Time, na.rm = TRUE), "%Y-%m-%d"),
    "to", format(max(WetStation4Stream$Time, na.rm = TRUE), "%Y-%m-%d"), "\n")
cat("Flow range (m³/s):", round(range(WetStation4Stream$Flow, na.rm = TRUE), 3), "\n")

# --- Plot: Rating curve with data ---
stage_seq <- seq(CONFIG$STAGE_SEQ_RANGE[1], CONFIG$STAGE_SEQ_RANGE[2], 
                 by = CONFIG$STAGE_SEQ_BY)
rating_predicted <- tibble(
  Stage_m = stage_seq,
  Flow_m3s = compute_flow(stage_seq, RATING_PARAMS$a, RATING_PARAMS$b)
)

print(ggplot() +
        geom_point(data = rating_data, aes(x = Stage_m, y = Flow_m3s),
                   color = "red", size = 3, alpha = 0.7) +
        geom_line(data = rating_predicted, aes(x = Stage_m, y = Flow_m3s),
                  color = "steelblue", linewidth = 1) +
        labs(title = "Rating Curve: Stage vs Discharge",
             x = "Stage (m)", y = "Discharge (m³/s)") +
        theme_mill())

# ============================================================
# SECTION 3: TIMEVIEW DISCHARGE RECONCILIATION
# ============================================================

cat("\n========== SECTION 3: TIMEVIEW DISCHARGE RECONCILIATION ==========\n")

# --- Load TimeView discharge data ---
TimeViewDischarge_raw <- load_csv_safely(DATA_URLS$wet_center_discharge, 
                                         "TimeView Discharge")

WetCenterDischarge <- TimeViewDischarge_raw %>%
  mutate(
    Time = mdy_hm(DateTime, tz = CONFIG$TIMEZONE),
    GaugeHeight_m = `Stage (m)`,
    Flow = compute_flow(GaugeHeight_m, RATING_PARAMS$a, RATING_PARAMS$b)
  ) %>%
  filter(!is.na(Time), !is.na(GaugeHeight_m)) %>%
  select(Time, GaugeHeight_m, Flow)

cat("\n\nWetCenter TimeView Summary:\n")
cat("Date range:", format(min(WetCenterDischarge$Time, na.rm = TRUE), "%Y-%m-%d"),
    "to", format(max(WetCenterDischarge$Time, na.rm = TRUE), "%Y-%m-%d"), "\n")
cat("Stage range (m):", round(range(WetCenterDischarge$GaugeHeight_m, na.rm = TRUE), 3), "\n")

# --- Define overlap period and filter outliers ---
overlap_start <- max(min(WetCenterDischarge$Time, na.rm = TRUE),
                     min(WetStation4Stream$Time, na.rm = TRUE))
overlap_end <- min(max(WetCenterDischarge$Time, na.rm = TRUE),
                   max(WetStation4Stream$Time, na.rm = TRUE))

cat("Overlap period:", format(overlap_start, "%Y-%m-%d"), "to", 
    format(overlap_end, "%Y-%m-%d"), "\n")

WetCenterDischarge <- WetCenterDischarge %>%
  filter(GaugeHeight_m < CONFIG$OVERLAP_FILTER_VALUE)

# --- Perform rolling nearest join for overlap matching ---
pt_dt <- WetStation4Stream %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  select(Time, GaugeHeight_m) %>%
  rename(H_pt = GaugeHeight_m) %>%
  as.data.table()

tv_dt <- WetCenterDischarge %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  select(Time, GaugeHeight_m) %>%
  as.data.table()

overlap_matched <- tv_dt[pt_dt, on = "Time", roll = "nearest"] %>%
  as_tibble() %>%
  rename(H_tv = GaugeHeight_m) %>%
  filter(!is.na(H_tv), !is.na(H_pt))

# --- Fit correction model ---
overlap_model <- lm(H_pt ~ H_tv, data = overlap_matched)
intercept <- coef(overlap_model)["(Intercept)"]
slope <- coef(overlap_model)["H_tv"]

cat("\n\nOverlap Correction Model:\n")
cat("H_corrected =", round(slope, 6), "* H_tv +", round(intercept, 4), "\n")
cat("Correlation:", round(cor(overlap_matched$H_tv, overlap_matched$H_pt), 3), "\n")
print(summary(overlap_model))

# --- Apply correction ---
WetCenterDischarge <- WetCenterDischarge %>%
  mutate(
    GaugeHeight_m = slope * GaugeHeight_m + intercept,
    Flow = compute_flow(GaugeHeight_m, RATING_PARAMS$a, RATING_PARAMS$b)
  )

# --- Plot: Overlap comparison before/after ---
print(ggplot(overlap_matched, aes(x = H_tv, y = H_pt)) +
        geom_point(alpha = 0.3, color = "steelblue") +
        geom_smooth(method = "lm", color = "red", se = TRUE) +
        labs(title = "TimeView vs Pressure Transducer — Overlap Period",
             x = "TimeView Stage (raw)", y = "Pressure Transducer Stage (m)") +
        theme_mill())

# --- Plot: Time series comparison corrected ---
print(ggplot() +
        geom_line(data = WetCenterDischarge %>%
                    filter(Time >= overlap_start, Time <= overlap_end),
                  aes(x = Time, y = GaugeHeight_m, color = "TimeView"),
                  linewidth = 0.8) +
        geom_line(data = WetStation4Stream %>%
                    filter(Time >= overlap_start, Time <= overlap_end),
                  aes(x = Time, y = GaugeHeight_m, color = "Pressure Transducer"),
                  linewidth = 0.8) +
        scale_color_manual(values = c("TimeView" = "darkorange",
                                      "Pressure Transducer" = "steelblue")) +
        labs(title = "Overlap Period — Corrected TimeView vs Pressure Transducer",
             x = "Date/Time", y = "Stage (m)") +
        theme_mill())

# ============================================================
# SECTION 4: COMBINED FLOW RECORD & VISUALIZATIONS
# ============================================================

cat("\n========== SECTION 4: COMBINED FLOW RECORD & VISUALIZATIONS ==========\n")

# --- Combine both datasets (avoid duplication) ---
WetCenterDischarge_trimmed <- WetCenterDischarge %>%
  filter(Time < overlap_start) %>%
  transmute(Time, GaugeHeight_m, Flow, Source = "TimeView")

WetStation4Stream_labeled <- WetStation4Stream %>%
  slice(-1) %>%
  transmute(Time, GaugeHeight_m, Flow, Source = "Pressure Transducer")

WetCenter_Flow_Combined <- bind_rows(WetCenterDischarge_trimmed,
                                     WetStation4Stream_labeled) %>%
  arrange(Time) %>%
  mutate(Time = as.POSIXct(Time, tz = CONFIG$TIMEZONE))

cat("\n\nCombined Flow Record:\n")
cat("Date range:", format(min(WetCenter_Flow_Combined$Time, na.rm = TRUE), "%Y-%m-%d"),
    "to", format(max(WetCenter_Flow_Combined$Time, na.rm = TRUE), "%Y-%m-%d"), "\n")
cat("TimeView records:", sum(WetCenter_Flow_Combined$Source == "TimeView"), "\n")
cat("Pressure Transducer records:", sum(WetCenter_Flow_Combined$Source == "Pressure Transducer"), "\n")

# --- Plot: Full combined flow record ---
print(ggplot(WetCenter_Flow_Combined, aes(x = Time, y = Flow, color = Source)) +
        geom_line(linewidth = 0.6) +
        scale_color_manual(values = c("TimeView" = "darkorange",
                                      "Pressure Transducer" = "steelblue")) +
        labs(title = "Mill River — Combined Flow Record",
             x = "Date/Time", y = "Discharge (m³/s)") +
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month", expand = c(0, 0)) +
        theme_mill())

# --- Define event periods ---
events <- create_events(tribble(
  ~label, ~start, ~end, ~type,
  "Storm 15", "2025-10-10 00:00", "2025-10-16 00:00", "Storm",
  "Storm 16", "2025-10-29 20:00", "2025-11-07 18:00", "Storm",
  "12", "2025-07-09 00:00", "2025-07-12 10:20", "Storm",
  "13", "2025-07-14 05:00", "2025-07-16 13:00", "Storm"
))

# --- Plot: Combined flow with event windows ---
print(ggplot() +
        geom_rect(data = events,
                  aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = type),
                  alpha = 0.25) +
        geom_line(data = WetCenter_Flow_Combined,
                  aes(x = Time, y = Flow, color = Source),
                  linewidth = 0.6) +
        geom_text(data = events,
                  aes(x = start + (end - start) / 2, y = Inf, label = label),
                  vjust = 1.5, size = 3.5, fontface = "bold") +
        scale_fill_manual(values = event_colors, name = "Event Type") +
        scale_color_manual(values = c("TimeView" = "darkorange",
                                      "Pressure Transducer" = "steelblue"),
                           name = "Source") +
        labs(title = "Mill River — Combined Flow Record with Sampling Events",
             x = "Date/Time", y = "Discharge (m³/s)") +
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month", expand = c(0, 0)) +
        theme_mill())

# SECTION 5: LOAD & PARSE PFAS INVENTORY (CORRECTED)

# ============================================================

cat("\n========== SECTION 5: LOAD & PARSE PFAS INVENTORY ==========\n")

# --- Load MDL (Method Detection Limit) data ---

MDL_TABLE <- tibble(
  
  Compound = PFAS_COMPOUNDS,
  
  MDL_ng_L = c(0.693028, 0.302376, 0.787358, 4.041153, 1.490055, 4.535903, 0.77321,
               
               0.692001, 0.765486, 0.789132, 0.497405, 0.288365, 1.459179, 0.348747,
               
               0.258046, 0.453422, 1.248901, 0.655795, 0.330108, 0.261853, 0.184027,
               
               0.40008, 0.430525, 0.303981, 0.300411, 0.365672, 0.189522, 0.56971,
               
               0.379403, 0.266126, 0.247594, 0.353461, 0.309123, 0.382599, 0.465255,
               
               0.344344, 0.247385, 0.175991, 0.267794, 0.248134),
  
  MRL_ng_L = c(3.94, 3.786, 3.202, 24.757, 5.173, 23.521, 4.143,
               
               4.278, 3.962, 3.927, 0.961, 0.941, 10.284, 2.245,
               
               1.215, 1.152, 10.649, 4.232, 1.042, 0.901, 1.011,
               
               0.897, 0.934, 2.021, 1.028, 1.024, 1.013, 1.204,
               
               1.9, 1.754, 1.026, 0.859, 1.049, 1.097, 0.645,
               
               2.01, 1.073, 1.032, 1.003, 0.891)
  
)

cat("Loaded MDL table with", nrow(MDL_TABLE), "PFAS compounds\n")

# --- Fixed DateTime Parser for <time> class ---

# The Time column from readr comes as <time> class, which needs special handling

parse_mill_datetime_fixed <- function(date_col, time_col = NULL, tz = CONFIG$TIMEZONE) {
  
  if (!is.null(time_col)) {
    
    # Handle <time> or <difftime> class by converting to character HH:MM format
    
    if (inherits(time_col, "hms") || inherits(time_col, "difftime") || class(time_col)[1] == "time") {
      
      time_str <- format(time_col, "%H:%M")
      
    } else {
      
      time_str <- format(time_col, "%H:%M")
      
    }
    
    paste_str <- paste(date_col, time_str)
    
    fmt <- "%m/%d/%Y %H:%M"
    
  } else {
    
    paste_str <- date_col
    
    fmt <- CONFIG$DATE_FORMAT_FULL
    
  }
  
  as.POSIXct(paste_str, format = fmt, tz = tz)
  
}

# --- Load comprehensive PFAS inventory ---

PFAS_inventory_raw <- load_csv_safely(DATA_URLS$pfas_inventory, "PFAS Inventory")

cat("PFAS Inventory loaded with", nrow(PFAS_inventory_raw), "samples\n\n")

# --- Parse datetime and clean inventory ---

PFAS_inventory_full <- PFAS_inventory_raw %>%
  
  mutate(
    
    DateTime = parse_mill_datetime_fixed(Date, Time), # Use fixed parser
    
    Sample_ID = `Sample Name on LC/MS` %||% paste0("S_", row_number()),
    
    PFAS_Collected_Flag = `PFAS Collected?` %in% c("Yes", "yes", "YES"),
    
    PFAS_Analyzed_Flag = `PFAS Analyzed?` %in% c("Yes", "yes", "YES")
    
  ) %>%
  
  filter(Type == "Stream", PFAS_Collected_Flag)

cat("After filtering to stream samples with PFAS collected:", nrow(PFAS_inventory_full), "samples\n")

# --- Check DateTime parsing ---

cat("\nDateTime parsing check:\n")

cat(" • Samples with valid DateTime:", sum(!is.na(PFAS_inventory_full$DateTime)), "\n")

cat(" • Samples with NA DateTime:", sum(is.na(PFAS_inventory_full$DateTime)), "\n")

cat(" • Date range:", format(min(PFAS_inventory_full$DateTime, na.rm = TRUE), "%Y-%m-%d"),
    
    "to", format(max(PFAS_inventory_full$DateTime, na.rm = TRUE), "%Y-%m-%d"), "\n\n")

# --- Extract PFAS concentration data and reshape to long format ---

PFAS_long <- PFAS_inventory_full %>%
  
  filter(!is.na(DateTime)) %>% # Remove rows with unparseable dates
  
  select(Sample_ID, DateTime, Location, Project, contains("Results")) %>%
  
  pivot_longer(
    
    cols = contains("Results"),
    
    names_to = "Compound_Result",
    
    values_to = "Concentration_Str"
    
  ) %>%
  
  mutate(
    
    Compound = str_replace(Compound_Result, " Results", ""),
    
    Concentration_ng_L = as.numeric(str_extract(Concentration_Str, "^[0-9.]+")),
    
    Concentration_Str = as.character(Concentration_Str)
    
  ) %>%
  
  filter(!is.na(Concentration_ng_L)) %>%
  
  select(Sample_ID, DateTime, Location, Project, Compound, Concentration_ng_L, Concentration_Str) %>%
  
  left_join(MDL_TABLE %>% select(Compound, MDL_ng_L, MRL_ng_L), by = "Compound") %>%
  
  mutate(
    
    BDL_Flag = str_detect(Concentration_Str, "bdl"),
    
    Detected = if_else(BDL_Flag, FALSE, Concentration_ng_L >= MDL_ng_L),
    
    Detection_Status = case_when(
      
      BDL_Flag ~ "BDL",
      
      Concentration_ng_L < MDL_ng_L ~ "Below MDL",
      
      Concentration_ng_L >= MRL_ng_L ~ "Above MRL",
      
      TRUE ~ "Between MDL and MRL"
      
    )
    
  ) %>%
  
  arrange(DateTime, Compound)

cat("\n\nPFAS Long Format Summary:\n")

cat("Total compound measurements:", nrow(PFAS_long), "\n")

cat("Unique samples:", n_distinct(PFAS_long$Sample_ID), "\n")

cat("Unique compounds:", n_distinct(PFAS_long$Compound), "\n")

cat("Unique locations:", n_distinct(PFAS_long$Location), "\n")

cat("\nDetection Summary:\n")

print(table(PFAS_long$Detection_Status))

cat("\n\nUnique Locations:\n")

print(unique(PFAS_long$Location))

# --- Create summary by sample (all compounds detected) ---

PFAS_by_sample <- PFAS_long %>%
  
  group_by(Sample_ID, DateTime, Location, Project) %>%
  
  summarise(
    
    N_Compounds_Measured = n(),
    
    N_Detected = sum(Detected),
    
    N_BDL = sum(BDL_Flag),
    
    Sum_Detected_Conc_ng_L = sum(Concentration_ng_L[Detected], na.rm = TRUE),
    
    Max_Single_Compound_ng_L = max(Concentration_ng_L, na.rm = TRUE),
    
    .groups = "drop"
    
  ) %>%
  
  mutate(DateTime = as.POSIXct(DateTime, tz = CONFIG$TIMEZONE)) %>%
  
  arrange(DateTime)

cat("\n\nPFAS Summary by Sample:\n")

cat("Total samples:", nrow(PFAS_by_sample), "\n")

cat("Date range:", format(min(PFAS_by_sample$DateTime, na.rm = TRUE), "%Y-%m-%d"),
    
    "to", format(max(PFAS_by_sample$DateTime, na.rm = TRUE), "%Y-%m-%d"), "\n\n")

print(PFAS_by_sample %>% select(DateTime, Location, Project, N_Detected, Sum_Detected_Conc_ng_L))

# ============================================================
# SECTION 6: KEY PFAS ANALYSIS
# ============================================================

cat("\n========== SECTION 6: KEY PFAS ANALYSIS ==========\n")

# --- Extract key compounds for detailed analysis ---
KEY_COMPOUNDS <- c("PFOA", "PFOS", "PFHxA", "PFHxS", "PFDA", "PFNA", "PFBA")

PFAS_key <- PFAS_long %>%
  filter(Compound %in% KEY_COMPOUNDS) %>%
  select(Sample_ID, DateTime, Location, Compound, Concentration_ng_L, 
         Detected, Detection_Status, MDL_ng_L)

cat("\n\nKey PFAS Compounds Analysis:\n")
cat("Samples with detected key compounds:\n")
print(
  PFAS_key %>%
    filter(Detected) %>%
    group_by(Compound) %>%
    summarise(
      N_Detections = n(),
      Mean_Conc = mean(Concentration_ng_L, na.rm = TRUE),
      Max_Conc = max(Concentration_ng_L, na.rm = TRUE),
      Min_Conc = min(Concentration_ng_L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(N_Detections))
)

print(ggplot(PFAS_key %>% filter(Detected),
             
             aes(x = Compound, y = Concentration_ng_L, fill = Location)) +
        
        geom_boxplot(alpha = 0.7) +
        
        scale_y_log10() +
        
        labs(title = "Detected Key PFAS Compounds by Location",
             
             y = "Concentration (ng/L, log scale)", x = "Compound") +
        
        theme_mill() +
        
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# ----SECTION 7: PFAS STATISTICS & SUMMARIES----------

# --- Summary statistics by location ---

location_summary <- PFAS_by_sample %>%
  
  group_by(Location) %>%
  
  summarise(
    
    N_Samples = n(),
    
    Avg_Total_PFAS = mean(Sum_Detected_Conc_ng_L, na.rm = TRUE),
    
    Max_Total_PFAS = max(Sum_Detected_Conc_ng_L, na.rm = TRUE),
    
    Avg_Compounds_Detected = mean(N_Detected, na.rm = TRUE),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(desc(Avg_Total_PFAS))

cat("\n\nPFAS Summary by Location:\n")

print(location_summary)

# --- Summary by project ---

project_summary <- PFAS_by_sample %>%
  
  group_by(Project) %>%
  
  summarise(
    
    N_Samples = n(),
    
    Avg_Total_PFAS = mean(Sum_Detected_Conc_ng_L, na.rm = TRUE),
    
    Median_Total_PFAS = median(Sum_Detected_Conc_ng_L, na.rm = TRUE),
    
    Max_Total_PFAS = max(Sum_Detected_Conc_ng_L, na.rm = TRUE),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(desc(Avg_Total_PFAS))

cat("\nPFAS Summary by Project:\n")

print(project_summary)

# --- Compound detection frequency ---

compound_freq <- PFAS_long %>%
  
  group_by(Compound) %>%
  
  summarise(
    
    N_Measured = n(),
    
    N_Detected = sum(Detected),
    
    Detection_Rate = round(sum(Detected) / n() * 100, 1),
    
    Mean_Conc_Detected = mean(Concentration_ng_L[Detected], na.rm = TRUE),
    
    Max_Conc = max(Concentration_ng_L, na.rm = TRUE),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(desc(Detection_Rate))

cat("\n\nPFAS Compound Detection Frequency (Top 15):\n")

print(compound_freq %>% filter(Detection_Rate > 0) %>% head(15))


cat("\n\nPFAS Analysis Summary:\n")

cat("Total samples analyzed:", n_distinct(PFAS_long$Sample_ID), "\n")

cat("Total compound measurements:", nrow(PFAS_long), "\n")

cat("Total detections (above MDL):", sum(PFAS_long$Detected), "\n\n")

# ============================================================

# SECTION 8: INTEGRATED FLOW + PFAS VISUALIZATIONS

# ============================================================

cat("\n========== SECTION 8: INTEGRATED FLOW + PFAS VISUALIZATIONS ==========\n")

# --- CRITICAL FIX: Ensure consistent datetime formats ---

WetCenter_Flow_Combined$Time <- as.POSIXct(WetCenter_Flow_Combined$Time, tz = CONFIG$TIMEZONE)

PFAS_by_sample$DateTime <- as.POSIXct(PFAS_by_sample$DateTime, tz = CONFIG$TIMEZONE)

# --- Verify overlap ---

cat("\n\n=== DATE OVERLAP CHECK ===\n")

cat("Flow data covers:", format(min(WetCenter_Flow_Combined$Time, na.rm = TRUE), "%Y-%m-%d"),
    
    "to", format(max(WetCenter_Flow_Combined$Time, na.rm = TRUE), "%Y-%m-%d"), "\n")

cat("PFAS samples collected:", format(min(PFAS_by_sample$DateTime, na.rm = TRUE), "%Y-%m-%d"),
    
    "to", format(max(PFAS_by_sample$DateTime, na.rm = TRUE), "%Y-%m-%d"), "\n")

# --- Filter PFAS samples to flow period and interpolate flow ---

samples_in_range <- PFAS_by_sample %>%
  
  filter(DateTime >= min(WetCenter_Flow_Combined$Time, na.rm = TRUE),
         
         DateTime <= max(WetCenter_Flow_Combined$Time, na.rm = TRUE)) %>%
  
  mutate(
    
    Flow = approx(x = WetCenter_Flow_Combined$Time,
                  
                  y = WetCenter_Flow_Combined$Flow,
                  
                  xout = DateTime,
                  
                  rule = 2)$y
    
  ) %>%
  
  filter(!is.na(Flow))

cat("PFAS samples within flow period:", nrow(samples_in_range), "\n\n")


#  Clean and filter PFAS data for Wet Center 

# Filter to Wet Center only and remove zero/NA values
WetCenter_PFAS <- PFAS_by_sample %>%
  filter(Location == "Wet Center",
         Sum_Detected_Conc_ng_L > 0,
         !is.na(DateTime),
         !is.na(Sum_Detected_Conc_ng_L)) %>%
  arrange(DateTime)

cat("\n\nWet Center PFAS Samples:\n")
cat("Total samples:", nrow(WetCenter_PFAS), "\n")
cat("Date range:", format(min(WetCenter_PFAS$DateTime, na.rm = TRUE), "%Y-%m-%d"),
    "to", format(max(WetCenter_PFAS$DateTime, na.rm = TRUE), "%Y-%m-%d"), "\n\n")
print(WetCenter_PFAS)

# Interpolate flow for Wet Center samples
WetCenter_PFAS <- WetCenter_PFAS %>%
  mutate(
    Flow = approx(x = WetCenter_Flow_Combined$Time,
                  y = WetCenter_Flow_Combined$Flow,
                  xout = DateTime,
                  rule = 2)$y
  ) %>%
  filter(!is.na(Flow))

cat("\nWet Center PFAS samples with matched flow:", nrow(WetCenter_PFAS), "\n\n")


# --- Plot 1: Flow with PFAS samples overlaid ---

print(ggplot() +
        
        geom_line(data = WetCenter_Flow_Combined,
                  
                  aes(x = Time, y = Flow, color = Source),
                  
                  linewidth = 0.6) +
        
        geom_point(data = samples_in_range,
                   
                   aes(x = DateTime, y = Flow, size = N_Detected, fill = Sum_Detected_Conc_ng_L),
                   
                   shape = 21, alpha = 0.7, color = "black") +
        
        scale_color_manual(values = c("TimeView" = "darkorange",
                                      
                                      "Pressure Transducer" = "steelblue"),
                           
                           name = "Flow Source") +
        
        scale_size_continuous(name = "# Compounds\nDetected", range = c(1, 5)) +
        
        scale_fill_viridis_c(name = "Total PFAS\n(ng/L)") +
        
        labs(title = "Mill River — Combined Flow Record with PFAS Samples",
             
             x = "Date/Time", y = "Discharge (m³/s)") +
        
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month") +
        
        theme_mill())

# --- Plot 2 : Flow with events and Wet Center PFAS samples ---
print(ggplot() +
        geom_rect(data = events,
                  aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = type),
                  alpha = 0.25) +
        geom_line(data = WetCenter_Flow_Combined,
                  aes(x = Time, y = Flow, color = Source),
                  linewidth = 0.6) +
        geom_point(data = WetCenter_PFAS,
                   aes(x = DateTime, y = Flow, size = Sum_Detected_Conc_ng_L),
                   color = "red", alpha = 0.7, shape = 21, fill = "darkred") +
        geom_text(data = events,
                  aes(x = start + (end - start) / 2, y = Inf, label = label),
                  vjust = 1.5, size = 3.5, fontface = "bold") +
        scale_fill_manual(values = event_colors, name = "Event Type") +
        scale_color_manual(values = c("TimeView" = "darkorange",
                                      "Pressure Transducer" = "steelblue"),
                           name = "Flow Source") +
        scale_size_continuous(name = "Total PFAS\n(ng/L)", range = c(3, 10)) +
        labs(title = "Wet Center — Flow Record with Sampling Events and PFAS Samples",
             x = "Date/Time", y = "Discharge (m³/s)") +
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month", expand = c(0, 0)) +
        theme_mill())

# --- Plot 3: PFAS by location boxplot ---

print(ggplot(PFAS_by_sample, aes(x = Location, y = Sum_Detected_Conc_ng_L, fill = Location)) +
        
        geom_boxplot(alpha = 0.7) +
        
        geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
        
        scale_y_log10() +
        
        labs(title = "Total Detected PFAS Concentration by Location",
             
             y = "Sum of Detected PFAS (ng/L, log scale)", x = "Location") +
        
        theme_mill() +
        
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        
        theme(legend.position = "none"))

# --- Plot 4 : PFAS temporal trend - Wet Center only ---
print(ggplot(WetCenter_PFAS, aes(x = DateTime, y = Sum_Detected_Conc_ng_L, size = N_Detected)) +
        geom_point(color = "darkred", alpha = 0.7) +
        geom_line(color = "darkred", alpha = 0.4, linewidth = 0.6) +
        scale_y_log10() +
        labs(title = "Wet Center — Temporal Trend: Total PFAS Concentration",
             x = "Date/Time", y = "Sum of Detected PFAS (ng/L, log scale)",
             size = "# Compounds\nDetected") +
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month") +
        theme_mill())

# --- Plot 5: PFAS by project ---

print(ggplot(PFAS_by_sample, aes(x = Project, y = Sum_Detected_Conc_ng_L, fill = Project)) +
        
        geom_boxplot(alpha = 0.7) +
        
        geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
        
        scale_y_log10() +
        
        labs(title = "Total Detected PFAS Concentration by Project",
             
             y = "Sum of Detected PFAS (ng/L, log scale)", x = "Project") +
        
        theme_mill() +
        
        theme(legend.position = "none"))

# --- Plot 6: Key compounds across all samples ---

print(ggplot(PFAS_key %>% filter(Detected),
             
             aes(x = DateTime, y = Concentration_ng_L, color = Compound, shape = Location)) +
        
        geom_point(size = 3, alpha = 0.7) +
        
        scale_y_log10() +
        
        scale_color_viridis_d() +
        
        labs(title = "Temporal Trend: Key PFAS Compounds",
             
             x = "Date/Time", y = "Concentration (ng/L, log scale)",
             
             color = "Compound", shape = "Location") +
        
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month") +
        
        theme_mill())


# ============================================================



# --- Identify top 5 compounds by maximum concentration at Wet Center ---

top_5_compounds <- PFAS_long %>%
  
  filter(str_detect(Location, regex("Wet Center$", ignore_case = TRUE)),
         
         Detected == TRUE) %>%
  
  group_by(Compound) %>%
  
  summarise(
    
    Max_Conc = max(Concentration_ng_L, na.rm = TRUE),
    
    Mean_Conc = mean(Concentration_ng_L, na.rm = TRUE),
    
    N_Detections = n(),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(desc(Max_Conc)) %>%
  
  head(5) %>%
  
  pull(Compound)

cat("\nTop 5 Compounds by Maximum Concentration:\n")

print(
  
  PFAS_long %>%
    
    filter(str_detect(Location, regex("Wet Center$", ignore_case = TRUE)),
           
           Detected == TRUE) %>%
    
    group_by(Compound) %>%
    
    summarise(
      
      Max_Conc = max(Concentration_ng_L, na.rm = TRUE),
      
      Mean_Conc = mean(Concentration_ng_L, na.rm = TRUE),
      
      N_Detections = n(),
      
      .groups = "drop"
      
    ) %>%
    
    arrange(desc(Max_Conc)) %>%
    
    head(5)
  
)

# --- Get ALL Wet Center measurements (detected and non-detected) for top 5 compounds ---

WetCenter_All_Top5 <- PFAS_long %>%
  
  filter(str_detect(Location, regex("Wet Center$", ignore_case = TRUE)),
         
         Compound %in% top_5_compounds,
         
         !is.na(DateTime)) %>%
  
  mutate(
    
    Flow = approx(x = WetCenter_Flow_Combined$Time,
                  
                  y = WetCenter_Flow_Combined$Flow,
                  
                  xout = DateTime,
                  
                  rule = 2)$y
    
  ) %>%
  
  filter(!is.na(Flow))

# Separate into detected and non-detected

WetCenter_Top5_Detected <- WetCenter_All_Top5 %>% filter(Detected == TRUE)

WetCenter_Top5_ND <- WetCenter_All_Top5 %>%
  
  filter(Detected == FALSE) %>%
  
  mutate(Concentration_ng_L = NA_real_) # Will plot at MDL level

cat("\nWet Center Top 5 Compounds Summary:\n")

cat(" • Total measurements:", nrow(WetCenter_All_Top5), "\n")

cat(" • Detected:", nrow(WetCenter_Top5_Detected), "\n")

cat(" • Non-detects (ND/BDL):", nrow(WetCenter_Top5_ND), "\n\n")

# --- Plot 1: Dual Y-Axis with Non-Detects as hollow circles at MDL level ---

print(ggplot() +
        
        # Flow line (primary y-axis)
        
        geom_line(data = WetCenter_Flow_Combined,
                  
                  aes(x = Time, y = Flow, color = "Flow"),
                  
                  linewidth = 0.8, alpha = 0.6) +
        
        # DETECTED PFAS compounds (filled points)
        
        geom_point(data = WetCenter_Top5_Detected,
                   
                   aes(x = DateTime,
                       
                       y = 10^(log10(Concentration_ng_L)) * scale_factor,
                       
                       color = Compound,
                       
                       shape = Compound,
                       
                       fill = Compound),
                   
                   size = 4, alpha = 0.85, stroke = 0.8) +
        
        # NON-DETECTED (hollow points at MDL level)
        
        geom_point(data = WetCenter_Top5_ND,
                   
                   aes(x = DateTime,
                       
                       y = 10^(log10(MDL_ng_L)) * scale_factor,
                       
                       color = Compound),
                   
                   shape = 1, size = 2.5, alpha = 0.4, stroke = 0.8) +
        
        # Primary y-axis labels
        
        scale_y_continuous(
          
          name = "Discharge (m³/s)",
          
          breaks = scales::pretty_breaks(n = 6),
          
          expand = c(0, 0),
          
          limits = c(0, max_flow),
          
          # Secondary y-axis (right side) - PFAS log scale
          
          sec.axis = sec_axis(
            
            trans = ~. / scale_factor,
            
            name = "PFAS Concentration (ng/L, log scale)",
            
            breaks = c(0.1, 1, 10, 100),
            
            labels = c("0.1", "1", "10", "100"),
            
            transform = function(x) log10(x)
            
          )
          
        ) +
        
        scale_color_manual(
          
          values = c("Flow" = "steelblue",
                     
                     setNames(scales::hue_pal()(5), top_5_compounds)),
          
          name = "Compound"
          
        ) +
        
        scale_shape_manual(
          
          values = setNames(c(19, 17, 15, 18, 8), top_5_compounds),
          
          name = "Compound"
          
        ) +
        
        scale_fill_manual(
          
          values = c(setNames(scales::hue_pal()(5), top_5_compounds)),
          
          guide = "none"
          
        ) +
        
        labs(
          
          title = "Wet Center — Flow vs Top 5 PFAS Compounds with Non-Detects",
          
          subtitle = "Filled points = detected; Hollow circles = non-detect (at MDL level)",
          
          x = "Date/Time",
          
          y = "Discharge (m³/s)"
          
        ) +
        
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 month", expand = c(0, 0)) +
        
        theme_mill() +
        
        theme(
          
          axis.title.y = element_text(color = "steelblue", size = 11),
          
          axis.title.y.right = element_text(color = "darkred", size = 11),
          
          axis.text.y = element_text(color = "steelblue"),
          
          axis.text.y.right = element_text(color = "darkred"),
          
          legend.position = "right",
          
          legend.ncol = 1
          
        ))

# --- Plot 2 (UPDATED): Faceted view with lines connecting observations ---

print(ggplot() +
        
        facet_wrap(~Compound, ncol = 1) +
        
        # Flow line (normalized)
        
        geom_line(data = WetCenter_Flow_Combined,
                  
                  aes(x = Time, y = Flow / max(WetCenter_Flow_Combined$Flow, na.rm = TRUE) * 100),
                  
                  color = "steelblue", linewidth = 0.8, alpha = 0.4) +
        
        # DETECTED points - connected with line
        
        geom_line(data = WetCenter_Top5_Detected %>% arrange(DateTime),
                  
                  aes(x = DateTime, y = Concentration_ng_L, color = Compound),
                  
                  linewidth = 0.6, alpha = 0.5) +
        
        geom_point(data = WetCenter_Top5_Detected,
                   
                   aes(x = DateTime, y = Concentration_ng_L, color = Compound, fill = Compound),
                   
                   size = 4, alpha = 0.85, shape = 19, stroke = 0.5) +
        
        # NON-DETECTED points - connected with dashed line at MDL
        
        geom_line(data = WetCenter_Top5_ND %>% arrange(DateTime),
                  
                  aes(x = DateTime, y = MDL_ng_L, color = Compound),
                  
                  linewidth = 0.4, alpha = 0.3, linetype = "dotted") +
        
        geom_point(data = WetCenter_Top5_ND,
                   
                   aes(x = DateTime, y = MDL_ng_L, color = Compound),
                   
                   size = 2, alpha = 0.4, shape = 1, stroke = 0.8) +
        
        # Horizontal MDL reference line per compound
        
        geom_hline(data = WetCenter_Top5_ND %>% select(Compound, MDL_ng_L) %>% distinct(),
                   
                   aes(yintercept = MDL_ng_L),
                   
                   linetype = "dashed", color = "gray50", linewidth = 0.4, alpha = 0.5) +
        
        scale_y_log10(
          
          breaks = c(0.1, 1, 10, 100),
          
          labels = c("0.1", "1", "10", "100"),
          
          name = "Concentration (ng/L, log scale)"
          
        ) +
        
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 month") +
        
        scale_color_manual(values = scales::hue_pal()(5), guide = "none") +
        
        scale_fill_manual(values = scales::hue_pal()(5), guide = "none") +
        
        labs(
          
          title = "Wet Center — Top 5 PFAS Compounds with Temporal Patterns",
          
          subtitle = "Solid line = detected detections; Dotted line = non-detects at MDL; Filled = detected, Hollow = ND",
          
          x = "Date/Time",
          
          y = "Concentration (ng/L, log scale)"
          
        ) +
        
        theme_mill() +
        
        theme(
          
          strip.text = element_text(size = 12, face = "bold"),
          
          panel.spacing.y = unit(1.2, "cm")
          
        ))

# --- Alternative: Single plot with lines and all compounds overlaid ---

print(ggplot() +
        
        # Flow line (secondary reference)
        
        geom_line(data = WetCenter_Flow_Combined,
                  
                  aes(x = Time, y = Flow / max(WetCenter_Flow_Combined$Flow, na.rm = TRUE) * 100),
                  
                  color = "steelblue", linewidth = 0.8, alpha = 0.3) +
        
        # DETECTED observations with connecting lines
        
        geom_line(data = WetCenter_Top5_Detected %>% arrange(DateTime),
                  
                  aes(x = DateTime, y = Concentration_ng_L, color = Compound),
                  
                  linewidth = 0.7, alpha = 0.6) +
        
        geom_point(data = WetCenter_Top5_Detected,
                   
                   aes(x = DateTime, y = Concentration_ng_L, color = Compound),
                   
                   size = 3.5, alpha = 0.85, shape = 19, stroke = 0.3) +
        
        # NON-DETECTED observations with dotted lines at MDL
        
        geom_line(data = WetCenter_Top5_ND %>% arrange(DateTime),
                  
                  aes(x = DateTime, y = MDL_ng_L, color = Compound),
                  
                  linewidth = 0.4, alpha = 0.25, linetype = "dotted") +
        
        geom_point(data = WetCenter_Top5_ND,
                   
                   aes(x = DateTime, y = MDL_ng_L, color = Compound),
                   
                   size = 2, alpha = 0.35, shape = 1, stroke = 0.7) +
        
        scale_y_log10(
          
          breaks = c(0.1, 1, 10, 100),
          
          labels = c("0.1", "1", "10", "100"),
          
          name = "Concentration (ng/L, log scale)",
          
          limits = c(0.1, 150)
          
        ) +
        
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month") +
        
        scale_color_manual(
          
          values = scales::hue_pal()(5),
          
          name = "Compound"
          
        ) +
        
        labs(
          
          title = "Wet Center — Top 5 PFAS Compounds Temporal Patterns",
          
          subtitle = "Solid lines = detected; Dotted lines = non-detects at MDL; Light blue = normalized flow",
          
          x = "Date/Time",
          
          y = "Concentration (ng/L, log scale)"
          
        ) +
        
        scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month", expand = c(0, 0)) +
        
        theme_mill() +
        
        theme(
          
          legend.position = "right",
          
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
          
          panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.1)
          
        ))

# --- Enhanced seasonal view: Connect by Month ---

cat("\n\n=== SEASONAL PATTERN VISUALIZATION ===\n\n")

WetCenter_Top5_Monthly <- WetCenter_All_Top5 %>%
  
  mutate(
    
    Year_Month = floor_date(DateTime, "month"),
    
    Month_Name = format(Year_Month, "%b %Y")
    
  ) %>%
  
  group_by(Compound, Year_Month) %>%
  
  summarise(
    
    N_Detected = sum(Detected),
    
    N_ND = sum(!Detected),
    
    Mean_When_Detected = mean(Concentration_ng_L[Detected == TRUE], na.rm = TRUE),
    
    Max_Conc = max(Concentration_ng_L, na.rm = TRUE),
    
    MDL = first(MDL_ng_L),
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    
    Plot_Value = if_else(N_Detected > 0, Mean_When_Detected, MDL),
    
    Point_Type = if_else(N_Detected > 0, "Detected", "Non-Detect")
    
  )

print(ggplot(WetCenter_Top5_Monthly, aes(x = Year_Month, y = Plot_Value, color = Compound)) +
        
        geom_line(linewidth = 0.8, alpha = 0.7) +
        
        geom_point(aes(shape = Point_Type, fill = Point_Type),
                   
                   size = 4.5, alpha = 0.8, stroke = 1.2) +
        
        scale_y_log10(
          
          breaks = c(0.1, 1, 10, 100),
          
          labels = c("0.1", "1", "10", "100"),
          
          name = "Mean Concentration (ng/L, log scale)"
          
        ) +
        
        scale_shape_manual(
          
          values = c("Detected" = 19, "Non-Detect" = 1),
          
          name = "Status"
          
        ) +
        
        scale_fill_manual(
          
          values = c("Detected" = "black", "Non-Detect" = NA),
          
          guide = "none"
          
        ) +
        
        scale_color_manual(
          
          values = scales::hue_pal()(5),
          
          name = "Compound"
          
        ) +
        
        labs(
          
          title = "Wet Center — Monthly Mean PFAS Concentrations",
          
          subtitle = "Filled circles = detected (mean); Hollow = non-detect (at MDL); Lines connect monthly values",
          
          x = "Month",
          
          y = "Mean Concentration (ng/L, log scale)"
          
        ) +
        
        scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month") +
        
        theme_mill() +
        
        theme(
          
          axis.text.x = element_text(angle = 45, hjust = 1),
          
          legend.position = "right",
          
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2)
          
        ))
# --- Summary Table: Detection Frequency ---

cat("\n\nDetection Frequency by Compound:\n")

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

detection_freq <- WetCenter_All_Top5 %>%
  
  group_by(Compound) %>%
  
  summarise(
    
    Total_Samples = n(),
    
    Detected = sum(Detected),
    
    Non_Detect = sum(!Detected),
    
    Detection_Rate = round(sum(Detected) / n() * 100, 1),
    
    MDL = first(MDL_ng_L),
    
    Mean_Conc_When_Detected = round(mean(Concentration_ng_L[Detected == TRUE], na.rm = TRUE), 2),
    
    Max_Conc = round(max(Concentration_ng_L, na.rm = TRUE), 2),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(desc(Detection_Rate))

print(detection_freq)

# --- Summary: Seasonality of Non-Detects ---
seasonal_nd <- WetCenter_All_Top5 %>%
  
  mutate(Month = month(DateTime),
         
         Season = case_when(
           
           Month %in% c(12, 1, 2) ~ "Winter",
           
           Month %in% c(3, 4, 5) ~ "Spring",
           
           Month %in% c(6, 7, 8) ~ "Summer",
           
           Month %in% c(9, 10, 11) ~ "Fall"
           
         )) %>%
  
  group_by(Compound, Season) %>%
  
  summarise(
    
    Total = n(),
    
    Detected = sum(Detected),
    
    ND_Rate = round((n() - sum(Detected)) / n() * 100, 1),
    
    .groups = "drop"
    
  ) %>%
  
  pivot_wider(names_from = Season, values_from = c(Total, Detected, ND_Rate))

cat("\nNon-Detect Rate (%) by Season:\n")
print(seasonal_nd)

# --- Key Observation: PFOSA seasonality ---

cat("KEY OBSERVATION: PFOSA Seasonality\n")

pfosa_seasonal <- WetCenter_All_Top5 %>%
  
  filter(Compound == "PFOSA") %>%
  
  mutate(Month = month(DateTime),
         
         Year_Month = floor_date(DateTime, "month")) %>%
  
  group_by(Year_Month) %>%
  
  summarise(
    
    Detected = sum(Detected),
    
    ND = sum(!Detected),
    
    Mean_Conc = mean(Concentration_ng_L[Detected == TRUE], na.rm = TRUE),
    
    Max_Conc = max(Concentration_ng_L, na.rm = TRUE),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(Year_Month)

print(pfosa_seasonal)

# ANALYSIS COMPLETE

cat("\n========== That's All Folks ==========\n")

cat("\nKey PFAS Detections:\n")
top_detected <- compound_freq %>% filter(Detection_Rate > 0) %>% head(5)

for (i in 1:nrow(top_detected)) {
  
  row <- top_detected[i, ]
  
  cat(" •", row$Compound, ":", row$N_Detected, "detections",
      
      "(", row$Detection_Rate, "% detection rate)\n")
  
}


# ============================================================
# PFAS LOADING CALCULATION
# ============================================================

cat("\n========== PFAS LOADING ANALYSIS ==========\n")

# --- Step 1: Create daily averages of flow and PFAS concentrations ---

# Daily flow average
Daily_Flow <- WetCenter_Flow_Combined %>%
  mutate(Date = floor_date(Time, "day")) %>%
  group_by(Date) %>%
  summarise(
    Flow_m3_s = mean(Flow, na.rm = TRUE),
    Flow_m3_day = Flow_m3_s * 86400,  # Convert to m³/day
    Flow_L_day = Flow_m3_day * 1000,  # Convert to L/day
    .groups = "drop"
  )

cat("Daily flow summary:\n")
print(summary(Daily_Flow$Flow_m3_day))

# --- Step 2: Create daily average concentration for each compound ---

Daily_PFAS_Conc <- WetCenter_All_Top5 %>%
  mutate(Date = floor_date(DateTime, "day")) %>%
  group_by(Compound, Date) %>%
  summarise(
    N_Obs = n(),
    Conc_ng_L = mean(Concentration_ng_L, na.rm = TRUE),
    Conc_Min = min(Concentration_ng_L, na.rm = TRUE),
    Conc_Max = max(Concentration_ng_L, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n\nDaily concentration samples:\n")
print(Daily_PFAS_Conc %>% head(20))

# --- Step 3: Join flow and concentration data ---

PFAS_Loading_Daily <- Daily_PFAS_Conc %>%
  left_join(Daily_Flow, by = "Date") %>%
  filter(!is.na(Flow_L_day), !is.na(Conc_ng_L)) %>%
  mutate(
    # Calculate daily mass loading
    Loading_ng_day = Conc_ng_L * Flow_L_day,
    Loading_ug_day = Loading_ng_day / 1000,
    Loading_g_day = Loading_ug_day / 1000000,
    Loading_kg_day = Loading_g_day / 1000
  ) %>%
  arrange(Date, Compound)

cat("\n\nDaily PFAS Loading:\n")
print(PFAS_Loading_Daily %>% 
        select(Date, Compound, Conc_ng_L, Flow_m3_day, Loading_g_day) %>%
        head(20))

# --- Step 4: Summary statistics ---

cat("\n\n=== LOADING SUMMARY STATISTICS ===\n")

Loading_Summary <- PFAS_Loading_Daily %>%
  group_by(Compound) %>%
  summarise(
    N_Days = n(),
    Avg_Daily_Loading_g = mean(Loading_g_day, na.rm = TRUE),
    Max_Daily_Loading_g = max(Loading_g_day, na.rm = TRUE),
    Min_Daily_Loading_g = min(Loading_g_day, na.rm = TRUE),
    Total_Loading_kg = sum(Loading_kg_day, na.rm = TRUE),
    Median_Flow_m3_s = median(Flow_m3_s, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Loading_kg))

print(Loading_Summary)

# --- Step 5: Cumulative loading over time ---

PFAS_Loading_Cumulative <- PFAS_Loading_Daily %>%
  group_by(Compound) %>%
  arrange(Date) %>%
  mutate(
    Cumulative_kg = cumsum(Loading_kg_day)
  ) %>%
  ungroup()

# --- Plot: Daily loading over time ---
print(ggplot(PFAS_Loading_Daily, aes(x = Date, y = Loading_g_day, color = Compound, fill = Compound)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.7) +
        scale_y_continuous(name = "Daily PFAS Loading (g/day)") +
        scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month") +
        scale_color_manual(values = scales::hue_pal()(5)) +
        scale_fill_manual(values = scales::hue_pal()(5)) +
        labs(
          title = "Daily PFAS Loading at Wet Center",
          subtitle = "Stacked by compound; Top 5 that contribute the most to mass flux",
          x = "Date"
        ) +
        theme_mill() +
        theme(
          legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1)
        ))

# --- Plot: Cumulative loading ---
print(ggplot(PFAS_Loading_Cumulative, aes(x = Date, y = Cumulative_kg, color = Compound)) +
        geom_line(linewidth = 1, alpha = 0.8) +
        scale_y_continuous(name = "Cumulative PFAS Loading (kg)") +
        scale_x_datetime(date_labels = "%b %y", date_breaks = "1 month") +
        scale_color_manual(values = scales::hue_pal()(5), name = "Compound") +
        labs(
          title = "Cumulative PFAS Loading at Wet Center",
          subtitle = "Shows total mass exported over monitoring period",
          x = "Date"
        ) +
        theme_mill() +
        theme(
          legend.position = "right",
          panel.grid.major = element_line(color = "gray90")
        ))

# --- Event-based loading: Calculate loading during high flow periods ---

# Define high flow as > 75th percentile
flow_75 = quantile(Daily_Flow$Flow_m3_s, 0.75, na.rm = TRUE)

cat("\n\n=== EVENT VS BASEFLOW LOADING ===\n")
cat("High flow threshold (75th percentile):", round(flow_75, 2), "m³/s\n\n")

Event_Loading <- PFAS_Loading_Daily %>%
  mutate(
    Flow_Type = if_else(Flow_m3_s > flow_75, "High Flow", "Baseflow")
  ) %>%
  group_by(Compound, Flow_Type) %>%
  summarise(
    N_Days = n(),
    Avg_Daily_Loading_g = mean(Loading_g_day, na.rm = TRUE),
    Total_Loading_kg = sum(Loading_kg_day, na.rm = TRUE),
    Pct_Total = round(Total_Loading_kg / sum(PFAS_Loading_Daily$Loading_kg_day[PFAS_Loading_Daily$Compound == Compound]) * 100, 1),
    .groups = "drop"
  )

print(Event_Loading)

# --- Plot: Loading by flow regime ---
print(ggplot(Event_Loading, aes(x = Compound, y = Total_Loading_kg, fill = Flow_Type)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        scale_y_continuous(name = "Total PFAS Loading (kg)") +
        scale_fill_manual(values = c("High Flow" = "steelblue", "Baseflow" = "green4"), name = "Flow Regime") +
        labs(
          title = "PFAS Loading: High Flow vs Baseflow",
          subtitle = "Shows which flow conditions export most mass",
          x = "Compound"
        ) +
        theme_mill() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1)
        ))

# --- Concentration vs Loading correlation ---

cat("\n\n=== CONCENTRATION VS LOADING CORRELATION ===\n")

Loading_Flow_Analysis <- PFAS_Loading_Daily %>%
  group_by(Compound) %>%
  summarise(
    Corr_Conc_Loading = cor(Conc_ng_L, Loading_g_day, use = "complete.obs"),
    Corr_Flow_Loading = cor(Flow_m3_s, Loading_g_day, use = "complete.obs"),
    .groups = "drop"
  )

print(Loading_Flow_Analysis)

cat("\nInterpretation:\n")
cat("  • High Corr_Conc_Loading = Loading driven by concentration (source control)\n")
cat("  • High Corr_Flow_Loading = Loading driven by flow (mobilization from environment)\n\n")

