library(tidyverse)
library(lubridate)
library(readr)

# --- 1. Load Data ---
water_abs_df <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation4_WaterPressure.csv")
air_abs_df   <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation3_AirPressure.csv")
gauge_df     <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenter%20GaugeHeight%20-%20Sheet1-2.csv")

# --- 3. Parse Timestamps at Full 10-Minute Resolution ---
water_ts <- water_abs_df %>%
  mutate(DateTime = parse_date_time(`Date-Time (EDT)`, orders = c("mdy HM", "mdy HMS"))) %>%
  rename(Water_Abs_kPa = `Absolute Pressure kPa`) %>%
  select(DateTime, Water_Abs_kPa)

air_ts <- air_abs_df %>%
  mutate(DateTime = parse_date_time(`Date-Time (EDT)`, orders = c("mdy HM", "mdy HMS"))) %>%
  rename(Air_Abs_kPa = `Absolute Pressure kPa`) %>%
  select(DateTime, Air_Abs_kPa)

# --- 4. Join Air + Water on Exact Timestamp ---
# rolling nearest-match join 
library(data.table)
setDT(water_ts); setDT(air_ts)
sensor_ts <- air_ts[water_ts, on = "DateTime", roll = "nearest"] %>%
   mutate(Water_Pressure_Only = Water_Abs_kPa - Air_Abs_kPa)

# --- 5. Parse Gauge Readings with Date + Time ---
# *** ACTION REQUIRED: Make sure your gauge CSV has a time column! ***
# If it currently only has dates, add the time of measurement in the field notes.
# Expected format: "4/15/2025 14:30"
 # --- 5. Parse Gauge Readings — combine Date + Time columns ---
gauge_clean <- gauge_df %>%
  filter(!is.na(`GaugeHeight (ft)`)) %>%
  mutate(
    DateTime = as.POSIXct(paste(Date, format(Time, "%H:%M")),
                          format = "%m/%d/%y %H:%M",
                          tz = "America/New_York")
  ) %>%
  select(DateTime, `GaugeHeight (ft)`)
 

# --- 6. Match Each Gauge Reading to the Nearest Sensor Record 
library(data.table)

sensor_dt <- as.data.table(sensor_ts)
gauge_dt  <- as.data.table(gauge_clean)

# Rolling nearest-time join
calib_points <- sensor_dt[gauge_dt, on = "DateTime", roll = "nearest"] %>%
  as_tibble() %>%
  filter(!is.na(Water_Pressure_Only), !is.na(`GaugeHeight (ft)`))

print(paste("Calibration points matched:", nrow(calib_points)))

# --- 7. Build Calibration Model ---
calib_model <- lm(`GaugeHeight (ft)` ~ Water_Pressure_Only, data = calib_points)
summary(calib_model)
# Check R-squared — should be > 0.95 for a good pressure-height relationship

# --- 8. Apply Model to Full 10-Minute Time Series ---
sensor_calibrated <- sensor_ts %>%
  mutate(
    Calibrated_Height_ft = predict(calib_model,
                                   newdata = data.frame(Water_Pressure_Only = Water_Pressure_Only))
  )

# --- 9. Plot ---
ggplot() +
  geom_line(data = sensor_calibrated,
            aes(x = DateTime, y = Calibrated_Height_ft, color = "10-min Sensor"),
            linewidth = 0.5, alpha = 0.8) +
  geom_point(data = calib_points,
             aes(x = DateTime, y = `GaugeHeight (ft)`, color = "Manual Gauge"),
             size = 3, shape = 21, fill = "red") +
  scale_color_manual(values = c("10-min Sensor" = "steelblue", "Manual Gauge" = "red")) +
  labs(title = "Stream Height — 10-Minute Resolution",
       subtitle = "Barometric-compensated pressure, calibrated to manual gauge readings",
       x = "Date/Time (EDT)", y = "Water Level (ft)",
       color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")



# --- Rating Curve Data ---
rating_data <- tibble(
  Stage_m = c(0.23, 0.27, 0.58, 0.58, 0.26, 0.26, 0.37, 0.36,
              0.29, 0.29, 0.49, 0.48),
  Flow_m3s = c(0.34, 0.42, 1.125, 1.146, 0.179, 0.176, 0.522,
               0.505, 0.302, 0.299, 0.755, 0.724)
) %>%
  # Average duplicate measurements taken at the same visit
  group_by(Stage_m) %>%
  summarise(Flow_m3s = mean(Flow_m3s), .groups = "drop")

# --- Fit Power Law via NLS ---
rating_model <- nls(Flow_m3s ~ a * Stage_m^b,
                    data = rating_data,
                    start = list(a = 2.74, b = 1.73))

summary(rating_model)

a <- coef(rating_model)["a"]
b <- coef(rating_model)["b"]
print(paste("Q =", round(a, 3), "* H^", round(b, 3)))

# --- Apply to WetStation4Stream ---
WetStation4Stream <- sensor_calibrated %>%
  rename(Time = DateTime,
         GaugeHeight = Calibrated_Height_ft) %>%
  mutate(
    GaugeHeight_m = GaugeHeight * 0.3048,
    Flow = a * (GaugeHeight_m ^ b)
  ) %>%
  arrange(Time)

# --- Plot Rating Curve for Inspection ---
stage_seq <- seq(0.20, 0.65, by = 0.01)
predicted <- tibble(
  Stage_m = stage_seq,
  Flow_m3s = a * stage_seq^b
)



ggplot() +
  geom_line(data = WetStation4Stream,
            aes(x = Time, y = Flow, color = "10-min Sensor"),
            linewidth = 0.5, alpha = 0.8) +
  geom_point(data = calib_points,
             aes(x = DateTime, y = a * (`GaugeHeight (ft)` * 0.3048)^b, color = "Manual Gauge"),
             size = 3, shape = 21, fill = "red") +
  scale_color_manual(values = c("10-min Sensor" = "steelblue", "Manual Gauge" = "red")) +
  labs(title = "Stream Flow — 10-Minute Resolution",
       subtitle = "Barometric-compensated pressure, calibrated to manual gauge readings",
       x = "Date/Time (EDT)", y = "Discharge (m³/s)",
       color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")



# ============================================================
# PFAS + HYDROGRAPH PIPELINE
# ============================================================

# --- Load Data ---
PFAS_Inventory_temp        <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PFAS%20Inventory%20-%20temp.csv")
Storm16_PFAS               <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/RV_PFAS_Storm16.csv")
PrecipitationData_Storms16 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PrecipitationData_Storms16.csv")

# --- Parse DateTimes ---
PrecipitationData_Storms16 <- PrecipitationData_Storms16 %>%
  mutate(DateTime = mdy_hm(DateTime, tz = "America/New_York"))

Storm16_PFAS <- Storm16_PFAS %>%
  mutate(DateTime = as.POSIXct(paste(Date, format(Time, "%H:%M")),
                               format = "%m/%d/%Y %H:%M",
                               tz = "America/New_York"))

PFAS_Inventory_temp <- PFAS_Inventory_temp %>%
  mutate(
    DateTime = as.POSIXct(paste(Date, format(Time, "%H:%M")),
                          format = "%m/%d/%Y %H:%M",
                          tz = "America/New_York"),
    StormID = str_extract(Comment, regex("Storm \\d+", ignore_case = TRUE)) %>%
      str_to_title()
  )

# --- Join PFOA concentrations by DateTime ---
PFAS_joined <- PFAS_Inventory_temp %>%
  filter(!is.na(StormID)) %>%
  left_join(Storm16_PFAS %>% select(DateTime, PFOA), by = "DateTime")

# Sanity check
print(paste("PFAS joined rows:", nrow(PFAS_joined)))
print(paste("Rows with PFOA:", sum(!is.na(PFAS_joined$PFOA))))

# --- Deduplicate WetStation4Stream ---
WetStation4Stream <- WetStation4Stream %>%
  group_by(Time) %>%
  summarise(Flow         = mean(Flow, na.rm = TRUE),
            GaugeHeight_m = mean(GaugeHeight_m, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(Time)

# --- Plotting Function ---
plot_storm_hydrograph <- function(storm_id, pfas_data, flow_data,
                                  precip_data = NULL, pad_hours = 24) {
  
  storm_samples <- pfas_data %>%
    filter(StormID == storm_id, !is.na(PFOA)) %>%
    mutate(DateTime = as.POSIXct(DateTime, tz = "America/New_York"))
  
  if (nrow(storm_samples) == 0) {
    message("No PFOA data found for ", storm_id)
    return(NULL)
  }
  
  t_min <- min(storm_samples$DateTime) - hours(pad_hours)
  t_max <- max(storm_samples$DateTime) + hours(pad_hours)
  
  flow_window <- flow_data %>%
    filter(Time >= t_min, Time <= t_max)
  
  storm_samples <- storm_samples %>%
    mutate(Flow = approx(
      x    = flow_window$Time,
      y    = flow_window$Flow,
      xout = DateTime)$y)
  
  max_flow   <- max(flow_window$Flow, na.rm = TRUE) * 1.1
  max_pfoa   <- max(storm_samples$PFOA, na.rm = TRUE) * 1.1
  scale_pfoa <- max_flow / max_pfoa
  
  p <- ggplot() +
    geom_line(data = flow_window,
              aes(x = Time, y = Flow),
              color = "blue3", linewidth = 1.5) +
    geom_point(data = storm_samples,
               aes(x = DateTime, y = PFOA * scale_pfoa, color = PFOA),
               size = 5) +
    geom_label(data = storm_samples,
               aes(x = DateTime, y = PFOA * scale_pfoa,
                   label = round(PFOA, 1)),
               vjust = -0.7, size = 3.5, fontface = "bold",
               fill = "white", label.size = 0.2) +
    scale_color_viridis_c(name = "PFOA (ng/L)", option = "plasma") +
    scale_y_continuous(
      name     = "Discharge (m³/s)",
      limits   = c(0, max_flow),
      sec.axis = sec_axis(~ . / scale_pfoa, name = "PFOA (ng/L)"),
      expand   = c(0, 0)) +
    scale_x_datetime(date_labels = "%b %d", date_breaks = "1 day",
                     expand = c(0, 0)) +
    labs(title = paste(storm_id, "Hydrograph + PFOA"),
         x = "Date/Time (EDT)") +
    theme_bw() +
    theme(
      plot.title      = element_text(size = 20, face = "bold"),
      axis.title      = element_text(size = 16),
      axis.text       = element_text(size = 13),
      legend.position = "bottom"
    )
  
  if (!is.null(precip_data)) {
    max_precip   <- max(precip_data$`Precip (in)`, na.rm = TRUE)
    scale_precip <- (max_flow * 0.3) / max_precip
    
    p <- p +
      geom_rect(
        data = precip_data %>% filter(DateTime >= t_min, DateTime <= t_max),
        aes(xmin = DateTime - seconds(7.5 * 60),
            xmax = DateTime + seconds(7.5 * 60),
            ymin = max_flow - `Precip (in)` * scale_precip,
            ymax = max_flow),
        fill = "steelblue", alpha = 0.7)
  }
  
  return(p)
}

# --- Auto-generate plots for all storms with PFOA data ---
storms_with_data <- PFAS_joined %>%
  filter(!is.na(PFOA)) %>%
  pull(StormID) %>%
  unique()

print(paste("Storms to plot:", paste(storms_with_data, collapse = ", ")))

storm_plots <- map(storms_with_data, plot_storm_hydrograph,
                   pfas_data   = PFAS_joined,
                   flow_data   = WetStation4Stream,
                   precip_data = PrecipitationData_Storms16)

# --- View ---
walk(storm_plots, print)


