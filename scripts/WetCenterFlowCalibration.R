library(tidyverse)
library(lubridate)
library(readr)
library(ggplot2)
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
  labs(title = "Stream Height",
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
  labs(title = "Stream Flow",
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

library(highcharter)

# --- First, join Flow to Storm16_PFAS by nearest timestamp ---
Storm16_PFAS <- Storm16_PFAS %>%
  mutate(Flow = approx(
    x    = WetStation4Stream$Time,
    y    = WetStation4Stream$Flow,
    xout = DateTime)$y)

# --- Highcharter Plot ---
highchart() %>%
  hc_title(text = "Storm 16 PFAS + Hydrograph") %>%
  hc_yAxis_multiples(
    list(title = list(text = "Flow (m³/s)"),    lineWidth = 3, lineColor = "blue",
         min = min(Storm16_PFAS$Flow,  na.rm = TRUE), max = max(Storm16_PFAS$Flow,  na.rm = TRUE)),
    list(title = list(text = "PFOA (ng/L)"),    lineWidth = 3, lineColor = "red",
         min = min(Storm16_PFAS$PFOA,  na.rm = TRUE), max = max(Storm16_PFAS$PFOA,  na.rm = TRUE)),
    list(title = list(text = "PFHxA (ng/L)"),   lineWidth = 3, lineColor = "orange",
         min = min(Storm16_PFAS$PFHxA, na.rm = TRUE), max = max(Storm16_PFAS$PFHxA, na.rm = TRUE))
  ) %>%
  hc_plotOptions(series = list(connectNulls = TRUE)) %>%
  hc_add_series(
    name = "Flow",
    data = list_parse2(Storm16_PFAS %>% transmute(x = datetime_to_timestamp(DateTime), y = Flow)),
    yAxis = 0, type = "line", lineWidth = 4, color = "blue") %>%
  hc_add_series(
    name = "PFOA",
    data = list_parse2(Storm16_PFAS %>% transmute(x = datetime_to_timestamp(DateTime), y = PFOA)),
    yAxis = 1, type = "line", lineWidth = 3, color = "red") %>%
#  hc_add_series(
#    name = "PFHxA",
#    data = list_parse2(Storm16_PFAS %>% transmute(x = datetime_to_timestamp(DateTime), y = PFHxA)),
#    yAxis = 2, type = "line", lineWidth = 3, color = "orange") %>%
  hc_xAxis(type = "datetime", title = list(text = "Date/Time (EDT)"),
           dateTimeLabelFormats = list(day = "%b %d")) %>%
  hc_tooltip(shared = TRUE, crosshairs = TRUE,
             valueDecimals = 2) %>%
  hc_legend(enabled = TRUE)


#adding flow together. 

PFAS_WetCenter_2025 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenterPFASResults.csv")
TimeViewDischarge <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/2025_1YEAR_WET%20Center%20Umass%20AmherstWETCENTERUMASSMILL%20RIVER%20LEVEL%20SENSOR1_yearchannel1.csv")
AllChem_WetCenter_SE15 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Storm_Event_15_all_results.csv")

# --- Reparse WetCenterDischarge$Time ---
TimeViewDischarge <- TimeViewDischarge %>%
  mutate(Time = mdy_hm(DateTime, tz = "America/New_York"))

# --- Check what failed to parse ---
TimeViewDischarge %>%
  filter(is.na(Time)) %>%
  select(DateTime)

# --- Fix column reference and compute flow ---
WetCenterDischarge <- TimeViewDischarge %>%
  mutate(
    Time          = mdy_hm(DateTime, tz = "America/New_York"),
    GaugeHeight_m = `Stage (m)`,          # backticks for column reference, not quotes
    Flow          = a * (GaugeHeight_m ^ b)
  ) %>%
  filter(!is.na(Time), !is.na(GaugeHeight_m))  # drop the 4 failed rows

# --- Verify ---
summary(WetCenterDischarge$GaugeHeight_m)
summary(WetCenterDischarge$Flow)
range(WetCenterDischarge$Time, na.rm = TRUE)

# --- Define overlap period ---
overlap_start <- max(min(WetCenterDischarge$Time, na.rm = TRUE),
                     min(WetStation4Stream$Time,  na.rm = TRUE))
overlap_end   <- min(max(WetCenterDischarge$Time, na.rm = TRUE),
                     max(WetStation4Stream$Time,  na.rm = TRUE))

print(paste("Overlap:", overlap_start, "to", overlap_end))

# --- Match by nearest timestamp during overlap ---
library(data.table)

# --- Remove outlier readings at 15m ---
WetCenterDischarge <- WetCenterDischarge %>%
  filter(GaugeHeight_m < 14.9)  # slightly below 15 to catch any floating point near-15 values

pt_dt <- WetStation4Stream %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  select(Time, GaugeHeight_m) %>%
  rename(H_pt = GaugeHeight_m) %>%
  mutate(Time = as.POSIXct(Time, tz = "America/New_York")) %>%  # fix timezone warning
  as.data.table()

# --- Rerun overlap matching with clean data ---
tv_dt <- WetCenterDischarge %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  select(Time, GaugeHeight_m) %>%
  as.data.table()

overlap_matched <- tv_dt[pt_dt, on = "Time", roll = "nearest"] %>%
  as_tibble() %>%
  rename(H_tv = GaugeHeight_m) %>%
  filter(!is.na(H_tv), !is.na(H_pt))

# --- Refit regression ---
overlap_model <- lm(H_pt ~ H_tv, data = overlap_matched)
summary(overlap_model)

intercept <- coef(overlap_model)["(Intercept)"]
slope     <- coef(overlap_model)["H_tv"]

print(paste("Correction: H_corrected =", round(slope, 6),
            "* H_tv +", round(intercept, 4)))

# --- Apply correction and recompute flow ---
WetCenterDischarge <- WetCenterDischarge %>%
  mutate(
    GaugeHeight_m = slope * GaugeHeight_m + intercept,
    Flow          = a * (GaugeHeight_m ^ b)
  )

# --- Verify corrected range looks sensible ---
summary(WetCenterDischarge$GaugeHeight_m)
summary(WetCenterDischarge$Flow)
# --- Visualize overlap before correction ---
ggplot(overlap_matched, aes(x = H_tv, y = H_pt)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "TimeView vs Pressure Transducer — Overlap Period",
       x = "TimeView Stage (raw)", y = "Pressure Transducer Stage (m)") +
  theme_bw()


# --- Check correlation direction ---
cor(overlap_matched$H_tv, overlap_matched$H_pt)


# --- Visual check of corrected overlap ---
ggplot() +
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
  theme_bw() +
  theme(legend.position = "bottom")

# --- Trim WetCenterDischarge to end where pressure transducer begins ---
# --- Rebuild combined dataset ---
WetCenterDischarge_trimmed <- WetCenterDischarge %>%
  filter(Time < overlap_start) %>%
  transmute(
    Time,
    GaugeHeight_m,
    Flow,
    Source = "TimeView"
  )

WetStation4Stream_labeled <- WetStation4Stream %>%
  slice(-1) %>%
  transmute(
    Time,
    GaugeHeight_m,
    Flow,
    Source = "Pressure Transducer"
  )

WetCenter_Flow_Combined <- bind_rows(WetCenterDischarge_trimmed,
                                     WetStation4Stream_labeled) %>%
  arrange(Time)

# --- Sanity checks ---
range(WetCenter_Flow_Combined$Time, na.rm = TRUE)
table(WetCenter_Flow_Combined$Source)
summary(WetCenter_Flow_Combined$Flow)

# --- Full record plot ---
ggplot(WetCenter_Flow_Combined, aes(x = Time, y = Flow, color = Source)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue")) +
  labs(title = "Mill River — Combined Flow Record",
       x = "Date/Time", y = "Discharge (m³/s)") +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month",
                   expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))



# --- Define your event periods ---
events <- tribble(
  ~label,        ~start,                    ~end,                      ~type,
  "Storm 15",   "2025-10-10 00:00",        "2025-10-16 00:00",        "Storm",
  "Storm 16",   "2025-10-29 20:00",        "2025-11-07 18:00",        "Storm",
  # add other storms here as you go
  "12", "2025-07-09 00:00",        "2025-07-12 10:20",        "Storm",
  "13", "2025-07-14 05:00",        "2025-07-16 13:00",        "Storm"
) %>%
  mutate(
    start = as.POSIXct(start, tz = "America/New_York"),
    end   = as.POSIXct(end,   tz = "America/New_York")
  )

# --- Color palette per event type ---
event_colors <- c(
  "Storm"    = "steelblue",
  "Baseflow" = "forestgreen",
  "Snowmelt" = "mediumpurple"
)

# --- Plot ---
ggplot() +
  # Shaded event windows — drawn first so they sit behind the hydrograph
  geom_rect(data = events,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = type),
            alpha = 0.25) +
  # Hydrograph line
  geom_line(data = WetCenter_Flow_Combined,
            aes(x = Time, y = Flow, color = Source),
            linewidth = 0.6) +
  # Event labels at top of each shaded window
  geom_text(data = events,
            aes(x = start + (end - start) / 2,  # centered in window
                y = Inf, label = label),
            vjust = 1.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = event_colors, name = "Event Type") +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue"),
                     name = "Source") +
  labs(title = "Mill River — Combined Flow Record with Sampling Events",
       x = "Date/Time", y = "Discharge (m³/s)") +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month",
                   expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position  = "bottom",
        axis.text.x      = element_text(angle = 45, hjust = 1))



# --- Load inventory ---
PFAS_inventory_full <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PFAS%20Inventory%20-%20temp.csv")

# --- Parse DateTime ---
PFAS_inventory_full <- PFAS_inventory_full %>%
  mutate(
    DateTime = as.POSIXct(paste(Date, format(Time, "%H:%M")),
                          format = "%m/%d/%y %H:%M",
                          tz = "America/New_York")
  )

# --- Filter to Wet Center stream samples only ---
WetCenter_samples <- PFAS_inventory_full %>%
  filter(str_detect(Location, regex("wet center$", ignore_case = TRUE)),
         Type == "Stream") %>%
  filter(!is.na(DateTime)) %>%
  # Snap to nearest flow value
  mutate(Flow = approx(
    x    = WetCenter_Flow_Combined$Time,
    y    = WetCenter_Flow_Combined$Flow,
    xout = DateTime,
    rule = 2)$y)

# --- Plot ---
ggplot() +
  geom_rect(data = events,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf,
                fill = type),
            alpha = 0.25) +
  geom_line(data = WetCenter_Flow_Combined,
            aes(x = Time, y = Flow, color = Source),
            linewidth = 0.6) +
  geom_point(data = WetCenter_samples,
             aes(x = DateTime, y = Flow),
             color = "red", size = 2.5, alpha = 0.8) +
  geom_text(data = events,
            aes(x = start + (end - start) / 2,
                y = Inf, label = label),
            vjust = 1.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = event_colors, name = "Event Type") +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue"),
                     name = "Source") +
  labs(title = "Mill River — Combined Flow Record with PFAS Samples",
       x = "Date/Time", y = "Discharge (m³/s)") +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month",
                   expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position  = "bottom",
        axis.text.x      = element_text(angle = 45, hjust = 1))



# ============================================================
# STORM DETECTION + HYDROGRAPH SEPARATION
# ============================================================

library(zoo)  # for rolling mean
storm_find <- function(flow_data,
                       time_col          = "Time",
                       flow_col          = "Flow",
                       baseline_days     = 3,
                       rise_factor       = 1.5,
                       min_duration_hrs  = 6,
                       recession_factor  = 1.4,
                       pre_event_days    = 3,
                       min_gap_hrs       = 12) {
  
  df <- flow_data %>%
    rename(Time = all_of(time_col), Flow = all_of(flow_col)) %>%
    arrange(Time)
  
  baseline_n <- baseline_days * 24 * 6
  
  df <- df %>%
    mutate(
      Baseline        = rollapply(Flow, width = baseline_n,
                                  FUN = median, na.rm = TRUE,
                                  align = "right", fill = NA),
      Above_threshold = Flow > (Baseline * rise_factor)
    )
  
  message(paste("Baseline range:",
                round(min(df$Baseline, na.rm=TRUE), 3), "to",
                round(max(df$Baseline, na.rm=TRUE), 3)))
  message(paste("Points above threshold:", sum(df$Above_threshold, na.rm=TRUE)))
  
  df <- df %>%
    mutate(block_id = cumsum(!Above_threshold | is.na(Above_threshold)))
  
  storm_blocks <- df %>%
    filter(Above_threshold) %>%
    group_by(block_id) %>%
    summarise(
      peak_time    = Time[which.max(Flow)],
      peak         = max(Flow),
      thresh_start = min(Time),
      thresh_end   = max(Time),
      duration_hrs = as.numeric(difftime(max(Time), min(Time), units = "hours")),
      .groups      = "drop"
    ) %>%
    filter(duration_hrs >= min_duration_hrs)
  
  if (nrow(storm_blocks) == 0) {
    message("No storms detected — try lowering rise_factor")
    return(list(storms = storm_blocks, flow_flagged = df))
  }
  
  # Merge storms too close together
  if (nrow(storm_blocks) > 1) {
    storm_blocks <- storm_blocks %>%
      arrange(peak_time) %>%
      mutate(
        gap_to_next = as.numeric(difftime(lead(thresh_start), thresh_end,
                                          units = "hours")),
        new_storm   = is.na(lag(gap_to_next)) | lag(gap_to_next) > min_gap_hrs
      ) %>%
      mutate(storm_id = cumsum(new_storm)) %>%
      group_by(storm_id) %>%
      summarise(
        peak_time    = peak_time[which.max(peak)],
        peak         = max(peak),
        thresh_start = min(thresh_start),
        thresh_end   = max(thresh_end),
        duration_hrs = as.numeric(difftime(max(thresh_end), min(thresh_start),
                                           units = "hours")),
        .groups      = "drop"
      )
  }
  
  storm_blocks <- storm_blocks %>%
    mutate(storm_id = row_number())
  
  # --- Trough detection between consecutive storms ---
  # For each pair of adjacent storms, find the minimum flow between their
  # peaks and use that time as a hard boundary
  n <- nrow(storm_blocks)
  
  trough_times <- map_dbl(seq_len(n), ~ {
    if (.x == n) return(NA)  # no trough after last storm
    
    between <- df %>%
      filter(Time > storm_blocks$peak_time[.x],
             Time < storm_blocks$peak_time[.x + 1])
    
    if (nrow(between) == 0) return(NA)
    
    # Trough = minimum flow between the two peaks
    as.numeric(between$Time[which.min(between$Flow)])
  }) %>%
    as.POSIXct(origin = "1970-01-01", tz = "America/New_York")
  
  storm_blocks <- storm_blocks %>%
    mutate(
      trough_after  = trough_times,  # trough between this storm and the next
      trough_before = lag(trough_after),  # trough between previous storm and this one
      
      baseline_flow = map_dbl(peak_time, ~ {
        df %>%
          filter(Time >= .x - days(baseline_days),
                 Time <  .x - days(baseline_days - 1)) %>%
          summarise(m = median(Flow, na.rm = TRUE)) %>%
          pull(m)
      }),
      
      # Start: either pre_event_days before peak OR trough before this storm
      # whichever is LATER — prevents overlap with previous storm
      start = map2_dbl(peak_time, trough_before, ~ {
        candidate <- .x - days(pre_event_days)
        if (!is.na(.y)) {
          as.numeric(max(candidate, .y))
        } else {
          as.numeric(candidate)
        }
      }) %>%
        as.POSIXct(origin = "1970-01-01", tz = "America/New_York"),
      
      # End: recession to 1.4x baseline OR trough after this storm
      # whichever is EARLIER — prevents overlap with next storm
      end = pmap_dbl(list(peak_time, baseline_flow, trough_after), ~ {
        peak         <- ..1
        base         <- ..2
        trough_next  <- ..3
        recession_threshold <- base * recession_factor
        
        after_peak <- df %>% filter(Time > peak)
        recovered  <- after_peak %>% filter(Flow <= recession_threshold)
        
        recession_end <- if (nrow(recovered) == 0) {
          as.numeric(max(after_peak$Time))
        } else {
          as.numeric(min(recovered$Time))
        }
        
        # If there's a trough before the next storm, use whichever comes first
        if (!is.na(trough_next)) {
          min(recession_end, as.numeric(trough_next))
        } else {
          recession_end
        }
      }) %>%
        as.POSIXct(origin = "1970-01-01", tz = "America/New_York")
    )
  
  return(list(storms = storm_blocks, flow_flagged = df))
}

# --- Rerun ---
result <- storm_find(
  WetCenter_Flow_Combined,
  baseline_days    = 3,
  rise_factor      = 1.1,
  min_duration_hrs = 6,
  pre_event_days   = 3,
  recession_factor = 1.05,
  min_gap_hrs      = 12
)

storms_detected <- result$storms
flow_flagged    <- result$flow_flagged

# --- Check for any remaining overlaps ---
storms_detected %>%
  arrange(start) %>%
  mutate(overlap = start < lag(end)) %>%
  select(storm_id, start, peak_time, end, overlap)
# --- 2. Plot detected storms on hydrograph ---
ggplot() +
  geom_rect(data = storms_detected,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf),
            fill = "steelblue", alpha = 0.2) +
  geom_line(data = WetCenter_Flow_Combined,
            aes(x = Time, y = Flow, color = Source),
            linewidth = 0.6) +
  geom_point(data = storms_detected,
             aes(x = peak_time, y = peak),
             color = "blue", size = 3) +
  geom_point(data = WetCenter_samples,
             aes(x = DateTime, y = Flow),
             color = "red", size = 2, alpha = 0.7) +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue")) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "1 month",
                   expand = c(0, 0)) +
  labs(title = "Detected Storms — Mill River",
       x = "Date/Time", y = "Discharge (m³/s)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

# --- 3. Check which storms have enough isotope observations ---
# You'll need to load your isotope data here
# Expected format: DateTime, d18O_stream, d2H_stream, d18O_precip, d2H_precip

check_storm_isotopes <- function(storms, isotope_data,
                                 min_observations = 8) {
  storms %>%
    mutate(
      n_isotope_obs = map_int(1:nrow(storms), ~ {
        isotope_data %>%
          filter(DateTime >= storms$start[.x],
                 DateTime <= storms$end[.x],
                 !is.na(d18O_stream)) %>%
          nrow()
      }),
      has_enough = n_isotope_obs >= min_observations
    )
}

# --- 4. Two-component hydrograph separation ---
# Q_total  = Q_event + Q_pre_event
# C_total  = (C_event * Q_event + C_pre_event * Q_pre_event) / Q_total
# Solving: f_event = (C_stream - C_pre_event) / (C_event - C_pre_event)
# where C = isotope value (d18O or d2H), f = fraction of flow

isotope_data = read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/New%20Isotope%20Inventory%20_%20April%2013%202026.csv")
# See all column names cleanly
colnames(isotope_data)

# Check what locations and types you have
unique(isotope_data$Location)
unique(isotope_data$Type)

# Check the isotope value columns specifically
isotope_data %>%
  select(contains("Delta"), contains("d18"), contains("d2"), 
         contains("18O"), contains("2H")) %>%
  head()

# Check how many Wet Center stream samples have isotope values
isotope_data %>%
  filter(str_detect(Location, regex("wet center", ignore_case = TRUE)),
         Type == "Stream") %>%
  select(Date, Time, contains("Delta")) %>%
  head(10)

# --- 1. Clean and parse isotope data ---
isotope_clean <- isotope_data %>%
  rename(
    d2H     = `Processed Delta 2H`,
    d2H_sd  = `Processed Delta 2H StDev`,
    d18O    = `Processed Delta 18O`,
    d18O_sd = `Processed Delta 18O StDev`
  ) %>%
  mutate(
    # --- Fix known problem dates ---
    Date_clean = case_when(
      # Range dates — take the start date
      str_detect(Date, "-") & str_count(Date, "/") >= 2 ~
        str_extract(Date, "^[^-]+") %>% trimws(),
      # Single dash only date
      Date == "-" ~ NA_character_,
      TRUE ~ Date
    ),
    # --- Fix known problem times ---
    Time_clean = case_when(
      is.na(Time) | trimws(Time) %in% c("-", "w", "AM", "PM") ~ "12:00:00",
      # Remove backtick typo
      str_detect(Time, "`") ~ str_remove(Time, "`"),
      TRUE ~ as.character(Time)
    ),
    datetime_str = paste(Date_clean, Time_clean),
    # --- Parse all formats ---
    DateTime = as.POSIXct(datetime_str,
                          format = "%m/%d/%y %I:%M:%S %p",
                          tz = "America/New_York"),
    DateTime = if_else(is.na(DateTime),
                       as.POSIXct(datetime_str, format = "%m/%d/%y %I:%M %p",
                                  tz = "America/New_York"), DateTime),
    DateTime = if_else(is.na(DateTime),
                       as.POSIXct(datetime_str, format = "%m/%d/%y %H:%M:%S",
                                  tz = "America/New_York"), DateTime),
    DateTime = if_else(is.na(DateTime),
                       as.POSIXct(datetime_str, format = "%m/%d/%y %H:%M",
                                  tz = "America/New_York"), DateTime),
    DateTime = if_else(is.na(DateTime),
                       as.POSIXct(datetime_str, format = "%m/%d/%Y %H:%M:%S",
                                  tz = "America/New_York"), DateTime),
    DateTime = if_else(is.na(DateTime),
                       as.POSIXct(datetime_str, format = "%m/%d/%Y %I:%M:%S %p",
                                  tz = "America/New_York"), DateTime)
  ) %>%
  # Drop rows with no date at all or invalid dates like 2/29/25
  filter(!is.na(Date_clean), !is.na(d18O), !is.na(d2H))

# --- Check remaining failures ---
failed <- isotope_clean %>% filter(is.na(DateTime))
print(paste("Failed to parse:", nrow(failed)))
failed %>% select(Date, Time, datetime_str) %>% distinct()

# --- Check Storm 15 window now ---
isotope_clean %>%
  filter(Location == "Wet Center",
         Type %in% c("Stream", "stream", "Surface"),
         DateTime >= as.POSIXct("2025-10-10", tz = "America/New_York"),
         DateTime <= as.POSIXct("2025-10-17", tz = "America/New_York")) %>%
  select(DateTime, d18O, d2H)

# --- 2. Average replicates ---
isotope_clean <- isotope_clean %>%
  group_by(Location, DateTime, Type) %>%
  summarise(
    d2H     = mean(d2H,    na.rm = TRUE),
    d2H_sd  = mean(d2H_sd, na.rm = TRUE),
    d18O    = mean(d18O,   na.rm = TRUE),
    d18O_sd = mean(d18O_sd, na.rm = TRUE),
    .groups = "drop"
  )


# --- 3. Split into stream and precipitation endmembers ---
# Wet Center stream samples
WetCenter_isotope <- isotope_clean %>%
  filter(str_detect(Location, regex("wet center|WC ", ignore_case = TRUE)),
         Type %in% c("Stream", "stream", "Surface")) %>%
  arrange(DateTime)

# Precipitation samples
# --- Precipitation endmember: Wet Center + Woodside only ---
Precip_isotope <- isotope_clean %>%
  filter(
    Type %in% c("Precipitation", "Rain  / Ice", "Snow"),
    Location %in% c("Wet Center", "Woodside", "Elab 2", "Doolittle")
  ) %>%
  arrange(DateTime)

print(paste("Precip isotope samples:", nrow(Precip_isotope)))
range(Precip_isotope$DateTime, na.rm = TRUE)

# Sanity checks
print(paste("Wet Center isotope samples:", nrow(WetCenter_isotope)))
print(paste("Precip isotope samples:", nrow(Precip_isotope)))
range(WetCenter_isotope$DateTime, na.rm = TRUE)
range(Precip_isotope$DateTime, na.rm = TRUE)



# --- 1. Average duplicate precip samples at same DateTime ---
Precip_isotope <- Precip_isotope %>%
  group_by(DateTime) %>%
  summarise(
    d2H     = mean(d2H,     na.rm = TRUE),
    d2H_sd  = mean(d2H_sd,  na.rm = TRUE),
    d18O    = mean(d18O,    na.rm = TRUE),
    d18O_sd = mean(d18O_sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(DateTime)

# --- 2. Average duplicate stream samples at same DateTime ---
WetCenter_isotope <- WetCenter_isotope %>%
  group_by(DateTime) %>%
  summarise(
    d2H     = mean(d2H,     na.rm = TRUE),
    d2H_sd  = mean(d2H_sd,  na.rm = TRUE),
    d18O    = mean(d18O,    na.rm = TRUE),
    d18O_sd = mean(d18O_sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(DateTime)

# --- 3. Function to get closest precip sample before storm start ---
get_precip_endmember <- function(storm_start, precip_data) {
  # Find closest precip sample to storm start
  candidates <- precip_data %>%
    mutate(time_diff = abs(as.numeric(difftime(DateTime, storm_start, 
                                               units = "hours")))) %>%
    arrange(time_diff)
  
  if (nrow(candidates) == 0) {
    return(list(d18O = NA, d2H = NA, DateTime = NA, time_diff_hrs = NA))
  }
  
  best <- candidates %>% slice(1)
  
  return(list(
    d18O          = best$d18O,
    d2H           = best$d2H,
    precip_time   = best$DateTime,
    time_diff_hrs = best$time_diff
  ))
}

# --- 4. Two-component hydrograph separation function ---
# --- Updated hydrograph_separate with fallback pre-event search ---
hydrograph_separate <- function(storm_id_val, storms, flow_data,
                                stream_isotope, precip_isotope,
                                isotope = "d18O",
                                min_obs = 8) {
  
  storm <- storms %>% filter(storm_id == storm_id_val)
  
  # --- Get stream isotope samples during storm window ---
  storm_iso <- stream_isotope %>%
    filter(DateTime >= storm$start, DateTime <= storm$end)
  
  if (nrow(storm_iso) < min_obs) {
    message(paste("Storm", storm_id_val, "has only", nrow(storm_iso),
                  "isotope observations — skipping (min =", min_obs, ")"))
    return(NULL)
  }
  
  # --- Pre-event endmember ---
  # First try: 3 days before storm start
  pre_event <- stream_isotope %>%
    filter(DateTime >= storm$start - days(3),
           DateTime <  storm$start)
  
  # Fallback 1: widen to 7 days before storm start
  if (nrow(pre_event) == 0) {
    pre_event <- stream_isotope %>%
      filter(DateTime >= storm$start - days(7),
             DateTime <  storm$start)
    if (nrow(pre_event) > 0) 
      message(paste("Storm", storm_id_val, "— using 7-day pre-event window"))
  }
  
  # Fallback 2: samples on storm start day before peak
  if (nrow(pre_event) == 0) {
    pre_event <- stream_isotope %>%
      filter(DateTime >= storm$start,
             DateTime <  storm$peak_time)
    if (nrow(pre_event) > 0)
      message(paste("Storm", storm_id_val, 
                    "— using pre-peak samples as pre-event endmember"))
  }
  
  # Fallback 3: first sample in storm window
  if (nrow(pre_event) == 0) {
    pre_event <- storm_iso %>% slice(1)
    message(paste("Storm", storm_id_val,
                  "— WARNING: using first storm sample as pre-event endmember"))
  }
  
  C_pre_event <- pre_event %>%
    summarise(d18O = mean(d18O, na.rm = TRUE),
              d2H  = mean(d2H,  na.rm = TRUE))
  
  # --- Event endmember: closest precip sample to storm start ---
  precip_em <- get_precip_endmember(storm$start, precip_isotope)
  
  if (is.na(precip_em$d18O)) {
    message(paste("Storm", storm_id_val, "— no precipitation isotope data"))
    return(NULL)
  }
  
  message(paste0("Storm ", storm_id_val,
                 " | Pre-event δ¹⁸O: ", round(C_pre_event$d18O, 2),
                 " | Precip δ¹⁸O: ",    round(precip_em$d18O, 2),
                 " | Precip sample: ",   format(precip_em$precip_time, "%Y-%m-%d"),
                 " (", round(precip_em$time_diff_hrs, 1), " hrs from storm start)"))
  
  # --- Interpolate flow at isotope sample times ---
  storm_iso <- storm_iso %>%
    mutate(
      Q_total = approx(
        x    = flow_data$Time,
        y    = flow_data$Flow,
        xout = DateTime,
        rule = 2)$y
    )
  
  # --- Two-component mixing model ---
  C_pre   <- if (isotope == "d18O") C_pre_event$d18O else C_pre_event$d2H
  C_event <- if (isotope == "d18O") precip_em$d18O   else precip_em$d2H
  iso_col <- if (isotope == "d18O") "d18O"            else "d2H"
  
  storm_iso <- storm_iso %>%
    mutate(
      C_stream    = .data[[iso_col]],
      f_event     = (C_stream - C_pre)   / (C_event - C_pre),
      f_pre_event = 1 - f_event,
      f_event     = pmax(0, pmin(1, f_event)),
      f_pre_event = pmax(0, pmin(1, f_pre_event)),
      Q_event     = f_event     * Q_total,
      Q_pre_event = f_pre_event * Q_total,
      storm_id    = storm_id_val,
      isotope     = isotope,
      C_pre       = C_pre,
      C_event     = C_event
    )
  
  return(storm_iso)
}



# --- 5. Check which detected storms have enough isotope observations ---
storms_isotope_check <- storms_detected %>%
  mutate(
    n_isotope_obs = map_int(storm_id, ~ {
      WetCenter_isotope %>%
        filter(DateTime >= storms_detected$start[.x],
               DateTime <= storms_detected$end[.x]) %>%
        nrow()
    }),
    has_enough = n_isotope_obs >= 8
  )

print(storms_isotope_check %>% 
        select(storm_id, start, end, peak, n_isotope_obs, has_enough))

# --- 6. Run separation for all qualifying storms ---
qualifying_storms <- storms_isotope_check %>%
  filter(has_enough) %>%
  pull(storm_id)

print(paste("Storms with enough isotope data:", 
            paste(qualifying_storms, collapse = ", ")))

# --- Rerun separation ---
separated_d18O <- map_dfr(qualifying_storms, hydrograph_separate,
                          storms         = storms_detected,
                          flow_data      = WetCenter_Flow_Combined,
                          stream_isotope = WetCenter_isotope,
                          precip_isotope = Precip_isotope,
                          isotope        = "d18O")

separated_d2H <- map_dfr(qualifying_storms, hydrograph_separate,
                         storms         = storms_detected,
                         flow_data      = WetCenter_Flow_Combined,
                         stream_isotope = WetCenter_isotope,
                         precip_isotope = Precip_isotope,
                         isotope        = "d2H")

# --- Check all four storms now have results ---
separated_d18O %>%
  group_by(storm_id) %>%
  summarise(
    C_pre        = first(C_pre),
    C_event      = first(C_event),
    n_obs        = n(),
    mean_f_event = mean(f_event, na.rm = TRUE)
  )


# --- 7. Plot separation for each storm ---
plot_separation <- function(sep_data, storm_id_val, isotope = "d18O") {
  d <- sep_data %>% filter(storm_id == storm_id_val)
  
  if (nrow(d) == 0) {
    message("No data for storm ", storm_id_val)
    return(NULL)
  }
  
  storm_info  <- storms_detected %>% filter(storm_id == storm_id_val)
  flow_window <- WetCenter_Flow_Combined %>%
    filter(Time >= storm_info$start, Time <= storm_info$end)
  
  ggplot() +
    geom_area(data = flow_window,
              aes(x = Time, y = Flow),
              fill = "grey80", color = "grey50",
              alpha = 0.4, linewidth = 0.8) +
    # Pre-event water — bottom of bar to Q_pre_event
    geom_segment(data = d,
                 aes(x = DateTime, xend = DateTime,
                     y = 0, yend = Q_pre_event),
                 color = "royalblue", linewidth = 1.2, alpha = 0.8) +
    # Event water — Q_pre_event to Q_total
    geom_segment(data = d,
                 aes(x = DateTime, xend = DateTime,
                     y = Q_pre_event, yend = Q_total),
                 color = "steelblue", linewidth = 1.2, alpha = 0.8) +
    # Isotope line on second axis
    geom_line(data = d,
              aes(x = DateTime,
                  y = (.data[[if(isotope=="d18O") "d18O" else "d2H"]] -
                         min(.data[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) *
                    (max(flow_window$Flow, na.rm=TRUE) /
                       diff(range(.data[[if(isotope=="d18O") "d18O" else "d2H"]],
                                  na.rm=TRUE))),
                  color = "Stream Isotope"),
              linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("Stream Isotope" = "darkgreen")) +
    scale_y_continuous(
      name     = "Discharge (m³/s)",
      sec.axis = sec_axis(
        ~ . * diff(range(d[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) /
          max(flow_window$Flow, na.rm=TRUE) +
          min(d[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE),
        name = paste0(isotope, " (‰)")
      )
    ) +
    labs(title    = paste("Storm", storm_id_val, "— Hydrograph Separation (", isotope, ")"),
         subtitle = paste("C_pre =", round(d$C_pre[1], 2),
                          "‰  |  C_event =", round(d$C_event[1], 2),
                          "‰  |  Mean event fraction:",
                          round(mean(d$f_event, na.rm=TRUE), 2)),
         x = "Date/Time", color = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# --- Storm ID mapping ---
storm_labels <- tibble(
  storm_id   = c(7, 9,  10,  17,  19),
  field_name = c("Storm 10", "Storm 12", "Storm 13", "Storm 15", "Storm 16")
)

# Then in plot_separation, replace the title line with:
field_name <- storm_labels %>%
  filter(storm_id == storm_id_val) %>%
  pull(field_name)

print(unique(separated_d18O$storm_id))

storms_detected %>%
  select(storm_id, start, peak_time, end, peak) %>%
  arrange(start) %>%
  print(n = 27)

# --- Plot all ---
plot_separation <- function(sep_data, storm_id_val, isotope = "d18O") {
  d <- sep_data %>% filter(storm_id == storm_id_val)
  
  if (nrow(d) == 0) {
    message("No data for storm ", storm_id_val)
    return(NULL)
  }
  
  storm_info  <- storms_detected %>% filter(storm_id == storm_id_val)
  flow_window <- WetCenter_Flow_Combined %>%
    filter(Time >= storm_info$start, Time <= storm_info$end)
  
  # --- Interpolate separation components to 10-min resolution ---
  separation_interp <- flow_window %>%
    mutate(
      Q_pre_event = approx(
        x    = d$DateTime,
        y    = d$Q_pre_event,
        xout = Time,
        rule = 2)$y,
      f_pre_event = approx(
        x    = d$DateTime,
        y    = d$f_pre_event,
        xout = Time,
        rule = 2)$y
    ) %>%
    mutate(
      # Clamp fraction to 0-1
      f_pre_event  = pmax(0, pmin(1, f_pre_event)),
      # Use actual flow as total, split by interpolated fraction
      Q_pre_event  = f_pre_event * Flow,
      Q_event      = (1 - f_pre_event) * Flow
    )
  field_name <- storm_labels %>%
    filter(storm_id == storm_id_val) %>%
    pull(field_name)
  ggplot() +
    # --- Pre-event water: bottom ribbon ---
    geom_ribbon(data = separation_interp,
                aes(x = Time, ymin = 0, ymax = Q_pre_event),
                fill = "steelblue", alpha = 0.6) +
    # --- Event water: top ribbon ---
    geom_ribbon(data = separation_interp,
                aes(x = Time, ymin = Q_pre_event, ymax = Flow),
                fill = "black", alpha = 0.6) +
    # --- Total flow line on top for reference ---
    geom_line(data = flow_window,
              aes(x = Time, y = Flow),
              color = "steelblue", linewidth = 0.8) +
    # --- Isotope points at sample times ---
    geom_line(data = d,
              aes(x = DateTime,
                  y = (.data[[if(isotope=="d18O") "d18O" else "d2H"]] -
                         min(.data[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) *
                    (max(flow_window$Flow, na.rm=TRUE) /
                       diff(range(.data[[if(isotope=="d18O") "d18O" else "d2H"]],
                                  na.rm=TRUE))),
                  color = "Stream Isotope"),
              linewidth = 2, linetype = "dotdash") +
    scale_color_manual(values = c("Stream Isotope" = "magenta")) +
    scale_y_continuous(
      name     = "Discharge (m³/s)",
      sec.axis = sec_axis(
        ~ . * diff(range(d[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) /
          max(flow_window$Flow, na.rm=TRUE) +
          min(d[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE),
        name = paste0(isotope, " (‰)")
      )
    ) +
    labs(title = paste(field_name, "— Hydrograph Separation (", isotope, ")"),
         subtitle = paste("C_pre =", round(d$C_pre[1], 2),
                          "‰  |  C_event =", round(d$C_event[1], 2),
                          "‰  |  Mean event fraction:",
                          round(mean(d$f_event, na.rm=TRUE), 2)),
         x = "Date/Time", color = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}

library(patchwork)

plots <- map(unique(separated_d18O$storm_id),
             ~ plot_separation(separated_d18O, .x, "d18O"))

wrap_plots(plots, ncol = 3) +
  plot_annotation(title = "Hydrograph Separation — All Storms (δ¹⁸O)")



# --- 9. Compare d18O vs d2H separation as uncertainty check ---
if (nrow(separated_d18O) > 0 & nrow(separated_d2H) > 0) {
  comparison <- separated_d18O %>%
    select(DateTime, storm_id, f_event_d18O = f_event) %>%
    left_join(separated_d2H %>% select(DateTime, storm_id, f_event_d2H = f_event),
              by = c("DateTime", "storm_id"))
  
  ggplot(comparison, aes(x = f_event_d18O, y = f_event_d2H, 
                         color = factor(storm_id))) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "δ¹⁸O vs δ²H Separation Comparison",
         subtitle = "Points on 1:1 line = consistent separation between tracers",
         x = "Event fraction (δ¹⁸O)", y = "Event fraction (δ²H)",
         color = "Storm ID") +
    theme_bw()
}



# --- Check what isotope samples fall in the Storm 15 field window ---
WetCenter_isotope %>%
  filter(DateTime >= as.POSIXct("2025-10-10", tz = "America/New_York"),
         DateTime <= as.POSIXct("2025-10-17", tz = "America/New_York")) %>%
  select(DateTime, d18O, d2H)

# --- Check what the detected storm 15 window looks like ---
storms_detected %>% filter(storm_id == 15)

# --- Also check storm 16 detected window ---
storms_detected %>% filter(storm_id == 16)





library(dygraphs)
library(xts)

# --- Convert to xts format which dygraphs requires ---
flow_xts <- WetCenter_Flow_Combined %>%
  arrange(Time) %>%
  select(Time, Flow) %>%
  as.data.frame() %>%
  { xts(.$Flow, order.by = .$Time) }

# --- Basic interactive hydrograph ---
dygraph(flow_xts, main = "Mill River — Combined Flow Record") %>%
  dyAxis("y", label = "Discharge (m³/s)") %>%
  dyRangeSelector() %>%        # slider at bottom to zoom/pan
  dyOptions(fillGraph = TRUE, fillAlpha = 0.2, colors = "steelblue") %>%
  dyHighlight(highlightSeriesOpts = list(strokeWidth = 2))

# --- Add PFAS samples as events ---
# --- Create a sparse xts series for sample points ---
# Sets Flow value at sample times, NA everywhere else
sample_xts <- WetCenter_Flow_Combined %>%
  arrange(Time) %>%
  mutate(SampleFlow = ifelse(Time %in% 
                               as.POSIXct(WetCenter_samples$DateTime, tz = "America/New_York"),
                             Flow, NA)) %>%
  select(Time, Flow, SampleFlow) %>%
  { xts(.[, c("Flow", "SampleFlow")], order.by = .$Time) }

p <- dygraph(sample_xts, main = "Mill River — Combined Flow Record") %>%
  dyAxis("y", label = "Discharge (m³/s)") %>%
  dyOptions(fillAlpha = 0.2) %>%
  dySeries("Flow",       label = "Discharge",   color = "steelblue",
           fillGraph = TRUE) %>%
  dySeries("SampleFlow", label = "PFAS Sample", color = "red",
           drawPoints = TRUE, pointSize = 5, strokeWidth = 0) %>%
  dyRangeSelector() %>%
  dyHighlight(highlightSeriesOpts = list(strokeWidth = 2))

# --- Add storm shading ---
for (i in seq_len(nrow(storms_detected))) {
  p <- p %>%
    dyShading(from = format(storms_detected$start[i], "%Y-%m-%dT%H:%M:%S"),
              to   = format(storms_detected$end[i],   "%Y-%m-%dT%H:%M:%S"),
              color = "#ADD8E650")
}
p

