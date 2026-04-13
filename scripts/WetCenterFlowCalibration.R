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
WetCenterDischarge <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenterDischarge.csv")
AllChem_WetCenter_SE15 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Storm_Event_15_all_results.csv")

# --- Reparse WetCenterDischarge$Time ---
WetCenterDischarge <- WetCenterDischarge %>%
  mutate(Time = mdy_hm(Time, tz = "America/New_York"))

# --- Recompute WetCenterDischarge flow with new NLS coefficients ---
WetCenterDischarge <- WetCenterDischarge %>%
  mutate(
    GaugeHeight_m = `Adjusted Gauge Height` * 0.3048,
    Flow          = a * (GaugeHeight_m ^ b)  # same NLS coefficients from rating curve
  )

# --- Trim WetCenterDischarge to end where pressure transducer begins ---
cutover <- min(WetStation4Stream$Time, na.rm = TRUE)

WetCenterDischarge_trimmed <- WetCenterDischarge %>%
  filter(Time < cutover) %>%
  transmute(
    Time,
    GaugeHeight_m,
    Flow,
    Source = "Staff Gauge"
  )

# --- Prepare WetStation4Stream with matching columns and remove first entry ---
WetStation4Stream_labeled <- WetStation4Stream %>%
  slice(-1) %>%  # removes first row
  transmute(
    Time,
    GaugeHeight_m,
    Flow,
    Source = "Pressure Transducer"
  )
# --- Combine ---
WetCenter_Flow_Combined <- bind_rows(WetCenterDischarge_trimmed,
                                     WetStation4Stream_labeled) %>%
  arrange(Time)

# --- Sanity checks ---
range(WetCenter_Flow_Combined$Time, na.rm = TRUE)
table(WetCenter_Flow_Combined$Source)

# --- Plot to visually inspect the join ---
ggplot(WetCenter_Flow_Combined, aes(x = Time, y = Flow, color = Source)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Staff Gauge" = "darkorange",
                                "Pressure Transducer" = "steelblue")) +
  labs(title = "Mill River — Combined Flow Record",
       x = "Date/Time", y = "Discharge (m³/s)") +
  theme_bw() +
  theme(legend.position = "bottom")

# so the time view data is underpredicting flow
overlap_start <- min(WetStation4Stream$Time, na.rm = TRUE)
overlap_end   <- as.POSIXct("2025-10-20 21:15:00", tz = "America/New_York")

tv_overlap <- WetCenterDischarge %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  summarise(mean_H = mean(GaugeHeight_m, na.rm = TRUE))

pt_overlap <- WetStation4Stream %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  summarise(mean_H = mean(GaugeHeight_m, na.rm = TRUE))

# The offset we need to add to TimeView to match pressure transducer
datum_offset <- pt_overlap$mean_H - tv_overlap$mean_H
print(paste("Datum offset (m):", round(datum_offset, 4)))


# --- Remove outliers and apply datum correction ---
WetCenterDischarge_clean <- WetCenterDischarge %>%
  filter(GaugeHeight_m < 0.65) %>%
  mutate(
    GaugeHeight_m = GaugeHeight_m + datum_offset,
    Flow          = a * (GaugeHeight_m ^ b)
  )

# --- Verify ---
summary(WetCenterDischarge_clean$GaugeHeight_m)
summary(WetCenterDischarge_clean$Flow)

# --- Rebuild combined dataset ---
WetCenterDischarge_trimmed <- WetCenterDischarge_clean %>%
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

# --- Plot to inspect ---
ggplot(WetCenter_Flow_Combined, aes(x = Time, y = Flow, color = Source)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue")) +
  labs(title = "Mill River — Combined Flow Record",
       x = "Date/Time", y = "Discharge (m³/s)") +
  theme_bw() +
  theme(legend.position = "bottom")



# --- Plot ONLY the overlap period to assess the correction ---
ggplot() +
  geom_line(data = WetCenterDischarge_clean %>% 
              filter(Time >= overlap_start, Time <= overlap_end),
            aes(x = Time, y = Flow, color = "TimeView"),
            linewidth = 0.8) +
  geom_line(data = WetStation4Stream_labeled %>%
              filter(Time >= overlap_start, Time <= overlap_end),
            aes(x = Time, y = Flow, color = "Pressure Transducer"),
            linewidth = 0.8) +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue")) +
  labs(title = "Overlap Period — Datum Correction Check",
       x = "Date/Time", y = "Discharge (m³/s)") +
  theme_bw() +
  theme(legend.position = "bottom")


# --- Compute a regression-based correction instead of simple offset ---
# Join the two sources by nearest timestamp during overlap
library(data.table)

tv_dt <- WetCenterDischarge_clean %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  select(Time, GaugeHeight_m) %>%
  as.data.table()

pt_dt <- WetStation4Stream_labeled %>%
  filter(Time >= overlap_start, Time <= overlap_end) %>%
  select(Time, GaugeHeight_m) %>%
  rename(H_pt = GaugeHeight_m) %>%
  as.data.table()

overlap_matched <- tv_dt[pt_dt, on = "Time", roll = "nearest"] %>%
  as_tibble() %>%
  rename(H_tv = GaugeHeight_m)

# Fit a linear model: PT height ~ TimeView height
overlap_model <- lm(H_pt ~ H_tv, data = overlap_matched)
summary(overlap_model)

# Extract corrected coefficients
intercept    <- coef(overlap_model)["(Intercept)"]
slope        <- coef(overlap_model)["H_tv"]

print(paste("Correction: H_corrected =", round(slope, 4), "* H_tv +", round(intercept, 4)))



# --- Apply regression-based correction to TimeView data ---
WetCenterDischarge_corrected <- WetCenterDischarge_clean %>%
  mutate(
    GaugeHeight_m = slope * GaugeHeight_m + intercept,
    Flow          = a * (GaugeHeight_m ^ b)
  )

# --- Check the overlap now ---
ggplot() +
  geom_line(data = WetCenterDischarge_corrected %>%
              filter(Time >= overlap_start, Time <= overlap_end),
            aes(x = Time, y = Flow, color = "TimeView"),
            linewidth = 0.8) +
  geom_line(data = WetStation4Stream_labeled %>%
              filter(Time >= overlap_start, Time <= overlap_end),
            aes(x = Time, y = Flow, color = "Pressure Transducer"),
            linewidth = 0.8) +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue")) +
  labs(title = "Overlap Period — Regression Correction Check",
       x = "Date/Time", y = "Discharge (m³/s)") +
  theme_bw() +
  theme(legend.position = "bottom")


#Alot better!!

# --- Rebuild combined dataset with corrected TimeView ---
WetCenterDischarge_trimmed <- WetCenterDischarge_corrected %>%
  filter(Time < overlap_start) %>%
  transmute(Time, GaugeHeight_m, Flow, Source = "TimeView")

WetCenter_Flow_Combined <- bind_rows(WetCenterDischarge_trimmed,
                                     WetStation4Stream_labeled) %>%
  arrange(Time)

# --- Full record plot ---
ggplot(WetCenter_Flow_Combined, aes(x = Time, y = Flow, color = Source)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("TimeView" = "darkorange",
                                "Pressure Transducer" = "steelblue")) +
  labs(title = "Mill River — Combined Flow Record",
       x = "Date/Time", y = "Discharge (m³/s)") +
  theme_bw() +
  theme(legend.position = "bottom")



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
  rise_factor      = 1.4,
  min_duration_hrs = 6,
  pre_event_days   = 3,
  recession_factor = 1.1,
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



hydrograph_separate <- function(storm_id_val, storms, flow_data,
                                isotope_data, isotope = "d18O") {
  
  storm      <- storms %>% filter(storm_id == storm_id_val)
  iso_stream <- paste0(isotope, "_stream")
  iso_precip <- paste0(isotope, "_precip")
  
  # Get storm isotope observations
  storm_iso <- isotope_data %>%
    filter(DateTime >= storm$start, DateTime <= storm$end,
           !is.na(.data[[iso_stream]]))
  
  if (nrow(storm_iso) < 2) {
    message("Not enough isotope data for storm ", storm_id_val)
    return(NULL)
  }
  
  # Pre-event endmember: mean isotope value of stream before storm
  C_pre <- isotope_data %>%
    filter(DateTime >= storm$start - days(3),
           DateTime < storm$start) %>%
    summarise(m = mean(.data[[iso_stream]], na.rm = TRUE)) %>%
    pull(m)
  
  # Event endmember: precipitation isotope value
  C_event <- isotope_data %>%
    filter(DateTime >= storm$start, DateTime <= storm$end) %>%
    summarise(m = mean(.data[[iso_precip]], na.rm = TRUE)) %>%
    pull(m)
  
  # Interpolate total flow at isotope sample times
  storm_iso <- storm_iso %>%
    mutate(
      Q_total   = approx(x = flow_data$Time, y = flow_data$Flow,
                         xout = DateTime)$y,
      # Mixing model
      f_event   = (.data[[iso_stream]] - C_pre) / (C_event - C_pre),
      f_preevent = 1 - f_event,
      Q_event   = f_event * Q_total,
      Q_preevent = f_preevent * Q_total,
      storm_id  = storm_id_val,
      isotope   = isotope,
      C_pre     = C_pre,
      C_event   = C_event
    )
  
  return(storm_iso)
}

# --- 5. Run separation for all storms with enough data ---
# Once isotope data is loaded, run like this:
# separated <- storms_with_enough %>%
#   pull(storm_id) %>%
#   map_dfr(hydrograph_separate,
#           storms     = storms_detected,
#           flow_data  = WetCenter_Flow_Combined,
#           isotope_data = your_isotope_df,
#           isotope    = "d18O")

# --- 6. Plot separation results for a single storm ---
plot_separation <- function(separation_data, storm_id_val) {
  d <- separation_data %>% filter(storm_id == storm_id_val)
  
  ggplot(d, aes(x = DateTime)) +
    geom_area(aes(y = Q_total, fill = "Total Flow"),
              alpha = 0.3) +
    geom_area(aes(y = Q_preevent, fill = "Pre-event Flow"),
              alpha = 0.6) +
    geom_area(aes(y = Q_event, fill = "Event Flow"),
              alpha = 0.6) +
    geom_point(aes(y = Q_total), size = 2) +
    scale_fill_manual(values = c("Total Flow"    = "steelblue",
                                 "Pre-event Flow" = "brown3",
                                 "Event Flow"    = "skyblue")) +
    labs(title = paste("Hydrograph Separation — Storm", storm_id_val),
         subtitle = paste("Endmembers: C_pre =", round(d$C_pre[1], 2),
                          "| C_event =", round(d$C_event[1], 2)),
         x = "Date/Time", y = "Discharge (m³/s)", fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}


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

