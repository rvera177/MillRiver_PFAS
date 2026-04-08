
library(tidyverse)
library(lubridate)
library(readr)

# 1. Load the three datasets
water_abs_df <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation4_WaterPressure.csv")
air_abs_df   <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation3_AirPressure.csv")
gauge_df     <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenter%20GaugeHeight%20-%20Sheet1-2.csv")

# 2. Process Water Pressure (Daily Average)
water_daily <- water_abs_df %>%
  mutate(Date = mdy(word(`Date-Time (EDT)`, 1))) %>%
  group_by(Date) %>%
  summarise(Avg_Water_Abs_kPa = mean(`Absolute Pressure kPa`, na.rm = TRUE))

# 3. Process Air Pressure (Daily Average)
air_daily <- air_abs_df %>%
  mutate(Date = mdy(word(`Date-Time (EDT)`, 1))) %>%
  group_by(Date) %>%
  summarise(Avg_Air_kPa = mean(`Absolute Pressure kPa`, na.rm = TRUE))

# 4. Process Gauge Height (Ground Truth)
gauge_clean <- gauge_df %>%
  filter(!is.na(`GaugeHeight (ft)`)) %>%
  mutate(Date = mdy(Date))

# 5. Combine and Calculate "Water Only" Pressure
# Formula: Water Pressure = Absolute Pressure - Air Pressure
calib_points <- gauge_clean %>%
  inner_join(water_daily, by = "Date") %>%
  inner_join(air_daily, by = "Date") %>%
  mutate(Water_Pressure_Only = Avg_Water_Abs_kPa - Avg_Air_kPa)

# 6. Create the Calibration Model
# This relates your manual ft readings to the compensated kPa
calib_model <- lm(`GaugeHeight (ft)` ~ Water_Pressure_Only, data = calib_points)

# View model accuracy (look for R-squared)
summary(calib_model)

# 7. Apply to the entire high-frequency dataset
# First, we match the high-freq water data with daily air pressure for compensation
final_dataset <- water_abs_df %>%
  mutate(Date = mdy(word(`Date-Time (EDT)`, 1))) %>%
  left_join(air_daily, by = "Date") %>%
  mutate(
    Compensated_Pressure_kPa = `Absolute Pressure kPa` - Avg_Air_kPa,
    # Apply the slope and intercept from our model
    Calibrated_Height_ft = predict(calib_model, newdata = data.frame(Water_Pressure_Only = Compensated_Pressure_kPa))
  )

# 8. Plot the Calibration Curve
ggplot(calib_points, aes(x = Water_Pressure_Only, y = `GaugeHeight (ft)`)) +
  geom_point(color = "darkgreen", size = 2) +
  geom_smooth(method = "lm", color = "black") +
  labs(title = "Barometric Compensated Calibration",
       x = "Water Pressure (Absolute - Air) [kPa]",
       y = "Manual Gauge Height (ft)") +
  theme_bw()



# 1. Load libraries
library(tidyverse)
library(lubridate)
library(readr)

# 1. Load the three datasets
water_abs_df <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation4_WaterPressure.csv")
air_abs_df   <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WETStation3_AirPressure.csv")
gauge_df     <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenter%20GaugeHeight%20-%20Sheet1-2.csv")

# 3. Clean column names (removes hidden spaces)
colnames(water_abs_df) <- trimws(colnames(water_abs_df))
colnames(air_abs_df)   <- trimws(colnames(air_abs_df))
colnames(gauge_df)     <- trimws(colnames(gauge_df))

# 4. Convert all data to a pure "Date" format (No Time)
# We use parse_date_time because it handles 2-digit years (25) and 4-digit years (2025)
water_daily <- water_abs_df %>%
  mutate(Date = as.Date(parse_date_time(`Date-Time (EDT)`, orders = "mdy HMS"))) %>%
  group_by(Date) %>%
  summarise(Avg_Water_Abs = mean(`Absolute Pressure kPa`, na.rm = TRUE))

air_daily <- air_abs_df %>%
  mutate(Date = as.Date(parse_date_time(`Date-Time (EDT)`, orders = "mdy HMS"))) %>%
  group_by(Date) %>%
  summarise(Avg_Air = mean(`Absolute Pressure kPa`, na.rm = TRUE))

gauge_clean <- gauge_df %>%
  filter(!is.na(`GaugeHeight (ft)`)) %>%
  mutate(Date = as.Date(parse_date_time(Date, orders = "mdy")))

# 5. Join to create Daily Compensated Pressure
# Water Pressure = Absolute Pressure - Air Pressure
daily_comp <- inner_join(water_daily, air_daily, by = "Date") %>%
  mutate(Water_Pressure_Only = Avg_Water_Abs - Avg_Air)

# 6. Create Calibration Points by joining Gauge to Daily Data
calib_points <- inner_join(gauge_clean, daily_comp, by = "Date")

print(paste("Matches found:", nrow(calib_points)))

# 7. Build the Model
calib_model <- lm(`GaugeHeight (ft)` ~ Water_Pressure_Only, data = calib_points)
summary(calib_model)

# 8. Apply the model to the Daily Sensor Data
daily_results <- daily_comp %>%
  mutate(
    Calibrated_Height_ft = predict(calib_model, newdata = data.frame(Water_Pressure_Only = Water_Pressure_Only))
  )

# 9. Graph by Date Only
ggplot() +
  # The "Sensor Line" (Daily Averages)
  geom_line(data = daily_results, aes(x = Date, y = Calibrated_Height_ft, color = "Daily Average (Sensor)"), size = 1) +
  # The "Manual Points"
  geom_point(data = calib_points, aes(x = Date, y = `GaugeHeight (ft)`, color = "Manual Reading"), size = 3) +
  scale_color_manual(values = c("Daily Average (Sensor)" = "steelblue", "Manual Reading" = "red")) +
  labs(title = "Daily Water Level Calibration (Date Only)",
       subtitle = "Time of day ignored; matches based on daily averages",
       x = "Date",
       y = "Water Level (ft)") +
  theme_minimal() +
  theme(legend.position = "bottom")

