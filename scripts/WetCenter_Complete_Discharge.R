

#Mill River PFAS Project
#Created November 28, 2025
#Updated March 29, 2025
#Raul Vera


#Graphing results for storm 15 and 16
#storm 15: October 12-13. 0.7 inches of rainfall.
          #samples from every other isco bottle for 12 total samples

#Storm 16: October 30-31. 2 inches of rainfall
          #samples from every isco bottle for 24 total samples

library(readr) #reads csv from github
library(lubridate) # this is for making the date and time columns readable for time-series plots
library(ggplot2) #for plots
library(dplyr) #for pipes and renaming columns
library(tidyr) # for pivot_longer

library(scales)
library(viridis)  # for a specific color. This can be switch for wes anderson
library(zoo) # for na.approx
library(purrr) #i dont remember

#-----Graphing storm 16------------------



#now for storm 16, which was the Halloween storm. 
# I don't have stage height. I have pressure data. 
#this needs to be converted to stage height,
#and then I'll use rating curve equation for flow


#----convert pressure data to discharge------
WetStation3Air <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetStation3Air.csv")
WetStation4Stream <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetStation4Stream.csv")

#  Wet Stations are in spanish jaja
#translate from spanish to english using rename. 

WetStation3Air <- WetStation3Air %>%
  rename(Time = `Fecha/hora (EDT/EST)`,
    Pressure_kPa = `Presión absoluta , kPa`, #changing to Pa in a bit
    Temp_C = `Temperatura , °C`
  ) %>%
  mutate(Time = mdy_hms(Time, tz = "America/New_York"),
    Pressure_Pa = Pressure_kPa * 1000) %>% # converting kPa to Pa
  select(Time, Pressure_Pa, Temp_C) 

WetStation4Stream <- WetStation4Stream %>%
  rename(Time = `Fecha/hora (EDT/EST)`,
    Pressure_kPa = `Presión absoluta , kPa`,
    Temp_C = `Temperatura , °C`) %>%
  mutate(Time = mdy_hms(Time, tz = "America/New_York"),
    Pressure_Pa = Pressure_kPa * 1000) %>% # converting kPa to Pa
  select(Time, Pressure_Pa, Temp_C) 

#known stage depth used for calibrating
cal_time <- as.POSIXct("2025-10-24 10:07")
cal_depth <- 0.366  #in meters

air_interp <- approx( #lining up air data with stream data chronologically
  x = WetStation3Air$Time,
  y = WetStation3Air$Pressure_Pa,
  xout = WetStation4Stream$Time)$y

WetStation4Stream$AirPressure <- air_interp #adds air pressure to the stream dataset
#next, the difference between the two pressures.
WetStation4Stream$Pressure_diff <- WetStation4Stream$Pressure_Pa - WetStation4Stream$AirPressure

P_diff_cal <- approx( #the difference in pressure during calibration time
  x = WetStation4Stream$Time,
  y = WetStation4Stream$Pressure_diff,
  xout = cal_time)$y

k <- cal_depth / (P_diff_cal / (1000 * 9.81))  # ρ = 1000 kg/m3, g = 9.81 m/s²
# k is the calibration of the pressure tranducer
#the transducer says it's a certain depth, but I measured something different. 
# I'm going with my measurement to be actual depth. K is 1.06. so not much of a difference anyway

WetStation4Stream$GaugeHeight <- (WetStation4Stream$Pressure_diff / (1000 * 9.81)) * k
WetStation4Stream$Flow <- 2.74 * (WetStation4Stream$GaugeHeight^1.73)

#bring in data from github. Not plotting precip on storm 15 at the moment. I will for AGU.
#PrecipitationData_Storms15 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PrecipitationData_Storms15.csv")

PrecipitationData_Storms16 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PrecipitationData_Storms16.csv")
PFAS_Inventory_temp <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/PFAS%20Inventory%20-%20temp.csv")
Storm16_PFAS <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/RV_PFAS_Storm16.csv")
#---mutating the DateTime format so it's consistent.

PrecipitationData_Storms16 <- PrecipitationData_Storms16 %>%
  mutate(DateTime = mdy_hm(DateTime, tz = "America/New_York"))

PFAS_Inventory_temp <- PFAS_Inventory_temp %>%
  mutate(DateTime = mdy_hm(DateTime, tz = "America/New_York"))


# joining sample times to WetStation4Stream to get corresponding Flow during sample time
PFAS_points <- PFAS_Inventory_temp %>%
  mutate(Flow = approx(
      x = WetStation4Stream$Time,
      y = WetStation4Stream$Flow,
      xout = DateTime)$y) 
PFAS_points <- PFAS_points %>%
  mutate(Sample = row_number())

max_flow <- 3.5  # Maximum discharge on the y axis for positioning the rainfall bars
max_precip <- 1  # Maximum rainfall on the y axis for scaling against flow
scale_factor <- max_flow / max_precip  #scaling again

#plot for storm 16.
ggplot() +
  geom_line(data = WetStation4Stream,
    aes(x = Time, y = Flow), color = "blue3",linewidth = 2.0) +
  geom_rect(
    data = PrecipitationData_Storms16, aes(
      xmin = DateTime - seconds(7.5*60), # widening rain bars for visibility
      xmax = DateTime + seconds(7.5*60),  
      ymin = max_flow - `Precip (in)` * scale_factor,
      ymax = max_flow),fill = "steelblue",alpha = 1) +
  geom_point(data = PFAS_points,
    aes(x = DateTime, y = Flow,), color = "red",size = 6) +
  geom_label(
    data = PFAS_points, aes(x = DateTime, y = Flow, label = Sample),
    vjust = -0.7, # positioning above the point
    fill = "white", # background label color
    color = "red",  # label text colour
    size = 4, label.size = 0.2, # border thickness around label
    fontface = "bold", alpha = 1)+
  scale_y_continuous(   # Y axis: left = Flow, right = Rainfall
    name = "Discharge (m³/s)",sec.axis = sec_axis(~ (max_flow - .) #left side
      / scale_factor, name = "Rainfall (in/hr)"), #right side
    expand = c(0,0)) + #expand(0,0) removes the whitespace buffer around the plots
  scale_x_datetime(
    limits = as.POSIXct(c("2025-10-29 20:00:00", "2025-11-07 18:00:00")),
    date_labels = "%b %d",
    date_breaks = "1 day", expand = c(0,0)) +
  labs(title = "Storm 16 Hydrograph (2 inches of total rainfall)", x = "Time") +
  
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 13))

ggplot() +
  geom_line(data = WetStation4Stream,
            aes(x = Time, y = Flow), color = "blue3",linewidth = 2.0)


library(dplyr)
library(lubridate)
library(highcharter)
library(purrr)

# date range (user-specified)
range_start <- mdy_hm("10/29/2025 23:59", tz = "EST")
range_end   <- mdy_hm("11/04/2025 23:59", tz = "EST")



WetStream2 <- WetStation4Stream %>%
  mutate(Time = as.POSIXct(Time, tz = "EST"))

# 2) Filter both datasets to the requested window
Storm16_PFAS_filt <- Storm16_PFAS %>%
  filter(DateTime >= range_start & DateTime <= range_end)

WetStream_filt <- WetStream2 %>%
  filter(Time >= range_start & Time <= range_end)

# 3) Interpolate Flow to PFAS sample times (so tooltips can show interpolated flow if desired)
if(nrow(WetStream_filt) >= 2 && nrow(Storm16_PFAS_filt) > 0){
  Storm16_PFAS_filt <- Storm16_PFAS_filt %>%
    mutate(
      Flow_interp = approx(
        x = as.numeric(WetStream_filt$Time),
        y = WetStream_filt$Flow,
        xout = as.numeric(DateTime),
        rule = 2
      )$y
    )
} else {
  Storm16_PFAS_filt <- Storm16_PFAS_filt %>% mutate(Flow_interp = NA_real_)
}

# 4) Helper to create axis min/max with small padding (avoids identical min==max)
axis_limits <- function(x){
  x <- x[!is.na(x)]
  if(length(x) == 0) return(c(0,1))
  mn <- min(x); mx <- max(x)
  if (mn == mx) {
    pad <- ifelse(mn == 0, 1, abs(mn) * 0.1)
    return(c(mn - pad, mx + pad))
  }
  pad <- (mx - mn) * 0.05
  c(mn - pad, mx + pad)
}

flow_lim   <- axis_limits(WetStream_filt$Flow)
pfoa_lim   <- axis_limits(Storm16_PFAS_filt$PFOA)
pfhxa_lim  <- axis_limits(Storm16_PFAS_filt$PFHxA)

# 5) Prepare series data for highcharter (datetime in ms)
flow_series <- WetStream_filt %>%
  transmute(x = as.numeric(Time) * 1000, y = Flow) %>%
  list_parse2()

pfoa_series <- Storm16_PFAS_filt %>%
  transmute(x = as.numeric(DateTime) * 1000, y = PFOA, sample = Sample, flow = round(Flow_interp,4)) %>%
  list_parse2()

pfhxa_series <- Storm16_PFAS_filt %>%
  transmute(x = as.numeric(DateTime) * 1000, y = PFHxA, sample = Sample, flow = round(Flow_interp,4)) %>%
  list_parse2()

# 6) Build the highchart (mimicking the style you posted)
hc <- highchart() %>%
  hc_title(text = "Storm 16: PFOA and Flow") %>%
  hc_xAxis(type = "datetime",
           min = as.numeric(range_start) * 1000,
           max = as.numeric(range_end) * 1000) %>%
  hc_yAxis_multiples(
    # y-axis 0: Flow (background area)
    list(title = list(text = "Flow (m3/s)"),
         lineWidth = 3, lineColor = "blue",
         min = flow_lim[1], max = flow_lim[2]),
    # y-axis 1: PFOA
    list(title = list(text = "PFOA (ng/L)"),
         lineWidth = 3, lineColor = "red",
         min = pfoa_lim[1], max = pfoa_lim[2])
    # y-axis 2: PFHxA
#    list(title = list(text = "PFHxA (ng/L)"),
#         lineWidth = 3, lineColor = "purple",
#         min = pfhxa_lim[1], max = pfhxa_lim[2])
  ) %>%
  # Flow as area in the background (yAxis = 0)
  hc_add_series(name = "Flow",
                data = flow_series,
                yAxis = 0,
                type = "area",
                lineWidth = 5,
                color = "blue",
                fillOpacity = 0.25,
                marker = list(enabled = FALSE)
  ) %>%
  hc_plotOptions(series = list(connectNulls = TRUE)) %>%
  # PFOA as thick line with markers (yAxis = 1)
  hc_add_series(name = "PFOA",
                data = pfoa_series,
                yAxis = 1,
                type = "scatter",
                lineWidth = 5,
                color = "red",
                marker = list(symbol = "circle", radius = 6)) %>%
  # PFHxA as thick line with markers (yAxis = 2)
 # hc_add_series(name = "PFHxA",
#                data = pfhxa_series,
#                yAxis = 2,
#                type = "scatter",
#               lineWidth = 5,
#                color = "purple",
#                marker = list(symbol = "diamond", radius = 6)
#  ) %>%
  hc_tooltip(shared = FALSE, crosshairs = TRUE,
             headerFormat = "<b>{point.key}</b><br/>") %>%
  hc_legend(enabled = TRUE)

# show
hc
#-------Hysteresis plot!----------
library(dplyr)
library(lubridate)
library(ggplot2)
library(viridis)

#make sure times are in POSIXct (it should be already, but just in case)
Storm15 <- Storm15 %>%
  mutate(Time = as.POSIXct(Time)) %>%
  arrange(Time)

#subset my datasets and mutate them 
PFAS <- PFAS_WetCenter_2025 %>%
  filter(Type == "Stream") %>%       #only stream samples, no precip
  mutate(DateTime = mdy_hm(DateTime)) %>%   # convert timing ex. 7/16/25 05:00
  arrange(DateTime)

#filtering PFAS samples inside my objects
storm_start <- min(Storm15$Time)
storm_end   <- max(Storm15$Time)

PFAS_storm15 <- PFAS %>%
  filter(DateTime >= storm_start - hours(6),   # buffer optional
         DateTime <= storm_end + hours(6))

#combine PFAS concentration to stream flow using left_join
PFAS_joined <- PFAS_storm15 %>%
  mutate(
    Time_Flow = Storm15$Time[findInterval(DateTime, Storm15$Time)]
  ) %>%
  left_join(Storm15, by = c("Time_Flow" = "Time"))

#pfas compounds that I want to make hysteresis loops for.
# there's probably a better way of making the list
pfas_cols <- c("PFOS Results",
               "PFOA Results",
               "PFHxS Results",
               "PFHxA Results",
               "PFNA Results",
               "PFBA Results",
               "PFNS Results",
               "PFBS Results", 
               "8:2FTS Results",
               "ADONA Results")

# then it's the hystersis and flushing index calculation function
#this normalizes flow and concentrations
calc_hysteresis <- function(flow, conc) {
  
  # Normalizing
  nq <- (flow - min(flow)) / (max(flow) - min(flow))
  nc <- (conc - min(conc)) / (max(conc) - min(conc))
  
  # Peak discharge index
  peak_idx <- which.max(flow)
  
  rising <- nc[1:peak_idx]
  falling <- nc[(peak_idx+1):length(nc)]
  
  L <- min(length(rising), length(falling))
  
  rising_resamp  <- approx(seq_along(rising), rising, seq(1, length(rising), length.out=L))$y
  falling_resamp <- approx(seq_along(falling), falling, seq(1, length(falling), length.out=L))$y
  
  list(
    nq = nq,
    nc = nc,
    HI = rising_resamp - falling_resamp,
    FI = nc[peak_idx] - nc[1])
}

out_dir <- "PFAS_Hysteresis_Storm15"
dir.create(out_dir, showWarnings = FALSE)
#for loop through each species you added to your compounds list.
for (compound in pfas_cols) {
  df <- PFAS_joined %>%
    select(DateTime, `Flow (m3/s)`, conc = all_of(compound)) %>%
    filter(!is.na(conc))
  #this is for a case where I don't have enough data across a storm.
  #hysteresis datapoints should be maximum 8 hours from each other
  # ^ that's according to a talk from AGU2024
  if (nrow(df) < 4) {
    message("Skipping ", compound, " (too few samples)")
    next
  }
  
  H <- calc_hysteresis(df$`Flow (m3/s)`, df$conc)
  plot_df <- data.frame(
    nq = H$nq,
    nc = H$nc,
    t = 1:length(H$nq))
  #p = plots
  p <- ggplot(plot_df, aes(nq, nc, color = t)) +
    geom_path(linewidth = 2.5) +
    scale_color_viridis_c() +
    labs(
      title = paste0("Hysteresis – ", compound, " (Storm 15)"),
      x = "Normalized Discharge",
      y = paste0("Normalized ", compound),
      color = "Time →") +
    theme_classic(base_size = 18) #could do theme minimal?
  #this ggsave save plots to your working directory
  #run this if you don't know where that is btw ->
  #getwd()
  ggsave(file.path(out_dir, paste0("Hysteresis_", gsub(" ", "_", compound), ".png")),
         p, width = 8, height = 7, dpi = 300)
  #this spits back the Flushing and Hysteresis index for the stuff on the compound list.
  message("\n", compound,
          "\n  FI = ", round(H$FI, 3),
          "\n  HI = ", round(mean(H$HI), 3), "\n")
}

