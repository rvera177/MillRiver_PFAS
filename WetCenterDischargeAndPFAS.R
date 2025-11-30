

#Mill River PFAS Project
#Created November 28, 2025
#Updated November 30, 2025
#Raul Vera


#Graphing results for storm 15 and 16
#storm 15: October 12-13. 0.7 inches of rainfall.
          #samples from every other isco bottle for 12 total samples

#Storm 16: October 30-31. 2 inches of rainfall
          #samples from every isco bottle for 24 total samples

library(readr)
library(lubridate)
library(ggplot2)
library(dplyr) #for pipes and renaming columns
library(tidyr) # for pivot_longer

#i wonder if there's a way to make my github a working directory.
#so i don't have to write long url's

PFAS_WetCenter_2025 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/WetCenterPFASResults.csv")
WetCenterDischarge <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/WetCenterDischarge.csv")

WetCenterDischarge$Time <- mdy_hm(WetCenterDischarge$Time)
#using rating curve equation to calculate flow based on gauge height (feet)
WetCenterDischarge <- WetCenterDischarge %>%
  mutate(Adjusted_Gauge_meters = `Adjusted Gauge Height` *0.3048,  # 1 ft = 0.3048 m
    `Flow (m3/s)` = 2.74*(Adjusted_Gauge_meters^1.73)) # rating curve in m3/s

PFAS_WetCenter_2025 <- PFAS_WetCenter_2025 %>%
  mutate(DateTime = parse_date_time(DateTime,
      orders = c("mdy HM p", "mdy HMS p"), tz = "America/New_York"))#changing the hour format so it works

Start.Time <-  ymd_hm("2025-10-10 00:00")
End.Time <-  ymd_hm("2025-10-16 00:00")
PFAS_storm15 <- PFAS_WetCenter_2025 %>%
  filter(Type == "Stream") %>%                 
  filter(DateTime >= Start.Time, DateTime <= End.Time)


#RV: remember to remove gauge height values over = 2.30952 ft from WetCenterDischarge at some point
#these are random outliers from the machine. 

Storm15 <- WetCenterDischarge %>%
  filter(Time >= Start.Time,
         Time <= End.Time)
scaling=0.02 #i'm using this to scale the pfAS against discharge.

#plotting multiple at one time
PFAS6 = c("PFOA Results","PFOS Results","PFHxS Results", "PFNA Results","HFPO-DA Results","PFBS Results")
PFAS_storm15 <- PFAS_storm15 %>%
  mutate(Sum_PFAS6 = rowSums(across(all_of(PFAS6)), na.rm = TRUE))

pfas_vars <- c("Sum_PFAS6", "PFOS Results", "PFOA Results", "PFHxA Results", "PFBA Results", "PFDA Results", "PFHpA Results")
PFAS_long <- PFAS_storm15 %>%
  pivot_longer(
    cols = all_of(pfas_vars),
    names_to = "Compound",
    values_to = "Concentration")
#so PFAS_long makes it so each sample time has n(variables) rows. 
#something with how pivot longer works, but it doesn't seem to be an issue

#quick plot
ggplot(PFAS_long, aes(x = DateTime, y = Concentration, color = Compound)) +
  geom_point(size = 3) + geom_line() + theme_classic()

#plot against flow
ggplot() +

geom_line( data = Storm15,
  aes(x = Time, y = `Flow (m3/s)`), color = "blue", linewidth = 1.2) +
geom_line( data = PFAS_long,
  aes(x = DateTime, y = Concentration * scaling, color = Compound), linewidth = 1) +
geom_point(data = PFAS_long,
    aes(x = DateTime, y = Concentration * scaling, color = Compound), size = 2) +
scale_y_continuous( name = "Discharge (m³/s)",
  sec.axis = sec_axis(~ . / scaling, name = "PFAS (ng/L)")) +
  labs( title = "Discharge and PFAS Analytes for Storm 15",
    x = "Time", color = "Analyte") +

  #this is for greying out the time i did spatial collection
  #just for context. might delete for the presentation
  geom_rect(
    aes(xmin = as.POSIXct("2025-10-11 08:00:00"), xmax = as.POSIXct("2025-10-11 15:00:00"),
      ymin = -Inf, ymax = Inf), fill = "grey30", alpha = 0.4) +
  geom_rect(
    aes(xmin = as.POSIXct("2025-10-15 07:20:00"),xmax = as.POSIXct("2025-10-15 15:00:00"),
      ymin = -Inf, ymax = Inf), fill = "grey30", alpha = 0.4) +
#Adding labels for the spatial sets. just for show
annotate("text", x = as.POSIXct("2025-10-11 12:00:00"),  #midpoint of S1
    y = max(Storm15$`Flow (m3/s)`, na.rm = TRUE) * 0.75,  
    label = "S1", size = 5, fontface = "bold") +
  annotate("text", x = as.POSIXct("2025-10-15 11:10:00"),  # midpoint of S2
    y = max(Storm15$`Flow (m3/s)`, na.rm = TRUE) * 0.75,
    label = "S2", size = 5,fontface = "bold") +
  theme_classic()+
  theme(plot.title = element_text(size = 15, face = "bold"),
  axis.title.x = element_text(size = 14), 
  axis.title.y = element_text(size = 14), 
  axis.text.x = element_text(size = 12),   # x-axis numbers
  axis.text.y = element_text(size = 12),  # y-axis numbers
  legend.title = element_text(size = 13), 
  legend.text = element_text(size = 12)) # stuff in legend



#-----Graphing storm 16------------------



#now for storm 16, which was the Halloween storm. 
# I don't have stage height. I have pressure data. 
#this needs to be converted to stage height,
#and then I'll use rating curve equation for flow


#---bring in the pressure data---
WetStation3Air <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/WetStation3Air.csv")
WetStation4Stream <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/WetStation4Stream.csv")

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
#PrecipitationData_Storms15 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/PrecipitationData_Storms15.csv")
PrecipitationData_Storms16 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/PrecipitationData_Storms16.csv")
PFAS_Inventory_temp <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/PFAS%20Inventory%20-%20temp.csv")

#---mutating the DateTime format so it's consistent.
#PrecipitationData_Storms15 <- PrecipitationData_Storms15 %>%
#  mutate(DateTime = mdy_hm(DateTime, tz = "America/New_York"))

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

max_flow <- 3  # Maximum discharge on the y axis for positioning the rainfall bars
max_precip <- 1.5  # Maximum rainfall on the y axis for scaling against flow
scale_factor <- max_flow / max_precip  #scaling again

#plot for storm 16.
ggplot() +
  geom_line(data = WetStation4Stream,
    aes(x = Time, y = Flow), color = "blue",linewidth = 1.2) +
  geom_rect(
    data = PrecipitationData_Storms16, aes(
      xmin = DateTime - seconds(7.5*60), # widening rain bars for visibility
      xmax = DateTime + seconds(7.5*60),  
      ymin = max_flow - `Precip (in)` * scale_factor,
      ymax = max_flow),fill = "blue3",alpha = 0.5) +
  geom_point(data = PFAS_points,
    aes(x = DateTime, y = Flow,), color = "red",size = 3) +
  geom_label(
    data = PFAS_points, aes(x = DateTime, y = Flow, label = Sample),
    vjust = -0.7, # positioning above the point
    fill = "white", # background label color
    color = "red",  # label text colour
    size = 3, label.size = 0.2, # border thickness around label
    fontface = "bold", alpha = 0.9)+
  scale_y_continuous(   # Y axis: left = Flow, right = Rainfall
    name = "Discharge (m³/s)",sec.axis = sec_axis(~ (max_flow - .) #left side
      / scale_factor, name = "Rainfall (in/hr)"), #right side
    expand = c(0,0)) + #expand(0,0) removes the whitespace buffer around the plots
  scale_x_datetime(
    limits = as.POSIXct(c("2025-10-29 20:00:00", "2025-11-07 18:00:00")),
    date_labels = "%b %d",
    date_breaks = "1 day", expand = c(0,0)) +
  labs(title = "Storm 16 Hydrograph (Halloween Storm)", x = "Time") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12))

