

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

#i wonder if there's a way to make my github a working directory.
#or a way to shorten urls. so i don't have to write long url's

PFAS_WetCenter_2025 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenterPFASResults.csv")
WetCenterDischarge <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenterDischarge.csv")
AllChem_WetCenter_SE15 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Storm_Event_15_all_results.csv")

#now, i want to add the discharge data into AllChem_WetCenter_SE15
#this next step is organizing and cleaning up my discharge data
WetCenterDischarge$Time <- mdy_hm(WetCenterDischarge$Time)
AllChem_WetCenter_SE15$DateTime <- mdy_hm(AllChem_WetCenter_SE15$DateTime)
#using rating curve equation to calculate flow based on gauge height (feet)
WetCenterDischarge <- WetCenterDischarge %>%
  mutate(Adjusted_Gauge_meters = `Adjusted Gauge Height` *0.3048,  # 1 ft = 0.3048 m
         `Flow (m3/s)` = 2.74*(Adjusted_Gauge_meters^1.73)) # rating curve in m3/s

#now i want to get the newly created Flow column and put it in AllChem_WetCenter_SE15

flow <- approx(
  x = WetCenterDischarge$Time,
  y = WetCenterDischarge$`Flow (m3/s)`,
  xout = AllChem_WetCenter_SE15$DateTime,
  rule = 2   # enables extrapolation
)$y

AllChem_WetCenter_SE15$Flow <- flow #adds flow to the stream chemistry dataset in a new "Flow" column

#this plot is all ofthe flow data i have for wet center. this doesn't use
plot(WetCenterDischarge$Time, WetCenterDischarge$`Flow (m3/s)`)
#there's some outliers.For some reason, flow will randomly be 1.49. Machine error definitely. removing them here
# i'm using a 0.0000001 approximation because it's the only thing that worked
WetCenterDischarge <- WetCenterDischarge %>%
  filter(abs(`Flow (m3/s)` - 1.4927567) > 0.000001)

#nice, fixed.
plot(WetCenterDischarge$Time, WetCenterDischarge$`Flow (m3/s)`)


plot(AllChem_WetCenter_SE15$DateTime, flow)

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
#so this isn't exactly pfas 6 for the EPA drinking water regulations.
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
  scale_x_datetime(expand = c(0,0))+
  
  #this is for greying out the time i did spatial collection
  #just for context. might delete for the presentation
  geom_rect(
    aes(xmin = as.POSIXct("2025-10-11 08:00:00"), xmax = as.POSIXct("2025-10-11 15:00:00"),
        ymin = -Inf, ymax = Inf), fill = "grey30", alpha = 0.4) +
  # geom_rect(
  #  aes(xmin = as.POSIXct("2025-10-15 07:20:00"),xmax = as.POSIXct("2025-10-15 15:00:00"),
  #   ymin = -Inf, ymax = Inf), fill = "grey30", alpha = 0.4) +
  #Adding labels for the spatial sets. just for show
  annotate("text", x = as.POSIXct("2025-10-11 12:00:00"),  #midpoint of S1
           y = max(Storm15$`Flow (m3/s)`, na.rm = TRUE) * 0.75,  
           label = "S1", size = 5, fontface = "bold") +
  # annotate("text", x = as.POSIXct("2025-10-15 11:10:00"),  # midpoint of S2
  #    y = max(Storm15$`Flow (m3/s)`, na.rm = TRUE) * 0.75,
  #    label = "S2", size = 5,fontface = "bold") +
  theme_classic()+
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12),   # x-axis numbers
        axis.text.y = element_text(size = 12),  # y-axis numbers
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 12)) # stuff in legend



#--------example from stack overflow----------------
library(highcharter)
set.seed(1)
n <- 100
x1 <- cumsum(rnorm(n))
x2 <- cumsum(runif(n)-0.5)+10
x3 <- cumsum(rnorm(n,0,20))+100
x4 <- cumsum(rnorm(n,0,20))+1000

highchart() %>% 
  hc_add_series(data = x1) %>% 
  hc_add_series(data = x2, yAxis = 1) %>% 
  hc_add_series(data = x3, yAxis = 2) %>%
  hc_add_series(data = x4, yAxis = 3) %>%
  hc_yAxis_multiples(
    list(lineWidth = 3, lineColor='#7cb5ec', title=list(text="First y-axis")),
    list(lineWidth = 3, lineColor="#434348", title=list(text="Second y-axis")),
    list(lineWidth = 3, lineColor="#90ed7d", title=list(text="Third y-axis")),
    list(lineWidth = 3, lineColor="#f7a35c", title=list(text="Fourth y-axis"))
  )



#--------my attempt at highcharts------------

highchart() %>%
  hc_title(text = "Storm 15 PFAS results") %>%
  hc_yAxis_multiples(
    list(title = list(text = "PFOA (ng/L)"), lineWidth = 3, lineColor = "red", min = min(AllChem_WetCenter_SE15$`PFOA Results`), max = max(AllChem_WetCenter_SE15$`PFOA Results`)),
    list(title = list(text = "Sodium"), lineWidth = 3, lineColor = "green",   min = min(AllChem_WetCenter_SE15$Sodium), max = max(AllChem_WetCenter_SE15$Sodium)),
    list(title = list(text = "PFHpA (ng/L)"), lineWidth = 3, lineColor = "darkred",   min = min(AllChem_WetCenter_SE15$`PFHpA Results`), max = max(AllChem_WetCenter_SE15$`PFHpA Results`)),
    list(title = list(text = "Flow (m3/s)"), lineWidth = 3, lineColor = "blue",   min = min(AllChem_WetCenter_SE15$Flow), max = max(AllChem_WetCenter_SE15$Flow)),
    list(title = list(text = "Chloride"), lineWidth = 3, lineColor = "grey",   min = min(AllChem_WetCenter_SE15$Chloride), max = max(AllChem_WetCenter_SE15$Chloride))
  ) %>%
  #hc_add_theme(hc_theme_sandsignika()) %>%
  #adding in the flow time series first so it's in the background
  hc_add_series(name = "Flow",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Flow)),
                yAxis = 3, type="area", lineWidth=5, color = 'blue') %>%
   hc_plotOptions(type = "line", series = list(connectNulls = TRUE)) %>% #this makes the lines plot across NA values

  hc_add_series(name = "PFOA",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `PFOA Results`)),
                yAxis = 0, lineWidth = 5, color = 'red') %>%
  hc_add_series(name = "Sodium",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Sodium)),
                yAxis = 1, lineWidth = 5,color = "green") %>%
  hc_add_series(name = "PFHpA",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `PFHpA Results`)),
                yAxis = 2, lineWidth = 5,color = "darkred") %>%
hc_add_series(name = "Chloride",
              data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Chloride)),
              yAxis = 4, lineWidth = 5, color = "grey")

cor_data <- read.csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Zyna_%20storm%2015%20sampling%20All%20Chem%20results.csv")
#correlation. Removing columns that are time related or only had zero's
cor_results <- cor_data %>%
  select(-c(data.point, Date, Time, Ammonium, Bromide, Bromate, Phosphate, Chlorite, Chlorate, Flouride)) %>%
  cor(method = "spearman", use = "complete.obs")
library(reshape2)  # this is for melting
cor_long <- melt(cor_results)

#correlation results. Red is supposed to be negative correlation, so that's fixed now
library(pheatmap)
pheatmap(cor_results,
         color = colorRampPalette(c("red", "white", "blue"))(50),
         display_numbers = TRUE,
         breaks = seq(-1, 1, length.out = 51),
         main = "Spearman Correlation Heatmap All Chem")

cor_results[lower.tri(cor_results, diag = TRUE)] <- NA
cor_long <- melt(cor_results, varnames = c("Variable1", "Variable2"), value.name = "Spearman_r", na.rm = TRUE)

# Reordering Variable2 so labels are aligned with color with tiles
cor_long$Variable2 <- factor(cor_long$Variable2, levels = rev(unique(cor_long$Variable2)))

ggplot(cor_long, aes(x = Variable1, y = Variable2, fill = Spearman_r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Spearman_r, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1))




highchart() %>%
  hc_title(text = "Storm 15 IC/TOC results") %>%
  hc_yAxis_multiples(
    list(title = list(text = "DN"), lineWidth = 3, lineColor = "darkgreen", min = min(AllChem_WetCenter_SE15$DN), max = max(AllChem_WetCenter_SE15$DN)),
    list(title = list(text = "DOC"), lineWidth = 3, lineColor = 'orange', min = min(AllChem_WetCenter_SE15$DOC), max = max(AllChem_WetCenter_SE15$DOC)),
    list(title = list(text = "TN)"), lineWidth = 3, lineColor = "forestgreen",   min = min(AllChem_WetCenter_SE15$TN), max = max(AllChem_WetCenter_SE15$TN)),
    list(title = list(text = "Nitrate"), lineWidth = 3, lineColor = "green",   min = min(AllChem_WetCenter_SE15$Nitrate), max = max(AllChem_WetCenter_SE15$Nitrate)),
    list(title = list(text = "6:2 FTS (ng/L)"), lineWidth = 3, lineColor = "red",   min = min(AllChem_WetCenter_SE15$`6:2FTS Results`), max = max(AllChem_WetCenter_SE15$`6:2FTS Results`)),
    list(title = list(text = "Flow (m3/s)"), lineWidth = 3, lineColor = "blue",   min = min(AllChem_WetCenter_SE15$Flow), max = max(AllChem_WetCenter_SE15$Flow))
  ) %>%
  #hc_add_theme(hc_theme_sandsignika()) %>%
  #adding in the flow time series first so it's in the background
  hc_add_series(name = "Flow",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Flow)),
                yAxis = 5, type="area", lineWidth=5, color = 'blue') %>%
  hc_plotOptions(type = "line", series = list(connectNulls = TRUE)) %>% #this makes the lines plot across NA values
  hc_add_series(name = "DN",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = DN)),
                yAxis = 0, lineWidth = 5, color = "darkgreen") %>%
  hc_add_series(name = "DOC",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = DOC)),
                yAxis = 1, lineWidth = 5, color = 'orange') %>%
  hc_add_series(name = "TN",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = TN)),
                yAxis = 2, lineWidth = 5,color = 'forestgreen') %>%
  hc_add_series(name = "Nitrate",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Nitrate)),
                yAxis = 3, lineWidth = 5,color = "green") %>%
  hc_add_series(name = "6:2 FTS",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `6:2FTS Results`)),
                yAxis = 4, lineWidth = 5, color = 'red')





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

