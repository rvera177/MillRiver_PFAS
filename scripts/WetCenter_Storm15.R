

#Mill River PFAS Project
#Created November 28, 2025
#Updated November 30, 2025
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

AllChem_WetCenter_SE15 <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Storm_Event_15_all_results.csv")
WetCenterDischarge <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenterDischarge.csv")

#first, i want to add the discharge data into AllChem_WetCenter_SE15
#this next step is organizing and cleaning up my discharge data
WetCenterDischarge$Time <- mdy_hm(WetCenterDischarge$Time)

#using rating curve equation to calculate flow based on gauge height (feet)
WetCenterDischarge <- WetCenterDischarge %>%
  mutate(Adjusted_Gauge_meters = `Adjusted Gauge Height` *0.3048,  # 1 ft = 0.3048 m
    `Flow (m3/s)` = 2.74*(Adjusted_Gauge_meters^1.73)) # rating curve in m3/s

#now i want to get the newly created Flow column and put it in AllChem_WetCenter_SE15

flow <- approx( #lining up air data with stream data chronologically
  x = WetCenterDischarge$Time,
  y = WetCenterDischarge$`Flow (m3/s)`,
  xout = AllChem_WetCenter_SE15$DateTime)$y

AllChem_WetCenter_SE15$Flow <- flow #adds flow to the stream chemistry dataset in a new "Flow" column


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
  hc_title(text = "Storm 15 results") %>%
  hc_yAxis_multiples(
    list(title = list(text = "PFHxS (ng/L)"), lineWidth = 3, lineColor = "orange", min = min(AllChem_WetCenter_SE15$`PFHxS Results`), max = max(AllChem_WetCenter_SE15$`PFHxS Results`)),
    list(title = list(text = "PFOA (ng/L)"), lineWidth = 3, lineColor = "darkorange", min = min(AllChem_WetCenter_SE15$`PFOA Results`), max = max(AllChem_WetCenter_SE15$`PFOA Results`)),
    list(title = list(text = "PFOS (ng/L)"), lineWidth = 3, lineColor = "red",   min = min(AllChem_WetCenter_SE15$`PFOS Results`), max = max(AllChem_WetCenter_SE15$`PFOS Results`)),
    list(title = list(text = "PFHpA (ng/L)"), lineWidth = 3, lineColor = "rosybrown",   min = min(AllChem_WetCenter_SE15$`PFHpA Results`), max = max(AllChem_WetCenter_SE15$`PFHpA Results`)),
    list(title = list(text = "PFDA (ng/L)"), lineWidth = 3, lineColor = "darkred",   min = min(AllChem_WetCenter_SE15$`PFDA Results`), max = max(AllChem_WetCenter_SE15$`PFDA Results`)),
    list(title = list(text = "Flow (m3/s)"), lineWidth = 3, lineColor = "darkblue",   min = min(AllChem_WetCenter_SE15$Flow), max = max(AllChem_WetCenter_SE15$Flow))
  ) %>%
  #hc_add_theme(hc_theme_sandsignika()) %>%
  #adding in the flow time series first so it's in the background
  hc_add_series(name = "Flow",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Flow)),
                yAxis = 5, type="area", lineWidth=5, color = 'darkblue') %>%
   hc_plotOptions(type = "line", series = list(connectNulls = TRUE)) %>% #this makes the lines plot across NA values
  hc_add_series(name = "PFHxS",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `PFHxS Results`)),
                 yAxis = 0, lineWidth = 5, color = "orange") %>%
  hc_add_series(name = "PFOA",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `PFOA Results`)),
                yAxis = 1, lineWidth = 5, color = 'darkorange') %>%
  hc_add_series(name = "PFOS",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `PFOS Results`)),
                yAxis = 2, lineWidth = 5,color = 'red') %>%
  hc_add_series(name = "PFHpA",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `PFHpA Results`)),
                yAxis = 3, lineWidth = 5,color = "rosybrown") %>%
  hc_add_series(name = "PFDA",
              data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = `PFDA Results`)),
               yAxis = 4, lineWidth = 5, color = 'darkred')

cor_data <- read.csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Zyna_%20storm%2015%20sampling%20All%20Chem%20results.csv")
#correlation
cor_results <- cor_data %>%
  select(-c(data.point, Date, Time, Ammonium, Bromide, Bromate, Phosphate, Chlorite, Chlorate, Flouride)) %>%
  cor(method = "spearman", use = "complete.obs")
#library()  # I forget which package is for melt
cor_long <- melt(cor_results)

#plotting the correlation results. 
library(pheatmap)
pheatmap(cor_results,
         color = colorRampPalette(c("red", "white", "blue"))(50),
         display_numbers = TRUE,
         breaks = seq(-1, 1, length.out = 51),
         main = "Spearman Correlation Heatmap All Chem")

ggplot(cor_long, aes(x = Variable1, y = Variable2, fill = Spearman_r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Spearman_r, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(hjust = 1)
  )




highchart() %>%
  hc_title(text = "Storm 15 IC/TOC results") %>%
  hc_yAxis_multiples(
    list(title = list(text = "DOC"), lineWidth = 3, lineColor = "lightgreen", min = min(AllChem_WetCenter_SE15$TOC), max = max(AllChem_WetCenter_SE15$TOC)),
    list(title = list(text = "TOC"), lineWidth = 3, lineColor = 'forestgreen', min = min(AllChem_WetCenter_SE15$DOC), max = max(AllChem_WetCenter_SE15$DOC)),
    list(title = list(text = "TN)"), lineWidth = 3, lineColor = "black",   min = min(AllChem_WetCenter_SE15$TN), max = max(AllChem_WetCenter_SE15$TN)),
    list(title = list(text = "Nitrate"), lineWidth = 3, lineColor = "maroon",   min = min(AllChem_WetCenter_SE15$Nitrate), max = max(AllChem_WetCenter_SE15$Nitrate)),
    list(title = list(text = "6:2 FTS (ng/L)"), lineWidth = 3, lineColor = "red",   min = min(AllChem_WetCenter_SE15$`6:2FTS Results`), max = max(AllChem_WetCenter_SE15$`6:2FTS Results`)),
    list(title = list(text = "Flow (m3/s)"), lineWidth = 3, lineColor = "blue",   min = min(AllChem_WetCenter_SE15$Flow), max = max(AllChem_WetCenter_SE15$Flow))
  ) %>%
  #hc_add_theme(hc_theme_sandsignika()) %>%
  #adding in the flow time series first so it's in the background
  hc_add_series(name = "Flow",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Flow)),
                yAxis = 5, type="area", lineWidth=5, color = 'blue') %>%
  hc_plotOptions(type = "line", series = list(connectNulls = TRUE)) %>% #this makes the lines plot across NA values
  hc_add_series(name = "TOC",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = TOC)),
                yAxis = 0, lineWidth = 5, color = "lightgreen") %>%
  hc_add_series(name = "DOC",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = DOC)),
                yAxis = 1, lineWidth = 5, color = 'forestgreen') %>%
  hc_add_series(name = "TN",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = TN)),
                yAxis = 2, lineWidth = 5,color = 'black') %>%
  hc_add_series(name = "Nitrate",
                data = list_parse2(AllChem_WetCenter_SE15 %>% transmute(x = DateTime, y = Nitrate)),
                yAxis = 3, lineWidth = 5,color = "maroon") %>%
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

max_flow <- 3.5  # Maximum discharge on the y axis for positioning the rainfall bars
max_precip <- 1  # Maximum rainfall on the y axis for scaling against flow
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
  labs(title = "Storm 16 Hydrograph (2 inch rainfall)", x = "Time") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 12))



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

