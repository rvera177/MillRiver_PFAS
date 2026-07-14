library(tidyverse)
library(lubridate)
library(readr)
library(ggplot2)
library(dplyr)

# STORM DETECTION + HYDROGRAPH SEPARATION
# Code by Raul Vera
# I created this in April 2026, and made slight updates on July 2026
# still very rough and preliminary. 
# the code does the following steps
#1) storm find function and finding storms based off of flow at the wet center
#2) event/preevent seperation function applied on the found storms using the isotope inventory results 
# (needs at least 8 stream isotope and 1 precip isotope to actualy work for a given storm)
#3) that's it. lol.
# there are a few isotope samples that have not been run, so there are a few storms that can be added to this. 


WetCenter_Flow_Combined <- read.csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/WetCenterFlow_Updated_April_14_2026.csv") %>%
  mutate(Time = parse_date_time(Time, orders = c("mdY HM", "mdY")))

sum(is.na(WetCenter_Flow_Combined$Time))
range(WetCenter_Flow_Combined$Time)

WetCenter_samples <- read.csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/RAW_Isotope_Inventory_Downloaded_July_14_26.csv")


library(zoo)  # for rolling mean
storm_find <- function(flow_data,
                       time_col          = "Time",
                       flow_col          = "Flow_m3_s",
                       baseline_days     = 3,
                       rise_factor       = 1.5,
                       min_duration_hrs  = 6,
                       recession_factor  = 1.4,
                       pre_event_days    = 3,
                       min_gap_hrs       = 12) {
  
  df <- flow_data %>%
    rename(Time = all_of(time_col), Flow = all_of(flow_col)) %>%
    mutate(Time = as.POSIXct(Time, format = "%m/%d/%Y %H:%M", tz = "America/New_York")) %>%
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
          filter(Time >= .x - ddays(baseline_days),
                 Time <  .x - ddays(baseline_days - 1)) %>%
          summarise(m = median(Flow, na.rm = TRUE)) %>%
          pull(m)
      }),
      
      # Start: either pre_event_days before peak OR trough before this storm
      # whichever is LATER — prevents overlap with previous storm
      start = map2_dbl(peak_time, trough_before, ~ {
        candidate <- .x - ddays(pre_event_days)
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
  baseline_days    = 2.5,
  rise_factor      = 1.14,
  min_duration_hrs = 8,
  pre_event_days   = 3,
  recession_factor = 1.2,
  min_gap_hrs      = 24
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
            aes(x = Time, y = Flow_m3_s, color = Source),
            linewidth = 0.6) +
  geom_point(data = storms_detected,
             aes(x = peak_time, y = peak),
             color = "blue", size = 3) +
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

isotope_data = read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/RAW_Isotope_Inventory_Downloaded_July_14_26.csv")
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
      Q_pre_event  = f_pre_event * Flow_m3_s,
      Q_event      = (1 - f_pre_event) * Flow_m3_s
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
                aes(x = Time, ymin = Q_pre_event, ymax = Flow_m3_s),
                fill = "blue", alpha = 0.6) +
    # --- Total flow line on top for reference ---
    geom_line(data = flow_window,
              aes(x = Time, y = Flow_m3_s),
              color = "steelblue", linewidth = 0.8) +
    # --- Isotope points at sample times ---
    geom_line(data = d,
              aes(x = DateTime,
                  y = (.data[[if(isotope=="d18O") "d18O" else "d2H"]] -
                         min(.data[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) *
                    (max(flow_window$Flow_m3_s, na.rm=TRUE) /
                       diff(range(.data[[if(isotope=="d18O") "d18O" else "d2H"]],
                                  na.rm=TRUE))),
                  color = "Stream Isotope"),
              linewidth = 2, linetype = "dotdash") +
    scale_color_manual(values = c("Stream Isotope" = "magenta")) +
    scale_y_continuous(
      name     = "Discharge (m³/s)",
      sec.axis = sec_axis(
        ~ . * diff(range(d[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) /
          max(flow_window$Flow_m3_s, na.rm=TRUE) +
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


#plotting total flow only

plot_separation <- function(sep_data, storm_id_val, isotope = "d18O") {
  d <- sep_data %>% filter(storm_id == storm_id_val)
  
  if (nrow(d) == 0) {
    message("No data for storm ", storm_id_val)
    return(NULL)
  }
  
  # --- Get field name ---
  field_name <- storm_labels %>%
    filter(storm_id == storm_id_val) %>%
    pull(field_name)
  
  storm_info  <- storms_detected %>% filter(storm_id == storm_id_val)
  flow_window <- WetCenter_Flow_Combined %>%
    filter(Time >= storm_info$start, Time <= storm_info$end)
  
  # --- Interpolate fractions to 10-min resolution ---
  separation_interp <- flow_window %>%
    mutate(
      f_pre_event = approx(
        x    = d$DateTime,
        y    = d$f_pre_event,
        xout = Time,
        rule = 2)$y
    ) %>%
    mutate(
      f_pre_event = pmax(0, pmin(1, f_pre_event)),
      Q_pre_event = f_pre_event * Flow_m3_s,
      Q_event     = (1 - f_pre_event) * Flow_m3_s
    )
  
  # --- Plot 1: Total flow only ---
  p_flow <- ggplot() +
    geom_area(data = flow_window,
              aes(x = Time, y = Flow_m3_s),
              fill = "steelblue", alpha = 0.5, color = "black",
              linewidth = 0.8) +
    geom_point(data = d,
               aes(x = DateTime, y = Q_total),
               color = "red", size = 2) +
    scale_y_continuous(name = "Discharge (m³/s)") +
    scale_x_datetime(date_labels = "%b %d", date_breaks = "1 day") +
    labs(title    = paste(field_name, "— Total Flow"),
         x = NULL) +
    theme_bw()
  
  # --- Plot 2: Hydrograph separation ---
  p_sep <- ggplot() +
    geom_ribbon(data = separation_interp,
                aes(x = Time, ymin = 0, ymax = Q_pre_event),
                fill = "steelblue", alpha = 0.6) +
    geom_ribbon(data = separation_interp,
                aes(x = Time, ymin = Q_pre_event, ymax = Flow_m3_s),
                fill = "blue", alpha = 0.6) +
    geom_line(data = flow_window,
              aes(x = Time, y = Flow_m3_s),
              color = "black", linewidth = 0.8) +
    geom_line(data = d,
              aes(x = DateTime,
                  y = (.data[[if(isotope=="d18O") "d18O" else "d2H"]] -
                         min(.data[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) *
                    (max(flow_window$Flow_m3_s, na.rm=TRUE) /
                       diff(range(.data[[if(isotope=="d18O") "d18O" else "d2H"]],
                                  na.rm=TRUE))),
                  color = "Stream Isotope"),
              linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("Stream Isotope" = "darkgreen")) +
    scale_y_continuous(
      name     = "Discharge (m³/s)",
      sec.axis = sec_axis(
        ~ . * diff(range(d[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE)) /
          max(flow_window$Flow_m3_s, na.rm=TRUE) +
          min(d[[if(isotope=="d18O") "d18O" else "d2H"]], na.rm=TRUE),
        name = paste0(isotope, " (‰)")
      )
    ) +
    scale_x_datetime(date_labels = "%b %d", date_breaks = "1 day") +
    labs(title    = paste(field_name, "— Hydrograph Separation (", isotope, ")"),
         subtitle = paste("C_pre =", round(d$C_pre[1], 2),
                          "‰  |  C_event =", round(d$C_event[1], 2),
                          "‰  |  Mean event fraction:",
                          round(mean(d$f_event, na.rm=TRUE), 2)),
         x = "Date/Time", color = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # --- Stack vertically with patchwork ---
  p_flow / p_sep
}

# --- Plot all storms ---
library(patchwork)

plots <- map(unique(separated_d18O$storm_id),
             ~ plot_separation(separated_d18O, .x, "d18O"))

# View each storm individually
for (p in plots) print(p)

# Or all together in a grid
wrap_plots(plots, ncol = 3) +
  plot_annotation(title = "Hydrograph Separation — All Storms (δ¹⁸O)")



library(dygraphs)
library(xts)

# --- Convert to xts format which dygraphs requires ---
flow_xts <- WetCenter_Flow_Combined %>%
  arrange(Time) %>%
  select(Time, Flow_m3_s) %>%
  as.data.frame() %>%
  { xts(.$Flow_m3_s, order.by = .$Time) }

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
  rowwise() %>%
  mutate(
    SampleFlow = if (any(abs(as.numeric(difftime(WetCenter_isotope$DateTime, Time, units = "mins"))) <= 7.5, na.rm = TRUE)) {
      Flow_m3_s
    } else {
      NA_real_
    }
  ) %>%
  ungroup() %>%
  select(Time, Flow_m3_s, SampleFlow) %>%
  { xts(.[, c("Flow_m3_s", "SampleFlow")], order.by = .$Time) }

# --- Define alternating colors ---
shading_colors <- c("#ADD8E650", "#FFA50050")  # light blue and light orange

p <- dygraph(sample_xts, main = "Mill River — Combined Flow Record") %>%
  dyAxis("y", label = "Discharge (m³/s)") %>%
  dyOptions(fillAlpha = 0.2) %>%
  dySeries("Flow_m3_s", label = "Discharge", color = "blue", fillGraph = TRUE) %>%
  dySeries("SampleFlow", label = "Isotope Sample", color = "purple",
           drawPoints = TRUE, pointSize = 5, strokeWidth = 0) %>%
  dyRangeSelector() %>%
  dyHighlight(highlightSeriesOpts = list(strokeWidth = 2))

# --- Add alternating storm shading + edge lines ---
for (i in seq_len(nrow(storms_detected))) {
  
  # Alternating fill color
  fill_color <- shading_colors[(i %% 2) + 1]
  
  # Shaded window
  p <- p %>%
    dyShading(from  = format(storms_detected$start[i], "%Y-%m-%dT%H:%M:%S"),
              to    = format(storms_detected$end[i],   "%Y-%m-%dT%H:%M:%S"),
              color = fill_color)
  
  # Dark vertical line at storm start
  p <- p %>%
    dyAnnotation(format(storms_detected$start[i], "%Y-%m-%dT%H:%M:%S"),
                 text       = "|",
                 tooltip    = paste("Storm", storms_detected$storm_id[i], "start"),
                 attachAtBottom = FALSE)
  
  # Dark vertical line at storm end  
  p <- p %>%
    dyAnnotation(format(storms_detected$end[i], "%Y-%m-%dT%H:%M:%S"),
                 text       = "|",
                 tooltip    = paste("Storm", storms_detected$storm_id[i], "end"),
                 attachAtBottom = FALSE)
}

p

