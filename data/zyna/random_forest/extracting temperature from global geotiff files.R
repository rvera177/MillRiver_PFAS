library(terra)
library(dplyr)
library(readxl)

# Set working directory
setwd("C:/Users/vera780/OneDrive - PNNL/Documents/ArcGIS/Projects/Mapping K and S/wc2.1_30s_tavg")

# Load rasters
f <- list.files(path = getwd(), pattern = "\\.tif$", full.names = TRUE)
rasters <- rast(f)

# Fill small missing gaps in rasters using focal mean (5x5 window ~25km²)
rasters_filled <- focal(rasters, w = matrix(1, 9, 9), fun = mean, na.policy = "only")


# Read coordinate data 
coords <- read_excel("RV Tea_bag_studies_Metadata.xlsx")
colnames(coords)[8:9] <- c("longitude", "latitude")  # Correct order

# Convert to spatial points and match CRS
points <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
points <- project(points, crs(rasters_filled))  # Ensure CRS match

# Extract precipitation values using IDW (Inverse Distance Weighting) interpolation
#precip_values <- extract(rasters_filled, points, method = "idw", radius = 30.0, power = 2)
temp_values <- extract(rasters_filled, points, method = "nearest", radius = 70.0)


# Replace any remaining extreme NoData values
temp_values[temp_values == -32768] <- NA  

# Merge extracted values with original dataset
temp_data <- cbind(coords, temp_values[, -1])  # Remove ID column

# Rename precipitation columns dynamically
num_temp_months <- ncol(temp_values) - 1  # Number of extracted months
colnames(temp_data)[(ncol(coords) + 1):ncol(temp_data)] <- paste0("Month_", 1:num_temp_months)

# Plot raster and points
# Plot only the first raster (January)
plot(rasters_filled[[1]], main = "January Temperature")

# Save to CSV
write.csv(temp_data, "extracted_temperature.csv", row.names = FALSE)


#-----Rainfall during incubation period------- 
library(readxl)
library(dplyr)
library(lubridate)

# Load the dataset
data <- read_excel("RV Temp Tea_bag_studies_Metadata.xlsx")  # Update with your actual file name

# Convert dates to proper format
data <- data %>%
  mutate(Start_Date = as.Date('Start of Experiment (mm/dd/yy)', format = "%m/%d/%y"),
         End_Date = as.Date('End of Experiment (mm/dd/yy)', format = "%m/%d/%y"))

# Function to calculate weighted precipitation
calculate_weighted_precip <- function(start_date, end_date, precip_values) {
  # Extract start and end months
  start_month <- month(start_date)
  end_month <- month(end_date)
  
  # Days in each month
  start_days <- as.numeric(as.Date(paste(year(start_date), start_month, days_in_month(start_date), sep = "-")) - start_date + 1)
  end_days <- as.numeric(end_date - as.Date(paste(year(end_date), end_month, "01", sep = "-")) + 1)
  
  # Total experiment duration
  total_days <- as.numeric(end_date - start_date + 1)
  
  # Assign precipitation values to months
  month_weights <- rep(0, 12)
  month_weights[start_month] <- start_days / total_days
  month_weights[end_month] <- end_days / total_days
  
  if (start_month != end_month) {
    month_weights[(start_month + 1):(end_month - 1)] <- days_in_month(start_month + 1:end_month - 1) / total_days
  }
  
  # Calculate weighted precipitation
  weighted_precip <- sum(precip_values * month_weights, na.rm = TRUE)
  return(weighted_precip)
}

# Apply the function across all rows
data <- data %>%
  rowwise() %>%
  mutate(Weighted_Precip = calculate_weighted_precip(Start_Date, End_Date, c_across(Precip_Jan:Precip_Dec)))

# Save results
write.csv(data, "weighted_precipitation_results.csv", row.names = FALSE)


