# Load necessary libraries. 
#Not all of these are needed at this point though
# Finalized RF models start at line 391
library(tidyverse)
library(readxl)
library(readr)
library(ranger)
library(ggplot2)
library(ggbiplot)
library(ggdist)
library(dplyr)
library(reshape2)
library(randomForest)
library(pdp)
library(gridExtra)
library(grid)
library(caret)
library(stringr)

file_path <- "C:/Users/vera780/OneDrive - PNNL/Documents/ArcGIS/Projects/Mapping K and S/CoastalTeaBagsPoints.xlsx"

# Reads the coastal tea bag excel sheet, which is filtered down from "Clean TeaBagStats_Metadata.xlsx"
# This excel sheet only includes sites that are within 100km of the coast
# Shoreline shapefile: https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/
# Shoreline filter was made by turning the shapefile polygon into a polyline,
# creating a 100km buffer around the line then using Select By Location to select tea bag studies in the buffer area,
# and then converting selected sites to the above excel file.

#-----Setting up data by filtering-----
k_data <- read_excel(file_path)
# Filter out NA for essential variables, K and Elevation
# elevation is deemed essential because a lack of elevation for a site means
# that the site is also missing all variables from rasters (example pH, Nitrogen, ect.)
k_data <- k_data %>%
  drop_na(K, Elevation_meters)
# Filter dataset for records where In_90_day_range_76_104_days is "YES"
k_data <- k_data %>%
  filter(`In_90_day_range_76_104_days` == "YES")
# Convert to tibble. This makes it easier to manipulate the data in the future
k_data <- as_tibble(k_data)

# Define the wetlands keywords and uplands pattern
wetlands <- c("Saltmarsh", "Flooded Grassland", "Oceanic raised bog", "wetlands")
uplands_pattern <- "forest|Forest"

# Find the Koppen Geiger climate classes associated with wetlands
climate_classes <- k_data %>%
  filter(ecosystem_type_reported %in% wetlands) %>%
  dplyr::select(Koppen_Geiger_Climate_Class) %>%
  drop_na() %>%
  distinct()

# Filter the dataset for wetlands and uplands within the identified climate classes
filtered_data <- k_data %>%
  filter(
    (ecosystem_type_reported %in% wetlands) |
      (str_detect(ecosystem_type_reported, uplands_pattern) &
         Koppen_Geiger_Climate_Class %in% climate_classes$Koppen_Geiger_Climate_Class)
  )

# Add a new column to categorize the ecosystem type
filtered_data <- filtered_data %>%
  mutate(Category = case_when(
    ecosystem_type_reported %in% wetlands ~ "Wetland",
    str_detect(ecosystem_type_reported, uplands_pattern) ~ "Upland",
    TRUE ~ "Other"
  ))

# Define selected columns for PCA and Random Forest modeling
selected_columns <- c("Mean_Annual_Air_Temp_C", "Mean_Annual_Precip_mm",
                      "Course_Fragments_Volumetric", "Clay", "Nitrogen",
                      "Organic_Carbon_Density", "Organic_Carbon_Stock",
                      "pH", "sand", "silt", "Soil_Organic_Carbon_Content",
                      "Water_Content_33kpa", "TotalNemotode_per100g",
                      "Bulk_Density_5_15cm",
                      "Cation_Exchange_Capacity_5_15", "Water_Content_10kpa",
                      "Water_Content_1500kpa", "K", "Elevation_meters", "Category")
# Ensure that filtered_data only includes columns specified in selected_columns
k_ready <- filtered_data %>%
  dplyr::select(all_of(selected_columns))

# Changing all of my column names individually for simplicity.
# Very long process. but needed for easier to visualize plots
print(colnames(k_ready)) # Initial names in k_ready

# making a named vector to make new variable names
new_names <- c("Mean Annual Temperature" = "Mean_Annual_Air_Temp_C", "Course Fragment Content" = "Course_Fragments_Volumetric", "Nitrogen" = "Nitrogen",
               "Organic Carbon Stock" = "Organic_Carbon_Stock", "Sand" = "sand", "Organic Carbon Content" = "Soil_Organic_Carbon_Content", 
               "Nemotodes Population Density" = "TotalNemotode_per100g", "Bulk Density" = "Bulk_Density_5_15cm", "Water Content (10 kpa)" = "Water_Content_10kpa",
               "Annual Precipitation" = "Mean_Annual_Precip_mm", "Organic Carbon Density" = "Organic_Carbon_Density", "Silt" = "silt",
               "Water Content (33kpa)" = "Water_Content_33kpa", "Cation Exchange Capacity"="Cation_Exchange_Capacity_5_15", "Elevation" ="Elevation_meters",
               "pH" = "pH", "K" = "K", "Clay" = "Clay", "Category"="Category", "Water Content (1500kpa)" = "Water_Content_1500kpa")
# Assigning new names to columns
current_names <- colnames(k_ready)
# Using named vector "new_names" to replace column names in "k_ready"
colnames(k_ready) <- names(new_names)[match(current_names, new_names)]
# make sure they look ok.
print(colnames(k_ready))
# Save separate CSV files for wetlands and uplands
# These are getting saved to your working directory. On mine, it's the documents folder.
# write_csv(filtered_data %>% filter(Category == "Wetland"), "wetlands_data.csv")
# write_csv(filtered_data %>% filter(Category == "Upland"), "uplands_data.csv")
UPLANDS <- k_ready %>% filter(Category == "Upland")
WETLANDS <- k_ready %>% filter(Category == "Wetland")

#-----Initial plots of the data--------
# Plot the data on a global map
world_map <- map_data("world")
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
  geom_point(data = filtered_data, aes(x = longitude, y = latitude, color = Category), size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Wetlands and Uplands Locations",
       x = "Longitude",
       y = "Latitude",
       color = "Category") +
  scale_color_manual(values = c("Wetland" = "blue", "Upland" = "forestgreen"))

# Boxplot of K values separated by Category aka TAI position
ggplot(filtered_data, aes(x = Category, y = K, fill = Category)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Comparison of K Values", x = "Category", y = "K Value") +
  scale_fill_manual(values = c("Wetland" = "blue", "Upland" = "forestgreen"))

# Density plot of K values separated by Category(AKA TAI position)
ggplot(filtered_data, aes(x = K, fill = Category)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(title = "Density Plot of K Values by TAI Position", x = "K Value", y = "Density") +
  scale_fill_manual(values = c("Wetland" = "blue", "Upland" = "forestgreen"))

#---Figure 5: Raincloud plots for data---
#Function repurposed from Kaizad's code 
# https://zenodo.org/records/7106554 -> Code found in filename: 3-functions-analysis.R

# Updated figure with adjustments for reduced whitespace
# but I don't think the reduced whitespace is actually showing. Not important.
#first, creating function for the plot
create_combined_plot <- function(data) {
  # Ensure Category remains a character and use it directly in ggplot
  # Category is TAI Position
  plot_data <- data
  
  # Calculate sample sizes
  sample_size <- plot_data %>%
    group_by(Category) %>%
    summarize(n = n())
  
  # Create ggplot for boxplot and density plot combined
  gg_combined <- plot_data %>%
    ggplot(aes(x = Category, y = K, fill = Category, color = Category)) +
    geom_jitter(aes(color = Category), width = 0.1, alpha = 0.7, size = 1.2) +
    geom_boxplot(aes(fill = Category), width = 0.3, color = "black", alpha = 0.6, size = 0.5,
                 outlier.colour = "black") +
    ggdist::stat_halfeye(aes(fill = Category, color = Category),
                         size = 4, alpha = 0.6,
                         position = position_nudge(x = 0.2), width = 0.3) +
    geom_text(data = sample_size, aes(x = Category, y = max(plot_data$K, na.rm = TRUE) - 0.01, label = paste("n =", n)),
              color = "black", vjust = -1, size = 6) +
    theme_classic() +
    labs(title = "Combined Boxplot and Density Plot with Data Points",
         x = "TAI Position",
         y = "K Value") +
    scale_fill_manual(values = c("Upland" = "darkred", "Wetland" = "blue")) +
    scale_color_manual(values = c("Upland" = "darkred", "Wetland" = "blue")) +
    theme(legend.position = "none",
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14),
          panel.grid.minor = element_blank())  # Remove all minor grid lines for better visual
  
  # Return the plot
  return(gg_combined)
}

# Use the function with the dataset
gg_combined_plot <- create_combined_plot(k_ready)
print(gg_combined_plot)


#------Figure 6: PCA--------------------
# There are two pca plots for the figure 6. One with points seperated into wetlands and uplands
# and one showing arrows colored by groups
#Take out the Category column from PCA calculation, 
pca_columns <- k_ready %>% dplyr::select(-Category)
# Drop rows with missing values across all pca_columns. 
# there aren't many missing values, so this doesn't affect data quantity that much
pca_columns_clean <- pca_columns %>% drop_na()
# Perform PCA on the cleaned dataset with no NA values
selected_data_clean <- k_ready %>%
  drop_na()
pca_result <- prcomp(pca_columns_clean, scale. = TRUE)
groups <- factor(selected_data_clean$Category, levels = c("Upland", "Wetland"))

p_biplot <- ggbiplot(
  pca_result,
  obs.scale = 1,
  var.scale = 1,
  groups = groups,
  ellipse = TRUE,
  circle = TRUE,
  varname.size = 0  # Set variable name size to zero to effectively hide them
) +
  scale_color_manual(
    values = c("Upland" = "darkred", "Wetland" = "blue"),
    name = "TAI Position") +
  scale_fill_manual(
    values = c("Upland" = "darkred", "Wetland" = "blue"),
    name = "TAI Position") +
  guides(
    fill = guide_legend(title = "TAI Position"),
    color = guide_legend(title = "TAI Position"),
    shape = "none") +
  theme_minimal() +
  labs(
    title = "PCA Biplot for TAI Positions",
    x = "PC1",
    y = "PC2") +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(face = "bold")
  )

print(p_biplot)

#----Figure 6, part 2---------
# arrows colored by variable class

# Create data frame and assign classes to each variable. They are in numerical order
variable_classes <- data.frame(
  variable = c("Mean Annual Temperature", "Annual Precipitation", "Course Fragment Content",
               "Clay", "Nitrogen", "Organic Carbon Density", "Organic Carbon Stock",
               "pH", "Sand", "Silt", "Organic Carbon Content", "Water Content (33kpa)",
               "Nemotodes Population Density", "Bulk Density", "Cation Exchange Capacity",
               "Water Content (10 kpa)", "Water Content (1500kpa)", "K", "Elevation"),
  class = c("Environment", "Environment", "Physical",
            "Physical", "Chemical", "Chemical", "Chemical",
            "Chemical", "Physical", "Physical", "Chemical", "Physical",
            "Biological", "Physical", "Chemical",
            "Physical", "Physical", "Chemical", "Environment"),
  number = 1:19
)


# Verify the cleaned dataset has matching dimensions with pca_result or it won't work
print(dim(pca_columns_clean))
pca_result <- prcomp(pca_columns_clean, scale. = TRUE)
print(dim(pca_result$x)) # This should be the same as pca_columns_clean
var.cor <- cor(pca_columns_clean, pca_result$x) # Calculating variable loadings

# Rescaling the loadings
scaling_factor <- sqrt(nrow(pca_columns_clean) - 1)
rescale_factor <- 0.09
var.cor <- sweep(var.cor, 2, scaling_factor * rescale_factor, '*')

# Creating the data frame for plotting arrows
var.cor_data <- as.data.frame(var.cor)
colnames(var.cor_data) <- paste0("PC", 1:ncol(var.cor_data))
var.cor_data$varname <- rownames(var.cor)

# Combinging class information to the variable names
var.cor_data <- var.cor_data %>%
  left_join(variable_classes, by = c("varname" = "variable"))

class_colors <- c("Environment" = "black", "Physical" = "orange3", "Chemical" = "purple",
                  "Biological" = "green2")

p_base <- ggplot() # Base ggplot setup without points
# Custom arrows for variable loadings, colored by class
p_arrows <- p_base +
  geom_segment(data = var.cor_data, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2, color = class),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 1.5, #thickness of arrow
               show.legend = TRUE) +  # Show color arrows in legend
  # Add numbered text labels to the arrows if you want. I removed it.
  #geom_text(data = var.cor_data,
   #         aes(x = PC1, y = PC2, label = number),
    #        vjust = 0.75, hjust = 0.0, size = 4.5,
      #      show.legend = FALSE) +
  scale_color_manual(values = class_colors, name = "Variable Class") +
  labs(title = "PCA Biplot with Classified Variables",
       x = paste("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(face = "bold"))
# Plot with minimal theme and only arrows
print(p_arrows)

#-----Seperation--------
#everything below is for random forest. 
#Don't rename any datasets above or below at the same time to make sure stuff works

#Following along on Random Forest practice using IKEA furniture
# Tutorial is made by , Julia Silge, a data scientist at RStudio
#Video name: Get started with random forest tuning and tidymodels using IKEA price data
#Video link: https://www.youtube.com/watch?v=BgWCuyrwD1s 

#-----AMP's RF edits----------------

#  title: "teabag RFs"
#author: "AMP"
  
library(tidyverse)

#1) loading data and getting it in the right format (code from Raul's other script)
file_path <- "C:/Users/vera780/OneDrive - PNNL/Documents/ArcGIS/Projects/Mapping K and S/CoastalTeaBagsPoints.xlsx"


# Reads the coastal tea bag excel sheet, which is filtered down from "Clean TeaBagStats_Metadata.xlsx"
# This excel sheet only includes sites that are within 100km of the coast
# Shoreline shapefile: https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/
# Shoreline filter was made by turning the shapefile polygon into a polyline,
# creating a 100km buffer around the line then using Select By Location to select tea bag studies in the buffer area,
# and then converting selected sites to the above excel file.

#-----Setting up data by filtering-----
k_data <- readxl::read_excel(file_path)
# Filter out NA for essential variables, K and Elevation
# elevation is deemed essential because a lack of elevation for a site means
# that the site is also missing all variables from rasters (example pH, Nitrogen, ect.)
k_data <- k_data %>%
  drop_na(K, Elevation_meters)
# Filter dataset for records where In_90_day_range_76_104_days is "YES"
k_data <- k_data %>%
  filter(`In_90_day_range_76_104_days` == "YES")
# Convert to tibble. This makes it easier to manipulate the data in the future
k_data <- as_tibble(k_data)

# Define the wetlands keywords and uplands pattern
wetlands <- c("Saltmarsh", "Flooded Grassland", "Oceanic raised bog", "wetlands")
uplands_pattern <- "forest|Forest"

# Find the Koppen Geiger climate classes associated with wetlands
climate_classes <- k_data %>%
  filter(ecosystem_type_reported %in% wetlands) %>%
  dplyr::select(Koppen_Geiger_Climate_Class) %>%
  drop_na() %>%
  distinct()

# Filter the dataset for wetlands and uplands within the identified climate classes
filtered_data <- k_data %>%
  filter(
    (ecosystem_type_reported %in% wetlands) |
      (str_detect(ecosystem_type_reported, uplands_pattern) &
         Koppen_Geiger_Climate_Class %in% climate_classes$Koppen_Geiger_Climate_Class)
  )

# Add a new column to categorize the ecosystem type
filtered_data <- filtered_data %>%
  mutate(Category = case_when(
    ecosystem_type_reported %in% wetlands ~ "Wetland",
    str_detect(ecosystem_type_reported, uplands_pattern) ~ "Upland",
    TRUE ~ "Other"
  ))

# Define selected columns for PCA and Random Forest modeling
selected_columns <- c("Mean_Annual_Air_Temp_C", "Mean_Annual_Precip_mm",
                      "Course_Fragments_Volumetric", "Clay", "Nitrogen",
                      "Organic_Carbon_Density", "Organic_Carbon_Stock",
                      "pH", "sand", "silt", "Soil_Organic_Carbon_Content",
                      "Water_Content_33kpa", "TotalNemotode_per100g",
                      "Bulk_Density_5_15cm",
                      "Cation_Exchange_Capacity_5_15", "Water_Content_10kpa",
                      "Water_Content_1500kpa", "K", "Elevation_meters", "Category")
# Ensure that filtered_data only includes columns specified in selected_columns
k_ready <- filtered_data %>%
  dplyr::select(all_of(selected_columns))

#get data ready for RF
selected_data_clean <- k_ready %>%
  drop_na()

head(selected_data_clean)

#-----r build inital RF--------
#Build model: 
set.seed(123)
#split into training and test data frames: 
k_split <- rsample::initial_split(selected_data_clean, strata = K)
#strata is our outcome. What we want to predict
k_train <- rsample::training(k_split)
k_test <- rsample::testing(k_split)
#k_test is the training data that will be used until the end

set.seed(234)
#k_fold is the second sampling/training set up we use to validate k_test
k_folds <- rsample::bootstraps(k_train, strata = K)
k_folds

#see what the starting recipe is: 
usemodels::use_ranger(K ~ ., data = k_train)
#this gives us a recipe to use for the model based on the training dataset to help with tuning the random forest model.

#Modified output from use_ranger step above :
#modify to include normalization and correlation removal steps:
ranger_recipe <- 
  recipes::recipe(formula = K ~ ., data = k_train) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>% #need to normalize/scale data since all variables are not normalized 
  recipes::step_corr(recipes::all_numeric_predictors()) # remove co-correlating variables to reduce feature redundancy 

#include parameters necessary for downstream analysis 
ranger_spec <- 
  parsnip::rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  parsnip::set_mode("regression") %>% #modify the output to match analysis we want to do, which is a regression model, not a classification model
  parsnip::set_engine("ranger", importance= "impurity") #adding in importance metric to be able to get variable importance later. 

#Create the workflow:
ranger_workflow <- 
  workflows::workflow() %>% 
  workflows::add_recipe(ranger_recipe) %>% 
  workflows::add_model(ranger_spec) 

# Re-tune the model
# This takes my laptop # minutes to run
set.seed(59146)
ranger_tune <-
  tune::tune_grid(ranger_workflow, resamples = k_folds)

#note, you can try to speed that step up with your steps from the other code, I just let it run while I did something else. 

#now, show what the best models from that are using different metrics of model performance: 
tune::show_best(ranger_tune, metric = "rmse") #root mean square error
tune::show_best(ranger_tune, metric = "rsq")  #R^2

#finalize workflow with best performing parameters:
final_rf_wf <- ranger_workflow %>%
  tune::finalize_workflow(tune::select_best(ranger_tune, metric = "rmse"))

final_rf_wf

# directly fit a ranger model with importance}

# Clean column names in the training and testing datasets
k_train_clean <- k_train %>%
  janitor::clean_names()

k_test_clean <- k_test %>%
  janitor::clean_names()

# Retrieve the best parameters from the tuning
best_params <- tune::select_best(ranger_tune, metric = "rmse")
best_mtry <- best_params$mtry
best_min_n <- best_params$min_n

# Remove low-variance features
low_variance_threshold <- 0.01
variance_df <- k_train_clean %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "variance")

selected_features <- variance_df %>%
  filter(variance > low_variance_threshold | feature == "k") %>%
  pull(feature)

k_train_low_variance <- k_train_clean %>%
  select(all_of(selected_features))

k_test_low_variance <- k_test_clean %>%
  select(all_of(selected_features))

# Step 1: Standardize the data (mean = 0, sd = 1)
standardize <- function(x) {
  return((x - mean(x)) / sd(x))
}

k_train_standardized <- k_train_low_variance %>% mutate(across(where(is.numeric), standardize))
k_test_standardized <- k_test_low_variance %>% mutate(across(where(is.numeric), standardize))

# Step 2: Remove highly correlated features
cor_matrix <- cor(k_train_standardized %>% dplyr::select(-k))  # Exclude the target variable from correlation matrix 
highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.9, verbose = TRUE)
k_train_reduced <- k_train_standardized %>% dplyr::select(-names(highly_correlated))
k_test_reduced <- k_test_standardized %>% dplyr::select(-names(highly_correlated))

# Fit the ranger model with cleaned column names, standardized data and best parameters
ranger_model <- ranger::ranger(
  formula = k ~ ., 
  data = k_train_standardized, 
  num.trees = 1000, 
  mtry = best_mtry,
  min.node.size = best_min_n,
  importance = 'impurity'
)
# Predict on test data
test_predictions <- predict(ranger_model, data = k_test_standardized)$predictions

# Calculate performance metrics manually
metrics <- yardstick::metrics(tibble(truth = k_test_clean$k, estimate = test_predictions), truth, estimate)
print(metrics)

# Extract and print feature importance
importance <- ranger_model$variable.importance
print(importance)

total_importance <- sum(importance)
importance_percent <- (importance/total_importance) * 100

#-----r plotting-----
# Convert to data frame for plotting
importance_df <- as_tibble(importance_percent, rownames = "feature") %>%
  rename(importance = value) %>%
  arrange(desc(importance))

unique(importance_df$feature) #get feature names

# Create bar plot of feature importance
importance_plot <- ggplot(importance_df, aes(x = reorder(feature, importance), y = importance)) +
  geom_bar(stat = "identity",  fill = "purple4") +
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(x = "",
       y = "Feature Importance (%)") +
  theme_minimal() +
  scale_x_discrete(labels=c("water_content_1500kpa" = "Water Content (1500 kpa)",
                            "sand" = "Sand",
                            "silt" = "Silt",
                            "clay" = "Clay",
                            "organic_carbon_density" = "Organic Carbon Density (hectogram/cubic meter)",
                            "elevation_meters" = "Elevation (m)",
                            "bulk_density_5_15cm" = "Bulk Density (5-15cm)",
                            "mean_annual_air_temp_c" = "Mean Annual Air Temperature (C)", 
                            "mean_annual_precip_mm" = "Mean Annual Precip (mm)",
                            "organic_carbon_stock" = "Organic Carbon Stock",
                            "water_content_10kpa" = "Water Content (10 kpa)",
                            "nitrogen" = "Nitrogen",
                            "p_h" = "pH",
                            "soil_organic_carbon_content" = "Soil Organic Carbon Content",
                            "course_fragments_volumetric" = "Course Fragments Volumetric",
                            "total_nemotode_per100g" = "Total Nemotode Content (/100g)",
                            "cation_exchange_capacity_5_15" = "Cation Exchange Capacity (5-15cm)",
                            "water_content_33kpa" = "Water Content (33kpa)"
  )) 

# Display the plot
print(importance_plot)

#----------Uplands RF using AMP example--------


library(tidyverse)

#1) loading data and getting it in the right format (code from Raul's other script)
file_path <- "C:/Users/vera780/OneDrive - PNNL/Documents/ArcGIS/Projects/Mapping K and S/CoastalTeaBagsPoints.xlsx"


# Reads the coastal tea bag excel sheet, which is filtered down from "Clean TeaBagStats_Metadata.xlsx"
# This excel sheet only includes sites that are within 100km of the coast
# Shoreline shapefile: https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/
# Shoreline filter was made by turning the shapefile polygon into a polyline,
# creating a 100km buffer around the line then using Select By Location to select tea bag studies in the buffer area,
# and then converting selected sites to the above excel file.

#-----Setting up data by filtering-----
#data for uplands is already filtered in previous steps
upland_rf_data <- UPLANDS #renamed because I want to edit the uplands data
upland_rf_data <- dplyr::select(upland_rf_data, -Category) 

head(upland_rf_data)

#-----r build inital RF--------
#Build model: 
set.seed(123)
#split into training and test data frames: 
k_split <- rsample::initial_split(upland_rf_data, strata = K)
#strata is our outcome. What we want to predict
k_train <- rsample::training(k_split)
k_test <- rsample::testing(k_split)
#k_test is the training data that will be used until the end

set.seed(234)
#k_fold is the second sampling/training set up we use to validate k_test
k_folds <- rsample::bootstraps(k_train, strata = K)
k_folds

#see what the starting recipe is: 
usemodels::use_ranger(K ~ ., data = k_train)
#this gives us a recipe to use for the model based on the training dataset to help with tuning the random forest model.

#Modified output from use_ranger step above :
#modify to include normalization and correlation removal steps:
ranger_recipe <- 
  recipes::recipe(formula = K ~ ., data = k_train) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>% #need to normalize/scale data since all variables are not normalized 
  recipes::step_corr(recipes::all_numeric_predictors()) # remove co-correlating variables to reduce feature redundancy 

#include parameters necessary for downstream analysis 
ranger_spec <- 
  parsnip::rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  parsnip::set_mode("regression") %>% #modify the output to match analysis we want to do, which is a regression model, not a classification model
  parsnip::set_engine("ranger", importance= "impurity") #adding in importance metric to be able to get variable importance later. 

#Create the workflow:
ranger_workflow <- 
  workflows::workflow() %>% 
  workflows::add_recipe(ranger_recipe) %>% 
  workflows::add_model(ranger_spec) 

# Re-tune the model
set.seed(59146)
ranger_tune <-
  tune::tune_grid(ranger_workflow, resamples = k_folds)

#now, show what the best models from that are using different metrics of model performance: 
tune::show_best(ranger_tune, metric = "rmse") #root mean square error
tune::show_best(ranger_tune, metric = "rsq")  #R^2

#finalize workflow with best performing parameters:
final_rf_wf <- ranger_workflow %>%
  tune::finalize_workflow(tune::select_best(ranger_tune, metric = "rmse"))

final_rf_wf
# Clean column names in the training and testing datasets
k_train_clean <- k_train %>%
  janitor::clean_names()

k_test_clean <- k_test %>%
  janitor::clean_names()

# Retrieve the best parameters from the tuning
best_params <- tune::select_best(ranger_tune, metric = "rmse")
best_mtry <- best_params$mtry
best_min_n <- best_params$min_n

# Remove low-variance features
low_variance_threshold <- 0.01
variance_df <- k_train_clean %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "variance")

selected_features <- variance_df %>%
  filter(variance > low_variance_threshold | feature == "k") %>%
  pull(feature)

k_train_low_variance <- k_train_clean %>%
  select(all_of(selected_features))

k_test_low_variance <- k_test_clean %>%
  select(all_of(selected_features))

# Step 1: Standardize the data (mean = 0, sd = 1)
standardize <- function(x) {
  return((x - mean(x)) / sd(x))
}

k_train_standardized <- k_train_low_variance %>% mutate(across(where(is.numeric), standardize))
k_test_standardized <- k_test_low_variance %>% mutate(across(where(is.numeric), standardize))

# Step 2: Remove highly correlated features
cor_matrix <- cor(k_train_standardized %>% dplyr::select(-k))  # Exclude the target variable from correlation matrix 
highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.9, verbose = TRUE)
k_train_reduced <- k_train_standardized %>% dplyr::select(-names(highly_correlated))
k_test_reduced <- k_test_standardized %>% dplyr::select(-names(highly_correlated))

# Fit the ranger model with cleaned column names, standardized data and best parameters
ranger_model <- ranger::ranger(
  formula = k ~ ., 
  data = k_train_standardized, 
  num.trees = 1000, 
  mtry = best_mtry,
  min.node.size = best_min_n,
  importance = 'impurity'
)
# Predict on test data
test_predictions <- predict(ranger_model, data = k_test_standardized)$predictions

# Calculate performance metrics manually
metrics <- yardstick::metrics(tibble(truth = k_test_clean$k, estimate = test_predictions), truth, estimate)
print(metrics)

# Extract and print feature importance
importance <- ranger_model$variable.importance
print(importance)

total_importance <- sum(importance)
importance_percent <- (importance/total_importance) * 100

#-----r plotting-----
# Convert to data frame for plotting
importance_df <- as_tibble(importance_percent, rownames = "feature") %>%
  rename(importance = value) %>%
  arrange(desc(importance))

unique(importance_df$feature) #get feature names

# Create bar plot of feature importance
importance_plot <- ggplot(importance_df, aes(x = reorder(feature, importance), y = importance)) +
  geom_bar(stat = "identity",  fill = "darkred") +
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(x = "",
       y = "Feature Importance (%)") +
  theme_minimal() +
  scale_x_discrete(labels=c("water_content_1500kpa" = "Water Content (1500 kpa)",
                            "sand" = "Sand",
                            "silt" = "Silt",
                            "clay" = "Clay",
                            "organic_carbon_density" = "Organic Carbon Density (hectogram/cubic meter)",
                            "elevation_meters" = "Elevation (m)",
                            "bulk_density_5_15cm" = "Bulk Density (5-15cm)",
                            "mean_annual_air_temp_c" = "Mean Annual Air Temperature (C)", 
                            "mean_annual_precip_mm" = "Mean Annual Precip (mm)",
                            "organic_carbon_stock" = "Organic Carbon Stock",
                            "water_content_10kpa" = "Water Content (10 kpa)",
                            "nitrogen" = "Nitrogen",
                            "p_h" = "pH",
                            "soil_organic_carbon_content" = "Soil Organic Carbon Content",
                            "course_fragments_volumetric" = "Course Fragments Volumetric",
                            "total_nemotode_per100g" = "Total Nemotode Content (/100g)",
                            "cation_exchange_capacity_5_15" = "Cation Exchange Capacity (5-15cm)",
                            "water_content_33kpa" = "Water Content (33kpa)"
  )) 

# Display the plot
print(importance_plot)



#------WETLANDS RF!--------
#----------Uplands RF using AMP example--------


library(tidyverse)

#1) loading data and getting it in the right format (code from Raul's other script)
file_path <- "C:/Users/vera780/OneDrive - PNNL/Documents/ArcGIS/Projects/Mapping K and S/CoastalTeaBagsPoints.xlsx"


# Reads the coastal tea bag excel sheet, which is filtered down from "Clean TeaBagStats_Metadata.xlsx"
# This excel sheet only includes sites that are within 100km of the coast
# Shoreline shapefile: https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/
# Shoreline filter was made by turning the shapefile polygon into a polyline,
# creating a 100km buffer around the line then using Select By Location to select tea bag studies in the buffer area,
# and then converting selected sites to the above excel file.

#-----Setting up data by filtering-----
#data for wetlands is already filtered in previous steps
wetland_rf_data <- WETLANDS #renamed to edit the data
wetland_rf_data <- dplyr::select(wetland_rf_data, -Category) 

head(wetland_rf_data)

#-----r build inital RF--------
#Build model: 
set.seed(123)
#split into training and test data frames: 
k_split <- rsample::initial_split(wetland_rf_data, strata = K)
#strata is our outcome. What we want to predict
k_train <- rsample::training(k_split)
k_test <- rsample::testing(k_split)
#k_test is the training data that will be used until the end

set.seed(234)
#k_fold is the second sampling/training set up we use to validate k_test
k_folds <- rsample::bootstraps(k_train, strata = K)
k_folds

#see what the starting recipe is: 
usemodels::use_ranger(K ~ ., data = k_train)
#this gives us a recipe to use for the model based on the training dataset to help with tuning the random forest model.

#Modified output from use_ranger step above :
#modify to include normalization and correlation removal steps:
ranger_recipe <- 
  recipes::recipe(formula = K ~ ., data = k_train) %>%
  recipes::step_normalize(recipes::all_numeric_predictors()) %>% #need to normalize/scale data since all variables are not normalized 
  recipes::step_corr(recipes::all_numeric_predictors()) # remove co-correlating variables to reduce feature redundancy 

#include parameters necessary for downstream analysis 
ranger_spec <- 
  parsnip::rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  parsnip::set_mode("regression") %>% #modify the output to match analysis we want to do, which is a regression model, not a classification model
  parsnip::set_engine("ranger", importance= "impurity") #adding in importance metric to be able to get variable importance later. 

#Create the workflow:
ranger_workflow <- 
  workflows::workflow() %>% 
  workflows::add_recipe(ranger_recipe) %>% 
  workflows::add_model(ranger_spec) 

# Re-tune the model
set.seed(59146)
ranger_tune <-
  tune::tune_grid(ranger_workflow, resamples = k_folds)

#now, show what the best models from that are using different metrics of model performance: 
tune::show_best(ranger_tune, metric = "rmse") #root mean square error
tune::show_best(ranger_tune, metric = "rsq")  #R^2

#finalize workflow with best performing parameters:
final_rf_wf <- ranger_workflow %>%
  tune::finalize_workflow(tune::select_best(ranger_tune, metric = "rmse"))

final_rf_wf
# Clean column names in the training and testing datasets
k_train_clean <- k_train %>%
  janitor::clean_names()

k_test_clean <- k_test %>%
  janitor::clean_names()

# Retrieve the best parameters from the tuning
best_params <- tune::select_best(ranger_tune, metric = "rmse")
best_mtry <- best_params$mtry
best_min_n <- best_params$min_n

# Remove low-variance features
low_variance_threshold <- 0.001
variance_df <- k_train_clean %>%
  summarise(across(everything(), var)) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "variance")

selected_features <- variance_df %>%
  filter(variance > low_variance_threshold | feature == "k") %>%
  pull(feature)

k_train_low_variance <- k_train_clean %>%
  select(all_of(selected_features))

k_test_low_variance <- k_test_clean %>%
  select(all_of(selected_features))

# Step 1: Standardize the data (mean = 0, sd = 1)
standardize <- function(x) {
  return((x - mean(x)) / sd(x))
}

k_train_standardized <- k_train_low_variance %>% mutate(across(where(is.numeric), standardize))
k_test_standardized <- k_test_low_variance %>% mutate(across(where(is.numeric), standardize))

# Step 2: Remove highly correlated features
cor_matrix <- cor(k_train_standardized %>% dplyr::select(-k))  # Exclude the target variable from correlation matrix 
highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.9, verbose = TRUE)
k_train_reduced <- k_train_standardized %>% dplyr::select(-names(highly_correlated))
k_test_reduced <- k_test_standardized %>% dplyr::select(-names(highly_correlated))

# Fit the ranger model with cleaned column names, standardized data and best parameters
ranger_model <- ranger::ranger(
  formula = k ~ ., 
  data = k_train_standardized, 
  num.trees = 1000, 
  mtry = best_mtry,
  min.node.size = best_min_n,
  importance = 'impurity'
)
# Predict on test data
test_predictions <- predict(ranger_model, data = k_test_standardized)$predictions

# Calculate performance metrics manually
metrics <- yardstick::metrics(tibble(truth = k_test_clean$k, estimate = test_predictions), truth, estimate)
print(metrics)

# Extract and print feature importance
importance <- ranger_model$variable.importance
print(importance)

total_importance <- sum(importance)
importance_percent <- (importance/total_importance) * 100

#-----r plotting-----
# Convert to data frame for plotting
importance_df <- as_tibble(importance_percent, rownames = "feature") %>%
  rename(importance = value) %>%
  arrange(desc(importance))

unique(importance_df$feature) #get feature names

# Create bar plot of feature importance
importance_plot <- ggplot(importance_df, aes(x = reorder(feature, importance), y = importance)) +
  geom_bar(stat = "identity",  fill = "blue") +
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(x = "",
       y = "Feature Importance (%)") +
  theme_minimal() +
  scale_x_discrete(labels=c("water_content_1500kpa" = "Water Content (1500 kpa)",
                            "sand" = "Sand",
                            "silt" = "Silt",
                            "clay" = "Clay",
                            "organic_carbon_density" = "Organic Carbon Density (hectogram/cubic meter)",
                            "elevation_meters" = "Elevation (m)",
                            "bulk_density_5_15cm" = "Bulk Density (5-15cm)",
                            "mean_annual_air_temp_c" = "Mean Annual Air Temperature (C)", 
                            "mean_annual_precip_mm" = "Mean Annual Precip (mm)",
                            "organic_carbon_stock" = "Organic Carbon Stock",
                            "water_content_10kpa" = "Water Content (10 kpa)",
                            "nitrogen" = "Nitrogen",
                            "p_h" = "pH",
                            "soil_organic_carbon_content" = "Soil Organic Carbon Content",
                            "course_fragments_volumetric" = "Course Fragments Volumetric",
                            "total_nemotode_per100g" = "Total Nemotode Content (/100g)",
                            "cation_exchange_capacity_5_15" = "Cation Exchange Capacity (5-15cm)",
                            "water_content_33kpa" = "Water Content (33kpa)"
  )) 

# Display the plot
print(importance_plot)

