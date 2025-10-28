# ============================================================================
# ============================================================================
# SOCastR: Soil Organic Carbon (SOC) prediction workflow with uncertainties 
# ============================================================================
# ============================================================================
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("SOCastR: Soil Organic Carbon (SOC) prediction workflow with uncertainties\n")
cat("=" , rep("=", 79), "\n\n", sep = "")


# ============================================================================
# INITIAL SETTINGS
# ============================================================================
# Load required packages with version checking
required_packages <- c("CAST", # Spatial prediction and validation
                       "caret",# Machine learning framework
                       "classInt", # Class intervals for mapping
                       "dplyr",# Data manipulation
                       "doParallel",# Parallelization
                       "terra",# Raster data handling
                       "tidyterra", # Visualise SpatRaster objects
                       "sf",# Vector data handling
                       "randomForest",# Random Forest algorithm
                       "RColorBrewer", #  Color schemes for maps 
                       "quantregForest",# Quantile regression forests for uncertainty
                       "ggplot2",# Visualization
                       "viridis",# Color palettes
                       "gridExtra",# Multiple plots
                       "grid")# Graphics
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing required package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# MAIN FUNCTION
# ============================================================================

SOCastR <- function(working_dir,
                    input_dir,
                    output_dir,
                    samples,
                    covariates,
                    soc_column,
                    n.tile = 4,
                    model_uncertainty=TRUE,
                    distance_uncertainty=TRUE){

# Set up the environment
set.seed(42)  # For reproducibility

# Set directories
setwd(working_dir)
getwd()

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE) }

# ============================================================================
# DATA LOADING AND PREPARATION
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("DATA LOADING AND PREPARATION\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Load sample data (assuming we have SOC samples as shapefile)
soc_samples <- st_read(file.path(getwd(),input_dir,samples))
cat(paste("Loaded", nrow(soc_samples), "sample points\n"))

# Load raster covariates (GeoTIFF stack)
covariates <- rast(file.path(getwd(),input_dir,covariates))
cat(paste("Loaded", nlyr(covariates), "covariate layers\n"))
cat(paste("Raster resolution:", res(covariates), "x", res(covariates), "m\n"))

# Rename layers
names(covariates) <- paste0("COV",1:length(names(covariates)))
print("Data loaded successfully")

# ============================================================================
# FUNCTIONS
# ============================================================================
# Plotting of raster data
PlotR <- function(raster_result,
                  color_schema,
                  legend_title,
                  revers = FALSE,
                  accuracy = 0.01,
                  ClassIntervallMethod = "quantile",
                  n_classes = 9){
  
  # Calculate breaks
  raster_values <- values(raster_result, na.rm = TRUE)
  breaks_result <- classIntervals(raster_values, n = n_classes, style = ClassIntervallMethod)
  breaks_values <- breaks_result$brks
  
  # Create reclassification matrix: from, to, new_value
  rcl <- cbind(breaks_values[-length(breaks_values)], 
               breaks_values[-1], 
               1:n_classes)
  
  # Classify raster
  raster_class <- classify(raster_result, rcl, include.lowest = TRUE)
  
  # Set factor levels with meaningful labels
  breaks_labels <- sprintf(paste0("%.", abs(log10(accuracy)), "f"), breaks_values)
  class_labels <- paste0(breaks_labels[-length(breaks_labels)], 
                         " - ", 
                         breaks_labels[-1])
  
  levels(raster_class) <- data.frame(value = 1:n_classes, 
                                     class = class_labels)
  
  # Get colors
  if(n_classes <= 9){
    colors <- brewer.pal(n_classes, color_schema)
  } else {
    base <- brewer.pal(9, color_schema)
    colors <- colorRampPalette(base)(n_classes)
  }
  
  if(revers){
    colors <- rev(colors)
  }
  
  # Plot
  ggplot() +
    geom_spatraster(data = raster_class) +
    scale_fill_manual(values = colors,
                      name = legend_title,
                      na.value = "lightgrey") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.grid.major = element_line(colour = "black", linewidth = 0.2),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.text = element_text(size = 11, colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 11))
}

# Plotting point data 
SFPlotR <- function(sf_points,
                   split_col = "split",
                   legend_title = "Dataset",
                   colors = c("Train" = "#0077BB", "Test" = "#CC3311")) {
  
  ggplot() +
    geom_sf(data = sf_points,
            aes_string(color = split_col),
            size = 2,
            alpha = 0.7) +
    scale_color_manual(values = colors,
                       name = legend_title,
                       na.value = "lightgrey") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "black", linewidth = 0.2),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.text = element_text(size = 11, colour = "black"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(size = 11),
      legend.position = "right"
    ) +
    labs(title = "Training and Test Samples",
         x = "Longitude", y = "Latitude")
}

# Raster sample extraction using 9 cells 
extract_8_neighbors <- function(raster_stack, points) {
  # Get cell numbers for point locations
  cell_numbers <- cellFromXY(raster_stack, st_coordinates(points))
  
  # Initialize results data frame
  result_list <- list()
  
  for(i in seq_along(cell_numbers)) {
    if(!is.na(cell_numbers[i])) {
      # Get adjacent cells (8 neighbors)
      adj_cells <- adjacent(raster_stack, cell_numbers[i], directions = 8)
      
      # Include the focal cell
      all_cells <- c(cell_numbers[i], adj_cells[,2])
      
      # Extract values from all cells for each layer
      extracted_values <- extract(raster_stack, all_cells)
      
      # Calculate mean values across the 9 cells (focal + 8 neighbors)
      median_values <- apply(extracted_values[, ], 2, median, na.rm = TRUE)
      result_list[[i]] <- median_values
    } else {
      result_list[[i]] <- rep(NA, nlyr(raster_stack))
    }
  }
  
  # Combine results into data frame
  result_matrix <- do.call(rbind, result_list)
  result_df <- data.frame(ID = 1:nrow(result_matrix), result_matrix)
  names(result_df) <- c("ID", names(raster_stack))
  
  return(result_df)
}


predict_quantreg_raster <- function(covariates, model, quantiles, n_cores = NULL) {
  library(terra)
  library(parallel)
  library(quantregForest)
  
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  temp_cov <- tempfile(fileext = ".tif")
  temp_out <- tempfile()
  dir.create(temp_out)
  writeRaster(covariates, temp_cov, overwrite = TRUE)
  cat(sprintf("Processing %d quantiles with %d cores\n", length(quantiles), n_cores))
  
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  
  # Load necessary libraries and variables on workers
  clusterEvalQ(cl, {
    library(terra)
    library(quantregForest)
  })
  clusterExport(cl, c("temp_cov", "model", "temp_out", "quantiles"), envir = environment())
  
  start_time <- Sys.time()
  
  result_files <- parLapply(cl, seq_along(quantiles), function(i) {
    library(terra)
    q <- quantiles[i]
    cov <- rast(temp_cov)
    qrf_fun <- function(model, data) predict(model, newdata = data, what = q)
    pred <- terra::predict(cov, model, fun = qrf_fun, na.rm = TRUE)
    out_file <- file.path(temp_out, sprintf("quantile_%02d.tif", i))
    writeRaster(pred, out_file, overwrite = TRUE)
    return(out_file)
  })
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat(sprintf("Completed in %.2f minutes (%.1f sec per quantile)\n",
              elapsed, elapsed * 60 / length(quantiles)))
  results <- lapply(result_files, rast)
  names(results) <- paste0("q", gsub("\\.", "", quantiles))
  unlink(temp_cov)
  unlink(temp_out, recursive = TRUE)
  return(results)
}


# ============================================================================
# EXTRACT VALUES INCLUDING SURROUNDING 8 CELLS
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("EXTRACT (FILTERED) VALUES\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Use 3x3p pixel extraction method
  extracted_covariates <- extract_8_neighbors(covariates, soc_samples)

# Combine SOC values with extracted covariates
model_data <- cbind(
  SOC = soc_samples[[paste(soc_column)]],
  extracted_covariates[, -1]  # Remove ID column
)

print("Covariate extraction completed")
print(paste("Final dataset dimensions:", paste(dim(model_data), collapse = " x ")))

# Remove NA values
complete_idx <- complete.cases(model_data)
model_data_clean <- model_data[complete_idx, ]
soc_samples_clean <- soc_samples[complete_idx, ]

cat(paste("\nOriginal samples:", nrow(model_data), "\n"))
cat(paste("Clean samples (after NA removal):", nrow(model_data_clean), "\n"))
cat(paste("Removed samples:", nrow(model_data) - nrow(model_data_clean), "\n\n"))

# Save basic data summary
data_summary <- data.frame(
  Metric = c("Total samples", "Clean samples", "Removed samples", 
             "Number of covariates", "SOC mean", "SOC sd", "SOC min", "SOC max"),
  Value = c(nrow(model_data), nrow(model_data_clean), 
            nrow(model_data) - nrow(model_data_clean),
            nlyr(covariates),
            round(mean(model_data_clean$SOC), 2),
            round(sd(model_data_clean$SOC), 2),
            round(min(model_data_clean$SOC), 2),
            round(max(model_data_clean$SOC), 2))
)
# Export
write.csv(data_summary, file.path(getwd(),output_dir, "ExtractValues_SampleDataSummary.csv"), row.names = FALSE)


# =============================================================================
# SPATIAL DATA PARTITION: STRATIFIED TRAIN/TEST SPLIT
# =============================================================================
set.seed(42)

# 1. Create spatial folds for stratified splitting (e.g., 5 folds)
spatial_folds_obj <- CreateSpacetimeFolds(
  soc_samples_clean,
  spacevar = "geometry",
  k = 5
)
# Use first fold as test set, others as training (you can cross-validate or rotate as desired)
test_idx_spatial <- spatial_folds_obj$indexOut[[1]]
train_idx_spatial <- spatial_folds_obj$index[[1]]

# 2. Prepare datasets
train_data <- model_data_clean[train_idx_spatial, ]
test_data <- model_data_clean[test_idx_spatial, ]
train_samples <- soc_samples_clean[train_idx_spatial, ]
test_samples <- soc_samples_clean[test_idx_spatial, ]

cat(paste("Spatial stratified split using CreateSpacetimeFolds\n"))
cat(paste("Training samples:", nrow(train_data), 
          "(", round(100 * length(train_idx_spatial) / nrow(model_data_clean), 1), "%)\n"))
cat(paste("Test samples:", nrow(test_data), 
          "(", round(100 * length(test_idx_spatial) / nrow(model_data_clean), 1), "%)\n"))


# 3. Visualization of the split
soc_samples_clean_epsg4326 <- st_transform(soc_samples_clean, crs = 4326)
coords_all <- st_coordinates(soc_samples_clean_epsg4326)
visual_split <- data.frame(
  x = coords_all[, 1],
  y = coords_all[, 2],
  SOC = model_data_clean$SOC,
  split = "Test"
)
visual_split$split[train_idx_spatial] <- "Train"
visual_split$split <- factor(visual_split$split, levels = c("Train", "Test"))


# Create sf object
visual_split_sf <- st_as_sf(visual_split,
                            coords = c("x", "y"),
                            crs = 4326)

# Usage example:
visual_split_sf <- st_as_sf(visual_split,
                            coords = c("x", "y"),
                            crs = 4326)

p_spatial <- SFPlotR(visual_split_sf)

# Export
ggsave(file.path(getwd(),output_dir,"SpatialDataPartition_MapTrainTestSplit.png"), 
       p_spatial, 
       width = 5, 
       height = 6, 
       dpi = 300)


# 4. Distribution assessment
p_soc_dist <- ggplot(visual_split, aes(x = SOC, fill = split)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("Train" = "#2E86AB", "Test" = "#A23B72")) +
  labs(
    title = "SOC Distribution: Train vs Test",
    x = "SOC (g/kg)",
    y = "Density",
    fill = "Set"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Export
ggsave(file.path(getwd(),output_dir,"SpatialDataPartition_SocDistSpatialBlock.png"), 
       p_soc_dist, 
       width = 7, 
       height = 6, 
       dpi = 300)


# 5. Statistical comparison
ks_result <- ks.test(train_data$SOC, test_data$SOC)

split_stats <- data.frame(
  Metric = c(
    "N_train", "N_test",
    "Train_SOC_mean", "Train_SOC_sd", "Train_SOC_min", "Train_SOC_max",
    "Test_SOC_mean", "Test_SOC_sd", "Test_SOC_min", "Test_SOC_max",
    "SOC_mean_difference", "SOC_mean_diff_percent", "KS_statistic", "KS_p_value"
  ),
  Value = c(
    nrow(train_data), nrow(test_data), 
    round(mean(train_data$SOC), 2), round(sd(train_data$SOC), 2),
    round(min(train_data$SOC), 2), round(max(train_data$SOC), 2),
    round(mean(test_data$SOC), 2), round(sd(test_data$SOC), 2),
    round(min(test_data$SOC), 2), round(max(test_data$SOC), 2),
    round(abs(mean(train_data$SOC) - mean(test_data$SOC)), 2),
    round(abs(mean(train_data$SOC) - mean(test_data$SOC)) / mean(train_data$SOC) * 100, 2),
    round(ks_result$statistic, 4), round(ks_result$p.value, 4)
  )
)
# Export
write.csv(split_stats, file.path(getwd(),output_dir,"SpatialDataPartition_TrainTestStatisticsSpatialBlock.csv"), row.names = FALSE)

cat("Spatial stratified split visualization and statistics completed.\n\n")


# ============================================================================
# SPATIAL CROSS-VALIDATION ON TRAINING DATA
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("SPATIAL CROSS-VALIDATION SETUP (Training Data Only)\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Create spatial folds using KNNDM on training data only
detectCores()
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
spatial_folds_train <- knndm(train_samples, covariates, k = 5)
stopCluster(cl)

cat(paste("Created", length(spatial_folds_train$indx_train), "spatial CV folds\n"))

# Geodist analysis
geodist_result <- geodist(
  x = train_samples,
  modeldomain = covariates,
  cvfolds = spatial_folds_train$indx_test
)

# Plot geodist ECDF
p_geodist <- plot(geodist_result, unit = "km", stat = "ecdf") +
  ggtitle("Spatial CV Validation: Distance Distributions") +
  theme_minimal()
p_geodens <- plot(geodist_result, unit = "km") +
  ggtitle("Spatial CV Validation: Distance densities") +
  theme_minimal()

# Export
ggsave(file.path(getwd(),output_dir, "SpatialCrossValidation_GeodistEcdf.png"), 
       p_geodist, 
       width = 7, 
       height = 5, 
       dpi = 300)
ggsave(file.path(getwd(),output_dir, "SpatialCrossValidation_GeodistDensity.png"), 
       p_geodens, 
       width = 7, 
       height = 5, 
       dpi = 300)


cat("Geodist analysis completed and saved.\n\n")

# Extract and save geodist metrics
# Group by Col_two and summarize
geodist_summary <- geodist_result %>%
  group_by(what) %>%
  summarise(
    n = n(),
    min = min(dist, na.rm = TRUE),
    q25 = quantile(dist, 0.25, na.rm = TRUE),
    median = median(dist, na.rm = TRUE),
    mean = mean(dist, na.rm = TRUE),
    q75 = quantile(dist, 0.75, na.rm = TRUE),
    max = max(dist, na.rm = TRUE),
    sd = sd(dist, na.rm = TRUE))

# Export
write.csv(geodist_summary, file.path(getwd(),output_dir, "SpatialCrossValidation_GeodistMetrics.csv"), 
          row.names = FALSE)


# ============================================================================
# FORWARD FEATURE SELECTION
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("FORWARD FEATURE SELECTION\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# 1. Prepare data
predictors <- names(train_data)[names(train_data) != "SOC"]
cat(paste("Total available predictors:", length(predictors), "\n"))

cat("Preparing data for FFS...\n")
cat(paste("Number of predictors:", length(predictors), "\n"))
cat(paste("Training samples:", nrow(train_data), "\n"))

# Check for issues: Remove zero variance predictors
nzv <- nearZeroVar(train_data[, predictors], saveMetrics = TRUE)
if(any(nzv$zeroVar)) {
  cat("Removing zero-variance predictors:", 
      paste(rownames(nzv)[nzv$zeroVar], collapse = ", "), "\n")
  predictors <- predictors[!nzv$zeroVar]
}

# Check for NA
if(any(is.na(train_data[, predictors]))) {
  stop("NA values detected in predictors!")
}

# Remove any list columns
train_data <- as.data.frame(train_data)
train_data <- train_data[, sapply(train_data, function(x) !is.list(x))]

# Set up training control with spatial CV
train_control <- trainControl(
  method = "cv",
  index = spatial_folds_train$indx_train,
  indexOut = spatial_folds_train$indx_test,
  savePredictions = "final",
  returnResamp = "final",
  allowParallel = FALSE)

# Perform forward feature selection
cat("Running forward feature selection\n")
ffs_model <- ffs(
  predictors = train_data[, predictors, drop = FALSE],
  response = train_data$SOC,
  method = "rf",
  metric = "RMSE",
  tuneLength = 3,
  trControl = train_control,
  ntree = 100,
  verbose = FALSE
)

cat("\nForward feature selection completed!\n")
cat(paste("Selected", length(ffs_model$selectedvars), "variables out of", 
          length(predictors), "\n"))
cat("Selected variables:\n")
cat(paste("  -", ffs_model$selectedvars, collapse = "\n"))
cat("\n\n")

# Save FFS results
ffs_results <- data.frame(
  Variable = ffs_model$selectedvars,
  Selection_order = 1:length(ffs_model$selectedvars)
)

# Export
write.csv(ffs_results, file.path(getwd(),output_dir, "ForwardFeatureSelection_SelectedVariables.csv"), 
          row.names = FALSE)


# ============================================================================
# FINAL MODEL TRAINING ON TRAINING DATA
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("FINAL MODEL TRAINING\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Train final model with selected features
final_model <- train(
  x = train_data[, ffs_model$selectedvars, drop = FALSE],
  y = train_data$SOC,
  method = "rf",
  tuneLength = 3,
  trControl = train_control,
  ntree = 500,
  importance = TRUE
)

cat("Final model training completed.\n")
cat("Training performance (Spatial CV):\n")
print(final_model$results)
cat("\n")

# Variable importance
var_importance <- varImp(final_model)
var_imp_df <- data.frame(
  Variable = rownames(var_importance$importance),
  Importance = var_importance$importance$Overall
)
var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]


# Export
write.csv(var_imp_df, file.path(getwd(),output_dir, "FinalModel_VariableImportance.csv"), 
          row.names = FALSE)
write.csv(data.frame(final_model$results), file.path(getwd(),output_dir, "FinalModel_Accuracy.csv"), 
          row.names = FALSE)


# Plot variable importance
p_varimp <- ggplot(var_imp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "#2E86AB") +
  coord_flip() +
  labs(
    title = "Variable Importance (Random Forest)",
    x = "Variable",
    y = "Importance"
  ) +
  theme_minimal()

# Export
ggsave(file.path(getwd(),output_dir,"ForwardFeatureSelection_VariableImportance.png"), 
       p_varimp, 
       width = 5, 
       height = 5, 
       dpi = 300)


# ============================================================================
# EXTERNAL VALIDATION
# ============================================================================

cat("=" , rep("=", 79), "\n", sep = "")
cat("EXTERNAL VALIDATION\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Predict on test set
test_predictions <- predict(final_model, 
                            newdata = test_data[, ffs_model$selectedvars])

# Calculate test set performance metrics
test_performance <- data.frame(
  Dataset = "Test",
  RMSE = RMSE(test_predictions, test_data$SOC),
  MAE = MAE(test_predictions, test_data$SOC),
  R2 = cor(test_predictions, test_data$SOC)^2,
  Bias = mean(test_predictions - test_data$SOC),
  RMSE_percent = RMSE(test_predictions, test_data$SOC) / mean(test_data$SOC) * 100
)

# Get training CV performance
cv_performance <- data.frame(
  Dataset = "Training_CV",
  RMSE = final_model$results$RMSE[which.min(final_model$results$RMSE)],
  MAE = final_model$results$MAE[which.min(final_model$results$RMSE)],
  R2 = final_model$results$Rsquared[which.min(final_model$results$RMSE)],
  Bias = 0,  # Not directly available from CV
  RMSE_percent = final_model$results$RMSE[which.min(final_model$results$RMSE)] / 
    mean(train_data$SOC) * 100
)

# Combine performance metrics
performance_comparison <- rbind(cv_performance, test_performance)
performance_comparison$Difference <- c(0, test_performance$RMSE - cv_performance$RMSE)

# Export
write.csv(performance_comparison, 
          file.path(getwd(),output_dir, "Validation_PerformanceComparison.csv"), 
          row.names = FALSE)

cat("=== PERFORMANCE COMPARISON ===\n")
print(performance_comparison)
cat("\n")

# Create validation scatter plot
validation_df <- data.frame(
  Observed = test_data$SOC,
  Predicted = test_predictions)

p_validation <- ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6, color = "#2E86AB", size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = "Independent Test Set Validation",
    subtitle = sprintf("R² = %.3f | RMSE = %.2f g/kg | MAE = %.2f g/kg | Bias = %.2f g/kg", 
                       test_performance$R2, test_performance$RMSE, 
                       test_performance$MAE, test_performance$Bias),
    x = "Observed SOC (g/kg)",
    y = "Predicted SOC (g/kg)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Export
ggsave(file.path(getwd(),output_dir, "Validation_ScatterPlot.png"), 
       p_validation, 
       width = 6, 
       height = 6, 
       dpi = 300)

# Residual plot
validation_df$Residual <- validation_df$Predicted - validation_df$Observed

p_residuals <- ggplot(validation_df, aes(x = Predicted, y = Residual)) +
  geom_point(alpha = 0.6, color = "#A23B72", size = 3) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", color = "blue", se = TRUE) +
  labs(
    title = "Residual Plot (Test Set)",
    x = "Predicted SOC (g/kg)",
    y = "Residuals (Predicted - Observed)"
  ) +
  theme_minimal()

# Export
ggsave(file.path(getwd(),output_dir, "Validation_ResidualPlot.png"), 
       p_residuals, 
       width = 6, 
       height = 6, 
       dpi = 300)

cat("Test set validation completed and visualizations saved.\n\n")

# ============================================================================
# ALL DATA PREDICTION: RANDOM FOREST
# ============================================================================
cat("=" , rep("=", 79), "\n", sep = "")
cat("RETRAIN ON ALL DATA FOR FINAL PREDICTION\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Stop all parallel processing
stopImplicitCluster()
registerDoSEQ()  # Register sequential back end

# Close any open connections
closeAllConnections()

# Create spatial folds for full sample data set
detectCores()
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
spatial_folds_full <- knndm(soc_samples_clean, covariates, k = 5)
stopCluster(cl)

# Train final model on all data with PARALLEL DISABLED
final_model_full <- train(
  x = model_data_clean[, ffs_model$selectedvars, drop = FALSE],
  y = model_data_clean$SOC,
  method = "rf",
  tuneLength = 3,
  trControl = trainControl(
    method = "cv",
    index = spatial_folds_full$indx_train,
    indexOut = spatial_folds_full$indx_test,
    savePredictions = "final",
    allowParallel = FALSE  
  ),
  ntree = 500,
  importance = TRUE
)

# Prediction
cat("=" , rep("=", 79), "\n", sep = "")
cat("SPATIAL PREDICTION\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Create SOC prediction map
soc_prediction <- terra::predict(covariates, final_model_full, na.rm = TRUE)
names(soc_prediction) <- "SOC_predicted"

cat("SOC prediction map created.\n")

# Export prediction raster
writeRaster(soc_prediction, 
            file.path(getwd(),output_dir, "FinalPrediction_SocRaster.tif"), 
            overwrite = TRUE)

cat("Prediction map saved.\n\n")

# Export 
plot_prediction <- PlotR(soc_prediction,"YlOrBr","SOC [%]")
ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapSocRF.png"), 
       plot_prediction, 
       width = 6, 
       height = 7, dpi = 300)


# ============================================================================
# MODEL-BASED UNCERTAINTY: PREDICTION INTERVALLS
# ============================================================================# ============================================================================

if(model_uncertainty==TRUE){
cat("=", rep("=", 79), "\n", sep = "")
cat("MODEL UNCERTAINTY: QUANTILE REGRESSION FOREST UNCERTAINTY LAYERS\n")
cat("=", rep("=", 79), "\n\n", sep = "")

# Prepare training data for QRF
train_x <- model_data_clean[, ffs_model$selectedvars, drop = FALSE]
train_y <- model_data_clean$SOC

# Fit quantile regression forest model
cat("Training quantile regression forest model...\n")
qrf_model <- quantregForest(
  x = train_x,
  y = train_y,
  ntree = 500,
  importance = TRUE
)

# Calculate accuracy metrics for median prediction
qrf_model.lm <- lm(train_y ~ qrf_model$predicted)
qrf_model.lm.R2 <- summary(qrf_model.lm)$r.squared
qrf_model.lm.RMSE <- sqrt(mean((train_y - qrf_model$predicted)^2))

# Export training results
qrf_model.accuracy <- data.frame(
  R_squared = qrf_model.lm.R2,
  RMSE = qrf_model.lm.RMSE,
  Method = "Quantile_Regression_Forest"
)

# Export
write.csv(qrf_model.accuracy, 
          file.path(getwd(), output_dir, "QuantilePrediction_Accuracy.csv"), 
          row.names = FALSE)

cat("QRF model training completed.\n\n")

# Define quantiles for prediction intervals
quantiles <- c(0.05, 0.5, 0.95)

# Predict quantiles on raster stack - CORRECTED VERSION
cat("Generating quantile predictions across spatial domain...\n")

# Create separate predictions for each quantile
pb <- txtProgressBar(min = 0, max = length(quantiles), style = 3)

for(i in seq_along(quantiles)){
  # Assign prediction results dynamically using assign()
  assign(paste0("soc_q", quantiles[i]*100), 
         terra::predict(
           covariates,
           qrf_model,
           na.rm = TRUE,
           fun = function(model, data) {
             predict(model, data, what = quantiles[i])
           }
         )
  )
  # Update progress bar
  setTxtProgressBar(pb, i)
}
# Close the progress bar
close(pb)


#ncores <- detectCores() - 1
#quantile_results <- predict_quantreg_raster(
#  covariates = covariates,
#  model = qrf_model,
#  quantiles = c(0.05, 0.5, 0.95),
#  n_cores = ncores
#)
#soc_q05 <- quantile_results$q005
#soc_q50 <- quantile_results$q05
#soc_q95 <- quantile_results$q095

# Combine into one SpatRaster stack
soc_qrf_pred <- c(soc_q5, soc_q50, soc_q95)
names(soc_qrf_pred) <- c("SOC_q05", "SOC_q50", "SOC_q95")


cat("Quantile predictions completed.\n\n")

# Calculate prediction interval width (95% - 5%)
soc_qrf_interval_width <- soc_qrf_pred$SOC_q95 - soc_qrf_pred$SOC_q05
names(soc_qrf_interval_width) <- "SOC_PI_width"

# Combine all layers
soc_qrf_stack <- c(soc_qrf_pred, soc_qrf_interval_width)

# Export
writeRaster(soc_qrf_stack,
            file.path(getwd(), output_dir, "FinalPrediction_SocQuantileLayers.tif"),
            overwrite = TRUE)

cat("Quantile prediction raster layers saved.\n")

# Summary statistics
pi_summary <- data.frame(
  Metric = c("Mean_PI_width", "Median_PI_width", "95th_percentile_width"),
  Value = c(
    round(mean(values(soc_qrf_interval_width), na.rm = TRUE), 3),
    round(median(values(soc_qrf_interval_width), na.rm = TRUE), 3),
    round(quantile(values(soc_qrf_interval_width), 0.95, na.rm = TRUE), 3)
  )
)
write.csv(pi_summary, 
          file.path(getwd(), output_dir, "QuantilePrediction_UncertaintySummary.csv"), 
          row.names = FALSE)

cat("Quantile-based uncertainty summary statistics saved.\n\n")

# Export

plot_prediction <- PlotR(soc_qrf_pred$SOC_q05,"YlOrBr","SOC [%]")
ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapQuantialSoc5Percentile.png"), 
       plot_prediction, 
       width = 6, 
       height = 7, dpi = 300)

plot_prediction <- PlotR(soc_qrf_pred$SOC_q50,"YlOrBr","SOC [%]")
ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapQuantialSoc50Percentile.png"), 
       plot_prediction, 
       width = 6, 
       height = 7, dpi = 300)

plot_prediction <- PlotR(soc_qrf_pred$SOC_q95,"YlOrBr","SOC [%]")
ggsave(file.path(getwd(),output_dir,"FinalPrediction_QuantialMapSoc95Percentile.png"), 
       plot_prediction, 
       width = 6, 
       height = 7, dpi = 300)

plot_PIW <- PlotR(soc_qrf_interval_width,"Spectral","PIW (90%)",revers=TRUE)
ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapPIW90.png"), 
       plot_PIW, 
       width = 6, 
       height = 7, dpi = 300)


cat("Quantile prediction maps saved.\n")
}

# ============================================================================
# DISTANCE-BASED UNCERTAINTY QUANTIFICATION: DISSIMILARITY INDEX
# ============================================================================

if(distance_uncertainty==TRUE){
cat("=" , rep("=", 79), "\n", sep = "")
cat("DISTANCE UNCERTAINTY QUANTIFICATION\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# Calculate Area of Applicability
cat("Calculating Area of Applicability with Tile-based processing ...\n")

# 1: Pre-compute trainDI
cat("Step 1: Computing trainDI threshold...\n")
trainDI_result <- trainDI(final_model_full)
saveRDS(trainDI_result, paste0(getwd(),output_dir, "trainDI_object.rds"))
cat(paste("TrainDI threshold:", round(trainDI_result$threshold, 4), "\n\n"))

# 2: Check if tiling is needed
raster_cells <- ncell(covariates)
cat(paste("Total raster cells:", format(raster_cells, big.mark = ","), "\n"))

if(raster_cells > 5e6) {
  cat("\nLarge raster detected. Using tile-based processing...\n\n")
  
  # Create tiles directory
  dir.create("tiles", showWarnings = FALSE)
  dir.create("aoa_tiles", showWarnings = FALSE)
  
  # Manual tile creation function
  create_tiles_manual <- function(raster, n_tiles_x = 4, n_tiles_y = 4, 
                                  output_dir = "tiles") {
    
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    ext_full <- ext(raster)
    x_min <- ext_full[1]
    x_max <- ext_full[2]
    y_min <- ext_full[3]
    y_max <- ext_full[4]
    
    tile_width <- (x_max - x_min) / n_tiles_x
    tile_height <- (y_max - y_min) / n_tiles_y
    
    cat(paste("Creating", n_tiles_x * n_tiles_y, "tiles\n"))
    
    tile_files <- character()
    
    for(i in 1:n_tiles_y) {
      for(j in 1:n_tiles_x) {
        
        tile_xmin <- x_min + (j - 1) * tile_width
        tile_xmax <- x_min + j * tile_width
        tile_ymin <- y_min + (i - 1) * tile_height
        tile_ymax <- y_min + i * tile_height
        
        tile_extent <- ext(tile_xmin, tile_xmax, tile_ymin, tile_ymax)
        tile_raster <- crop(raster, tile_extent)
        
        tile_filename <- file.path(output_dir, 
                                   sprintf("tile_%02d_%02d.tif", i, j))
        writeRaster(tile_raster, tile_filename, overwrite = TRUE)
        
        tile_files <- c(tile_files, tile_filename)
      }
    }
    
    cat(paste("Created", length(tile_files), "tiles\n"))
    return(tile_files)
  }
  
  # Create tiles
  cat("Creating raster tiles...\n")
  tile_files <- create_tiles_manual(
    raster = covariates,
    n_tiles_x = n.tile,
    n_tiles_y = n.tile,
    output_dir = "tiles"
  )
  
  cat("\n")
  
  # Detect available cores
  ncores <- detectCores() - 1
  cat(paste("Using", ncores, "cores for parallel processing\n\n"))
  
  # Process tiles in parallel
  cat("Processing AOA for each tile...\n")
  start_time <- Sys.time()
  
  # For Unix/Linux/Mac
  if(.Platform$OS.type == "unix") {
    
    tile_results <- mclapply(seq_along(tile_files), function(i) {
      
      tile_path <- tile_files[i]
      
      # Load tile
      tile_raster <- rast(tile_path)
      
      # Calculate AOA for tile
      aoa_tile <- aoa(
        newdata = tile_raster,
        trainDI = trainDI_result,
        LPD = FALSE,
        verbose = FALSE
      )
      
      # Save results
      tile_basename <- basename(tile_path)
      tile_name <- sub("\\.tif$", "", tile_basename)
      
      di_file <- file.path("aoa_tiles", paste0("DI_", tile_name, ".tif"))
      aoa_file <- file.path("aoa_tiles", paste0("AOA_", tile_name, ".tif"))
      
      writeRaster(aoa_tile$DI, di_file, overwrite = TRUE)
      writeRaster(aoa_tile$AOA, aoa_file, overwrite = TRUE)
      
      return(list(di = di_file, aoa = aoa_file))
      
    }, mc.cores = ncores)
    
  } else {
    # For Windows
    library(doParallel)
    library(foreach)
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    clusterExport(cl, c("trainDI_result"))
    clusterEvalQ(cl, {
      library(CAST)
      library(terra)
    })
    
    tile_results <- foreach(i = seq_along(tile_files), 
                            .packages = c("CAST", "terra")) %dopar% {
                              
                              tile_path <- tile_files[i]
                              tile_raster <- rast(tile_path)
                              
                              aoa_tile <- aoa(
                                newdata = tile_raster,
                                trainDI = trainDI_result,
                                LPD = FALSE,
                                verbose = FALSE
                              )
                              
                              tile_basename <- basename(tile_path)
                              tile_name <- sub("\\.tif$", "", tile_basename)
                              
                              di_file <- file.path("aoa_tiles", paste0("DI_", tile_name, ".tif"))
                              aoa_file <- file.path("aoa_tiles", paste0("AOA_", tile_name, ".tif"))
                              
                              writeRaster(aoa_tile$DI, di_file, overwrite = TRUE)
                              writeRaster(aoa_tile$AOA, aoa_file, overwrite = TRUE)
                              
                              list(di = di_file, aoa = aoa_file)
                            }
    
    stopCluster(cl)
  }
  
  end_time <- Sys.time()
  processing_time <- difftime(end_time, start_time, units = "mins")
  
  cat(paste("\nTile processing completed in", 
            round(processing_time, 1), "minutes\n\n"))
  
  # Merge tiles
  cat("Merging tiles into final rasters...\n")
  
  di_files <- list.files("aoa_tiles", pattern = "^DI_.*\\.tif$", 
                         full.names = TRUE)
  aoa_files <- list.files("aoa_tiles", pattern = "^AOA_.*\\.tif$", 
                          full.names = TRUE)
  
  # Create mosaic (terra::mosaic is more robust than merge)
  di_rasters <- lapply(di_files, rast)
  aoa_rasters <- lapply(aoa_files, rast)
  
  di_merged <- do.call(mosaic, di_rasters)
  aoa_merged <- do.call(mosaic, aoa_rasters)
  
  # Export
  writeRaster(di_merged, 
              file.path(getwd(),output_dir, "FinalPrediction_DissimilarityIndex.tif"), 
              overwrite = TRUE)
  writeRaster(aoa_merged, 
              file.path(getwd(),output_dir, "FinalPrediction_AOAmask.tif"), 
              overwrite = TRUE)
  
  aoa_result <- list(DI = di_merged, AOA = aoa_merged)
  
  cat("Tile merging completed.\n\n")
  
  # Optional: Clean up temporary files
  # unlink("tiles", recursive = TRUE)
  # unlink("aoa_tiles", recursive = TRUE)
  
} else {
  
  # Small raster - process directly
  cat("Processing AOA directly (no tiling needed)...\n")
  
  aoa_result <- aoa(
    newdata = covariates,
    trainDI = trainDI_result,
    LPD = FALSE,
    verbose = TRUE
  )
  
  # Export
  writeRaster(di_merged, 
              file.path(getwd(),output_dir, "FinalPrediction_DissimilarityIndex.tif"), 
              overwrite = TRUE)
  writeRaster(aoa_merged, 
              file.path(getwd(),output_dir, "FinalPrediction_AOAmask.tif"), 
              overwrite = TRUE)
}

# Calculate summary statistics
cat("Calculating AOA summary statistics...\n")

aoa_summary <- data.frame(
  Metric = c("Total_prediction_cells", "Cells_within_AOA", 
             "Percent_within_AOA", "Mean_DI", "Max_DI"),
  Value = c(
    sum(!is.na(values(aoa_result$AOA))),
    sum(values(aoa_result$AOA) == 1, na.rm = TRUE),
    round(100 * sum(values(aoa_result$AOA) == 1, na.rm = TRUE) / 
            sum(!is.na(values(aoa_result$AOA))), 2),
    round(mean(values(aoa_result$DI), na.rm = TRUE), 4),
    round(max(values(aoa_result$DI), na.rm = TRUE), 4)
  )
)

# Export
write.csv(aoa_summary, file.path(getwd(),output_dir, "FinalPrediction_AoaSummaryStatistics.csv"), 
          row.names = FALSE)


cat("\n=== AOA COMPUTATION COMPLETED ===\n")
print(aoa_summary)
cat("\n")

plot_prediction <- PlotR(di_merged,
                         "Spectral",
                         "DI",
                         revers = TRUE,
                         accuracy = 0.01,
                         ClassIntervallMethod = "quantile",
                         n_classes = 9)

ggsave(file.path(getwd(),output_dir,"FinalPrediction_MapDI.png"), 
       plot_prediction, 
       width = 6, 
       height = 7, dpi = 300)

}
}














