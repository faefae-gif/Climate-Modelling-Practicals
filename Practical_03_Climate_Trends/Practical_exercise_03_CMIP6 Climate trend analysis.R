# =============================================================================
# ZAMBIA CMIP6 CLIMATE TREND ANALYSIS - R SCRIPT
# =============================================================================

rm(list=ls()) # Cleaning the shit off
gc()

# Purpose
# -------
# Reproduces the climate-grid workflow in R/Tidyverse + terra/sf for:
#   1. Daily maximum temperature (tasmax)
#   2. Daily minimum temperature (tasmin)
#   3. Daily precipitation (pr)
#
# Historical period: 1981-2010
# Model: MIROC6
#
# Features:
#   - Loads and validates gadm41_ZMB_2.shp (Level 2 districts)
#   - Crops and masks gridded data to shapefile boundaries
#   - Converts units (Kelvin -> °C, kg m-2 s-1 -> mm/day)
#   - Annual aggregations (mean for temperature, sum for precipitation)
#   - Pixel-wise Sen's slope (per decade) and Mann-Kendall p-values
#   - Exports CSV, NetCDF, GeoTIFF, and publication-ready ggplot2 maps
#
# Required R Packages:
#   install.packages(c("tidyverse", "terra", "sf", "trend", "ncdf4", "viridis"))
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(sf)
  library(trend)
  library(ncdf4)
})

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================

WORKDIR <- "C:/cs592_practicals_2026"

TASMAX_FILE <- file.path(WORKDIR, "tasmax_day_MIROC6_historical_r1i1p1f1_gn_19810101-20101231.nc")
TASMIN_FILE <- file.path(WORKDIR, "tasmin_day_MIROC6_historical_r1i1p1f1_gn_19810101-20101231.nc")
PR_FILE     <- file.path(WORKDIR, "pr_day_MIROC6_historical_r1i1p1f1_gn_19810101-20101231.nc")

# Shapefile for Zambia Level 2 districts
SHAPEFILE_PATH <- file.path(WORKDIR, "gadm41_ZMB_shp/gadm41_ZMB_2.shp")

OUTPUT_DIR <- file.path(WORKDIR, "climate_trend_outputs_r")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

START_YEAR <- 1981
END_YEAR   <- 2010
ALPHA      <- 0.05

# Export toggles
EXPORT_ANNUAL_CSV <- FALSE
EXPORT_NETCDF     <- TRUE
EXPORT_GEOTIFF    <- TRUE


# =============================================================================
# 2. BOUNDARY LOADING AND PREPARATION
# =============================================================================

load_zambia_boundary <- function(shp_path) {
  if (!file.exists(shp_path)) {
    stop(paste("Shapefile not found at:", shp_path))
  }
  
  message("\nLoading shapefile: ", shp_path)
  boundary <- st_read(shp_path, quiet = TRUE)
  
  if (st_crs(boundary)$epsg != 4326) {
    boundary <- st_transform(boundary, 4326)
  }
  
  return(boundary)
}

# Load boundary once
zambia_sf <- load_zambia_boundary(SHAPEFILE_PATH)
zambia_vect <- vect(zambia_sf)


# =============================================================================
# 3. CORE PROCESSING FUNCTIONS
# =============================================================================

process_climate_variable <- function(file_path, var_name, var_type, boundary_vect) {
  if (!file.exists(file_path)) {
    message("\nSkipping ", var_name, ": File not found -> ", file_path)
    return(NULL)
  }
  
  message("\n", paste0("=", rep("=", 78), collapse = ""))
  message("PROCESSING: ", toupper(var_name), " (", var_type, ")")
  message(paste0("=", rep("=", 78), collapse = ""))
  
  # Load NetCDF as a terra SpatRaster
  r <- rast(file_path)
  message("Loaded SpatRaster dimensions: ", nrow(r), " rows, ", ncol(r), " cols, ", nlyr(r), " layers")
  
  # Extract time vector from layer names or time dimension
  time_vals <- time(r)
  if (is.null(time_vals)) {
    # Fallback to parsing dates from layer names if metadata time is absent
    time_vals <- as.Date(sub(".*_([0-9]{8}).*", "\\1", names(r)), format="%Y%m%d")
  }
  
  years <- as.numeric(format(time_vals, "%Y"))
  
  # Subset to target years
  year_indices <- which(years >= START_YEAR & years <= END_YEAR)
  r <- r[[year_indices]]
  years <- years[year_indices]
  
  # Unit conversion
  if (var_type == "temperature") {
    # Check if values are in Kelvin (> 200)
    sample_val <- min(values(r[[1]]), na.rm = TRUE)
    if (sample_val > 200) {
      message("Converting temperature from Kelvin to Celsius...")
      r <- r - 273.15
    }
    units_label <- "°C"
  } else if (var_type == "precipitation") {
    # CMIP6 precipitation flux is typically kg m-2 s-1 -> convert to mm/day
    sample_val <- max(values(r[[1]]), na.rm = TRUE)
    if (sample_val < 0.01) {
      message("Converting precipitation flux (kg m-2 s-1) to mm/day...")
      r <- r * 86400.0
    }
    units_label <- "mm"
  }
  
  # Spatial Masking & Cropping to Shapefile
  message("Clipping raster to gadm41_ZMB_2.shp boundaries...")
  r <- crop(r, boundary_vect)
  r <- mask(r, boundary_vect)
  
  # Annual Aggregation
  message("Computing annual aggregations...")
  if (var_type == "temperature") {
    annual_r <- tapp(r, index = years, fun = mean, na.rm = TRUE)
    names(annual_r) <- paste0("yr_", unique(years))
  } else {
    annual_r <- tapp(r, index = years, fun = sum, na.rm = TRUE)
    names(annual_r) <- paste0("yr_", unique(years))
  }
  
  unique_years <- unique(years)
  
  # Optional Annual CSV export
  if (EXPORT_ANNUAL_CSV) {
    df_annual <- as.data.frame(annual_r, xy = TRUE) %>%
      pivot_longer(cols = starts_with("yr_"), names_to = "year", values_to = "value") %>%
      mutate(year = as.numeric(sub("yr_", "", year)))
    
    csv_out <- file.path(OUTPUT_DIR, paste0(var_name, "_annual_", START_YEAR, "_", END_YEAR, ".csv"))
    write_csv(df_annual, csv_out)
    message("Annual CSV saved: ", csv_out)
  }
  
  # Trend Calculation (Sen's Slope & Mann-Kendall)
  message("Calculating Sen's slopes and Mann-Kendall p-values per pixel...")
  
  n_cells <- ncell(annual_r)
  slope_vals <- rep(NA, n_cells)
  pval_vals  <- rep(NA, n_cells)
  
  vals_matrix <- values(annual_r) # rows = cells, cols = years
  
  for (i in seq_len(n_cells)) {
    series <- vals_matrix[i, ]
    valid_idx <- which(!is.na(series))
    
    if (length(valid_idx) >= 3) {
      sub_series <- series[valid_idx]
      sub_years  <- unique_years[valid_idx]
      
      # Sen's slope calculation (scaled to per decade: * 10)
      ss_fit <- tryCatch({
        trend::sens.slope(sub_series)
      }, error = function(e) NULL)
      
      mk_fit <- tryCatch({
        trend::mk.test(sub_series)
      }, error = function(e) NULL)
      
      if (!is.null(ss_fit)) {
        slope_vals[i] <- ss_fit$estimates * 10.0
      }
      if (!is.null(mk_fit)) {
        pval_vals[i] <- mk_fit$p.value
      }
    }
  }
  
  # Construct output trend raster
  trend_r <- annual_r[[1:2]]
  names(trend_r) <- c("slope", "p_value")
  values(trend_r[["slope"]]) <- slope_vals
  values(trend_r[["p_value"]]) <- pval_vals
  
  trend_r[["significant"]] <- trend_r[["p_value"]] < ALPHA
  
  # Print Summary Statistics
  valid_slopes <- slope_vals[!is.na(slope_vals)]
  sig_count <- sum(pval_vals < ALPHA, na.rm = TRUE)
  
  message("\n--- ", toupper(var_name), " TREND SUMMARY ---")
  message("Valid grid cells: ", length(valid_slopes))
  if (length(valid_slopes) > 0) {
    message("Mean slope      : ", round(mean(valid_slopes), 4), " ", units_label, "/decade")
    message("Significant cells: ", sig_count, " (", round(100 * sig_count / length(valid_slopes), 2), "%)")
  }
  
  # Save Outputs
  # 1. CSV
  df_trend <- as.data.frame(trend_r, xy = TRUE) %>%
    rename(lon = x, lat = y)
  write_csv(df_trend, file.path(OUTPUT_DIR, paste0(var_name, "_trend_", START_YEAR, "_", END_YEAR, ".csv")))
  
  # 2. NetCDF
  if (EXPORT_NETCDF) {
    nc_out <- file.path(OUTPUT_DIR, paste0(var_name, "_trend_", START_YEAR, "_", END_YEAR, ".nc"))
    writeRaster(trend_r, nc_out, overwrite = TRUE)
    message("NetCDF saved: ", nc_out)
  }
  
  # 3. GeoTIFF
  if (EXPORT_GEOTIFF) {
    tif_out <- file.path(OUTPUT_DIR, paste0(var_name, "_trend_", START_YEAR, "_", END_YEAR, ".tif"))
    writeRaster(trend_r[["slope"]], tif_out, overwrite = TRUE)
    message("GeoTIFF saved: ", tif_out)
  }
  
  # 4. Map Generation
  generate_trend_map(trend_r, var_name, units_label, zambia_sf)
}


# =============================================================================
# 4. MAP GENERATION ROUTINE
# =============================================================================

generate_trend_map <- function(trend_r, var_name, units_label, boundary_sf) {
  map_df <- as.data.frame(trend_r, xy = TRUE) %>%
    rename(lon = x, lat = y)
  
  cmap <- if (grepl("tas", var_name)) "RdBu" else "BrBG"
  if (grepl("tas", var_name)) cmap <- paste0(cmap, "_rev") # Cool-Warm reversed for temperature
  
  p <- ggplot(map_df) +
    geom_raster(aes(x = lon, y = lat, fill = slope)) +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linewidth = 0.4) +
    geom_point(data = subset(map_df, significant == TRUE), 
               aes(x = lon, y = lat), color = "black", size = 0.3, alpha = 0.5) +
    scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(11, if(grepl("tas", var_name)) "RdBu" else "BrBG"),
      name = paste0("Trend\n(", units_label, "/dec)")
    ) +
    labs(
      title = paste0("MIROC6 Historical ", toupper(var_name), " Trend (", START_YEAR, "-", END_YEAR, ")"),
      subtitle = "Masked to gadm41_ZMB_2.shp (Zambia Districts)",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  map_out <- file.path(OUTPUT_DIR, paste0(var_name, "_trend_map_", START_YEAR, "_", END_YEAR, ".png"))
  ggsave(map_out, p, width = 9, height = 8, dpi = 300)
  message("Map saved: ", map_out)
}


# =============================================================================
# 5. EXECUTION PIPELINE
# =============================================================================

main <- function() {
  message("\n", paste0("=", rep("=", 78), collapse = ""))
  message("STARTING R CLIMATE TREND ANALYSIS PIPELINE")
  message(paste0("=", rep("=", 78), collapse = ""))
  
  # Run pipeline for each variable
  process_climate_variable(TASMAX_FILE, "tasmax", "temperature", zambia_vect)
  process_climate_variable(TASMIN_FILE, "tasmin", "temperature", zambia_vect)
  process_climate_variable(PR_FILE,     "pr",     "precipitation", zambia_vect)
  
  message("\n", paste0("=", rep("=", 78), collapse = ""))
  message("ANALYSIS COMPLETE. Outputs saved to:\n", OUTPUT_DIR)
  message(paste0("=", rep("=", 78), collapse = ""))
}

# Run execution
main()
