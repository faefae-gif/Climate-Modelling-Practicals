# ============================================================
# PRACTICAL EXERCISE 3
# SPATIAL CLIMATE TREND ANALYSIS USING CMIP6 DATA
# ============================================================

# ------------------------------------------------------------
# 0. SETUP
# install.packages(c("tidyverse", "terra", "sf", "trend", "ncdf4", "viridis"))
# ------------------------------------------------------------

rm(list =ls())
gc()

library(terra)
library(tidyverse)
library(sf)
library(trend)
library(ncdf4)
library(viridis)

setwd("C:/Users/User/Desktop/Climate_Practical_3")
list.files()
getwd()

# ------------------------------------------------------------
# TASK 1: DATA EXPLORATION AND UNIT STANDARDIZATION
# ------------------------------------------------------------

# 1A. Load tasmax

tasmax <- rast("tasmax_day_MIROC6_historical_r1i1p1f1_gn_19810101-20101231.nc")
#Inspect
tasmax
crs(tasmax)
ext(tasmax)
time(tasmax)

# Check range
global(tasmax[[1]], range, na.rm =TRUE)

# Convert Kelvin to Celsius
tasmax_c <- tasmax - 273.15

global(tasmax_c[[1]], range, na.rm =TRUE)

# 1B. Load tasmin
tasmin <- rast("tasmin_day_MIROC6_historical_r1i1p1f1_gn_19810101-20101231.nc")
#Inspect
tasmin
crs(tasmin)
ext(tasmin)
time(tasmin)

# Check range
global(tasmin[[1]], range, na.rm =TRUE)

# Convert Kelvin to Celsius
tasmin_c <- tasmin - 273.15

global(tasmin_c[[1]], range, na.rm =TRUE)

# 1C. Load precipitation
pr <- rast("pr_day_MIROC6_historical_r1i1p1f1_gn_19810101-20101231.nc")

# Inspect
pr
crs(pr)
ext(pr)
time(pr)

# Check range
global(pr[[1]], range, na.rm =TRUE)

# Convert precipitation to mm/day
pr_mmday <- pr * 86400

global(pr_mmday[[1]], range, na.rm =TRUE)

# ------------------------------------------------------------
# TASK 2: SPATIAL MASKING AND CROPPING
# ------------------------------------------------------------

# 2A. Load Zambia Level 2 administrative boundaries
zambia <- vect("gadm41_ZMB_2.shp")

#Inspect boundary
zambia
crs(zambia)
ext(zambia)

# 2B. Check CRS compatibility with climate data

crs(tasmax_c)
crs(tasmin_c)
crs(pr_mmday)

same.crs(tasmax_c, zambia)
same.crs(tasmin_c, zambia)
same.crs(pr_mmday, zambia)

# 2C. Crop climate rasters to Zambia bounding box
tasmax_crop <- crop(tasmax_c, zambia)
tasmin_crop <- crop(tasmin_c, zambia)
pr_crop <- crop(pr_mmday, zambia)

# Inspect cropped extents
ext(tasmax_crop)
ext(tasmin_crop)
ext(pr_crop)
ext(zambia)

# 2D. Mask cells outside Zambia boundaries

tasmax_mask <- mask(tasmax_crop, zambia)
tasmin_mask <- mask(tasmin_crop, zambia)
pr_mask     <- mask(pr_crop, zambia)

# Inspect masked rasters
tasmax_mask
tasmin_mask
pr_mask


# 2E. Visual check

plot(
  tasmax_mask[[1]],
  main = "Spatial Masking of Daily Maximum Temperature over Zambia"
)
lines(zambia)

plot(
  tasmin_mask[[1]],
  main = "Spatial Masking of Daily Minimum Temperature over Zambia"
)
lines(zambia)

plot(
  pr_mask[[1]],
  main = "Spatial Masking of Daily Precipitation over Zambia"
)
lines(zambia)

# ------------------------------------------------------------
# TASK 3: TEMPORAL AGGREGATION
# ------------------------------------------------------------

# 3A. Create year index from daily time dimension
years <- format(time(tasmax_mask), "%Y")

# Check year index
head(years)
tail(years)
length(years)

# 3B. Count valid daily observations for each year
tasmax_count <- tapp(
  !is.na(tasmax_mask),
  years,
  fun = sum
)

tasmin_count <- tapp(
  !is.na(tasmin_mask),
  years,
  fun = sum
)

pr_count <- tapp(
  !is.na(pr_mask),
  years,
  fun = sum
)

# 3C. Calculate annual mean temperatures

tasmax_mean <- tapp(
  tasmax_mask,
  years,
  fun = mean,
  na.rm = TRUE
)

tasmin_mean <- tapp(
  tasmin_mask,
  years,
  fun = mean,
  na.rm = TRUE
)

# 3D. Calculate annual precipitation totals

pr_total <- tapp(
  pr_mask,
  years,
  fun = sum,
  na.rm = TRUE
)

# 3E. Apply completeness threshold
# Require at least 300 valid daily observations per year

tasmax_annual <- ifel(
  tasmax_count >= 300,
  tasmax_mean,
  NA
)

tasmin_annual <- ifel(
  tasmin_count >= 300,
  tasmin_mean,
  NA
)

pr_annual <- ifel(
  pr_count >= 300,
  pr_total,
  NA
)

# 3F. Check results
nlyr(tasmax_annual)
nlyr(tasmin_annual)
nlyr(pr_annual)

tasmax_annual
tasmin_annual
pr_annual

# ------------------------------------------------------------
# TASK 4: TREND AND SIGNIFICANCE DETECTION
# ------------------------------------------------------------


# ------------------------------------------------------------
# 4A. Function to calculate Sen's slope per decade
# ------------------------------------------------------------

sen_fun <- function(x) {
  
  # Remove missing values
  x <- x[!is.na(x)]
  
  # Require enough valid annual values
  if (length(x) < 10) {
    return(NA)
  }
  
  # Calculate Sen's slope
  slope <- sens.slope(x)$estimates
  
  # Convert from per year to per decade
  slope_decade <- slope * 10
  
  return(slope_decade)
}


# ------------------------------------------------------------
# 4B. Function to calculate Mann-Kendall p-value
# ------------------------------------------------------------

mk_fun <- function(x) {
  
  # Remove missing values
  x <- x[!is.na(x)]
  
  # Require enough valid annual values
  if (length(x) < 10) {
    return(NA)
  }
  
  # Mann-Kendall test
  p_value <- mk.test(x)$p.value
  
  return(p_value)
}


# ------------------------------------------------------------
# 4C. Apply Sen's slope to Tmax
# ------------------------------------------------------------

tasmax_slope <- app(
  tasmax_annual,
  fun = sen_fun
)


# Mann-Kendall p-values for Tmax

tasmax_p <- app(
  tasmax_annual,
  fun = mk_fun
)


# ------------------------------------------------------------
# 4D. Apply Sen's slope to Tmin
# ------------------------------------------------------------

tasmin_slope <- app(
  tasmin_annual,
  fun = sen_fun
)


# Mann-Kendall p-values for Tmin

tasmin_p <- app(
  tasmin_annual,
  fun = mk_fun
)


# ------------------------------------------------------------
# 4E. Apply Sen's slope to precipitation
# ------------------------------------------------------------

pr_slope <- app(
  pr_annual,
  fun = sen_fun
)


# Mann-Kendall p-values for precipitation

pr_p <- app(
  pr_annual,
  fun = mk_fun
)


# ------------------------------------------------------------
# 4F. Classify statistical significance
# alpha = 0.05
# ------------------------------------------------------------

tasmax_sig <- tasmax_p < 0.05
tasmin_sig <- tasmin_p < 0.05
pr_sig     <- pr_p < 0.05


# ------------------------------------------------------------
# 4G. Inspect results
# ------------------------------------------------------------

tasmax_slope
tasmax_p
tasmax_sig

tasmin_slope
tasmin_p
tasmin_sig

pr_slope
pr_p
pr_sig


# ------------------------------------------------------------
# 4H. Check value ranges
# ------------------------------------------------------------

global(tasmax_slope, range, na.rm = TRUE)
global(tasmin_slope, range, na.rm = TRUE)
global(pr_slope, range, na.rm = TRUE)

global(tasmax_p, range, na.rm = TRUE)
global(tasmin_p, range, na.rm = TRUE)
global(pr_p, range, na.rm = TRUE)


# ------------------------------------------------------------
# 4I. Quick visual checks
# ------------------------------------------------------------

plot(
  tasmax_slope,
  main = "Sen's Slope - Tmax (°C/decade)"
)
lines(zambia)

plot(
  tasmin_slope,
  main = "Sen's Slope - Tmin (°C/decade)"
)
lines(zambia)

plot(
  pr_slope,
  main = "Sen's Slope - Precipitation (mm/decade)"
)
lines(zambia)

# ------------------------------------------------------------
# TASK 5: VISUALIZATION
# ------------------------------------------------------------

# ------------------------------------------------------------
# 5A. PREPARE DAILY MAXIMUM TEMPERATURE DATA FOR MAPPING
# ------------------------------------------------------------
tasmax_slope_map <- mask(tasmax_slope, zambia)

tasmax_df <- as.data.frame(
  tasmax_slope_map,
  xy = TRUE,
  na.rm = TRUE
)

names(tasmax_df)[3] <- "slope"

zambia_sf <- st_as_sf(zambia)

tasmax_sig_points <- as.points(
  tasmax_sig,
  values = TRUE,
  na.rm = TRUE
)

tasmax_sig_sf <- st_as_sf(tasmax_sig_points)


# ------------------------------------------------------------
# 5B. DAILY MAXIMUM TEMPERATURE TREND MAP
# ------------------------------------------------------------

tasmax_map <- ggplot() +
  
  geom_raster(
    data = tasmax_df,
    aes(x = x, y = y, fill = slope)
  ) +
  
  geom_sf(
    data = zambia_sf,
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_sf(
    data = tasmax_sig_sf,
    color = "black",
    shape = 16,
    size = 0.8,
  ) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "°C/decade"
  ) +
  
  coord_sf() +
  
  labs(
    title = "Spatial Trend in Annual Maximum Temperature",
    subtitle = "Zambia, 1981–2010",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      size = 10
    ),
    axis.text = element_text(
      size = 9
    ),
    legend.position = "bottom",
    panel.grid.major = element_line(
      linewidth = 0.2
    ),
    panel.grid.minor = element_blank()
  )


# Display map
tasmax_map


# Save map
ggsave(
  "tasmax_trend_map.png",
  plot = tasmax_map,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 5C. PREPARE DAILY MINIMUM TEMPERATURE DATA FOR MAPPING
# ------------------------------------------------------------
tasmin_slope_map <- mask(tasmin_slope, zambia)

tasmin_df <- as.data.frame(
  tasmin_slope_map,
  xy = TRUE,
  na.rm = TRUE
)

names(tasmin_df)[3] <- "slope"

tasmin_sig_points <- as.points(
  tasmin_sig,
  values = TRUE,
  na.rm = TRUE
)

tasmin_sig_sf <- st_as_sf(tasmin_sig_points)


# ------------------------------------------------------------
# 5D. DAILY MINIMUM TEMPERATURE TREND MAP
# ------------------------------------------------------------

tasmin_map <- ggplot() +
  
  geom_raster(
    data = tasmin_df,
    aes(x = x, y = y, fill = slope)
  ) +
  
  geom_sf(
    data = zambia_sf,
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_sf(
    data = tasmin_sig_sf,
    color = "black",
    shape = 16,
    size = 0.8,
  ) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "°C/decade"
  ) +
  
  coord_sf() +
  
  labs(
    title = "Spatial Trend in Annual Minimum Temperature",
    subtitle = "Zambia, 1981–2010",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      size = 10
    ),
    axis.text = element_text(
      size = 9
    ),
    legend.position = "bottom",
    panel.grid.major = element_line(
      linewidth = 0.2
    ),
    panel.grid.minor = element_blank()
  )

#Display map
tasmin_map

#Save map
ggsave(
  "tasmin_trend_map.png",
  plot = tasmin_map,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 5E. PREPARE PRECIPITATION DATA FOR MAPPING
# ------------------------------------------------------------
pr_slope_map <- mask(pr_slope, zambia)

pr_df <- as.data.frame(
  pr_slope_map,
  xy = TRUE,
  na.rm = TRUE
)

names(pr_df)[3] <- "slope"

pr_sig_points <- as.points(
  pr_sig,
  values = TRUE,
  na.rm = TRUE
)

pr_sig_sf <- st_as_sf(pr_sig_points)


# ------------------------------------------------------------
# 5F. PRECIPITATION TREND MAP
# ------------------------------------------------------------

pr_map <- ggplot() +
  
  geom_raster(
    data = pr_df,
    aes(x = x, y = y, fill = slope)
  ) +
  
  geom_sf(
    data = zambia_sf,
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_sf(
    data = pr_sig_sf,
    color = "black",
    shape = 16,
    size = 0.8,
  ) +
  
  scale_fill_gradient2(
    low = "brown",
    mid = "white",
    high = "darkgreen",
    midpoint = 0,
    name = "mm/decade"
  ) +
  
  coord_sf() +
  
  labs(
    title = "Spatial trend in Annual Precipitation",
    subtitle = "Zambia, 1981–2010",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.title = element_text(
      size = 10
    ),
    axis.text = element_text(
      size = 9
    ),
    legend.position = "bottom",
    panel.grid.major = element_line(
      linewidth = 0.2
    ),
    panel.grid.minor = element_blank()
  )

# Display Map
pr_map

# Save Map
ggsave(
  "precipitation_trend_map.png",
  plot = pr_map,
  width = 7,
  height = 6,
  dpi = 300
)

