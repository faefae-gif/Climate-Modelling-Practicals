# ============================================================
# PRACTICAL 4: ETCCDI EXTREME CLIMATE INDICES
# Mount Makulu, Zambia (1981–2010)
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# 0. LOAD PACKAGES
# ------------------------------------------------------------

library(climdex.pcic)
library(PCICt)
library(tidyverse)
library(trend)
library(patchwork)

# Install only once if needed:
# install.packages("patchwork")

# ------------------------------------------------------------
# 1. SET WORKING DIRECTORY
# ------------------------------------------------------------

setwd("C:/Users/User/Desktop/Climate_Practical_04")

getwd()
list.files()

# ------------------------------------------------------------
# 2. LOAD AND INSPECT DAILY OBSERVATIONS
# ------------------------------------------------------------

df <- read.csv("mount_makulu_obs_daily_1981_2010.csv")

names(df)
head(df)
str(df)
summary(df)

# ------------------------------------------------------------
# 3. DATA PREPARATION AND QUALITY CONTROL
# ------------------------------------------------------------

# Convert date from MM/DD/YYYY to Date
df$date <- as.Date(df$date, format = "%m/%d/%Y")

# Check date range
range(df$date)

# Check missing values
colSums(is.na(df))

# Check impossible temperature values
sum(df$tmax < df$tmin, na.rm = TRUE)

# Check negative precipitation
sum(df$precip < 0, na.rm = TRUE)

# Check duplicated dates
sum(duplicated(df$date))

# ------------------------------------------------------------
# 4. CREATE CLIMDEX INPUT OBJECT
# ------------------------------------------------------------

clim_dates <- as.PCICt(
  as.POSIXct(df$date, tz = "UTC"),
  cal = "gregorian"
)

ci <- climdexInput.raw(
  tmax = df$tmax,
  tmin = df$tmin,
  prec = df$precip,
  tmax.dates = clim_dates,
  tmin.dates = clim_dates,
  prec.dates = clim_dates,
  base.range = c(1981, 2010)
)

# ------------------------------------------------------------
# 5. PREPARE DATA FOR CLIMPACT
# ------------------------------------------------------------

climpact_data <- df %>%
  mutate(
    year  = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    day   = as.integer(format(date, "%d"))
  ) %>%
  select(year, month, day, precip, tmax, tmin)

head(climpact_data)

write.csv(
  climpact_data,
  "mount_makulu_climpact.csv",
  row.names = FALSE,
  na = "-99.9"
)

# ------------------------------------------------------------
# 6. IMPORT CLIMPACT ANNUAL INDEX OUTPUTS
# ------------------------------------------------------------

txx <- read.csv("mount_makulu_climpact_txx_ANN.csv")
su <- read.csv("mount_makulu_climpact_su_ANN.csv")
rx1day <- read.csv("mount_makulu_climpact_rx1day_ANN.csv")
r95ptot <- read.csv("mount_makulu_climpact_r95ptot_ANN.csv")

# ------------------------------------------------------------
# 7. CLEAN CLIMPACT OUTPUTS
# ------------------------------------------------------------

clean_climpact <- function(x) {
  
  value_col <- names(x)[2]
  
  x %>%
    filter(grepl("^[0-9]{4}$", Description..)) %>%
    transmute(
      year = as.integer(Description..),
      value = as.numeric(.data[[value_col]])
    ) %>%
    
    # Climpact missing-value codes
    mutate(
      value = ifelse(value <= -99, NA, value)
    )
}

txx_clean <- clean_climpact(txx)
su_clean <- clean_climpact(su)
rx1day_clean <- clean_climpact(rx1day)
r95ptot_clean <- clean_climpact(r95ptot)

# Inspect cleaned data
head(txx_clean)
head(su_clean)
head(rx1day_clean)
head(r95ptot_clean)

# Check number of annual observations
nrow(txx_clean)
nrow(su_clean)
nrow(rx1day_clean)
nrow(r95ptot_clean)

# Check realistic ranges after removing missing codes
range(txx_clean$value, na.rm = TRUE)
range(su_clean$value, na.rm = TRUE)
range(rx1day_clean$value, na.rm = TRUE)
range(r95ptot_clean$value, na.rm = TRUE)

# ------------------------------------------------------------
# 8. MANN-KENDALL AND SEN'S SLOPE
# ------------------------------------------------------------

analyse_trend <- function(data, index_name) {
  
  # Remove missing observations before testing
  x <- na.omit(data$value)
  
  mk <- mk.test(x)
  sen <- sens.slope(x)
  
  cat("\n====================================\n")
  cat("Index:", index_name, "\n")
  cat("Valid years:", length(x), "\n")
  cat("Mann-Kendall p-value:",
      round(mk$p.value, 4), "\n")
  cat("Sen's slope:",
      round(as.numeric(sen$estimates), 4),
      "per year\n")
  
  if (mk$p.value < 0.05) {
    cat("Result: Statistically significant trend (p < 0.05)\n")
  } else {
    cat("Result: No statistically significant trend (p >= 0.05)\n")
  }
}

analyse_trend(txx_clean, "TXx")
analyse_trend(su_clean, "SU")
analyse_trend(rx1day_clean, "RX1day")
analyse_trend(r95ptot_clean, "R95pTOT")

# ------------------------------------------------------------
# 9. PREPARE DATA FOR VISUALISATION
# ------------------------------------------------------------

thermal_df <- bind_rows(
  txx_clean %>% mutate(index = "TXx"),
  su_clean %>% mutate(index = "SU")
)

precip_df <- bind_rows(
  rx1day_clean %>% mutate(index = "RX1day"),
  r95ptot_clean %>% mutate(index = "R95pTOT")
)

# ------------------------------------------------------------
# 10. THERMAL EXTREMES PLOT
# ------------------------------------------------------------

p1 <- ggplot(
  thermal_df,
  aes(x = year, y = value, colour = index)
) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  labs(
    title = "Thermal Extreme Indices",
    x = "Year",
    y = "Index Value",
    colour = "Index"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ------------------------------------------------------------
# 11. HEAVY PRECIPITATION PLOT
# ------------------------------------------------------------

p2 <- ggplot(
  precip_df,
  aes(x = year, y = value, colour = index)
) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  labs(
    title = "Heavy Precipitation Indices",
    x = "Year",
    y = "Index Value",
    colour = "Index"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ------------------------------------------------------------
# 12. DUAL-PANEL FIGURE
# ------------------------------------------------------------

dual_panel <- p1 / p2 +
  plot_annotation(
    title = "ETCCDI Extreme Climate Indices at Mount Makulu (1981–2010)"
  )

print(dual_panel)

# ------------------------------------------------------------
# 13. SAVE FINAL FIGURE
# ------------------------------------------------------------

ggsave(
  filename = "Mount_Makulu_ETCCDI_dual_panel.png",
  plot = dual_panel,
  width = 10,
  height = 8,
  dpi = 300
)