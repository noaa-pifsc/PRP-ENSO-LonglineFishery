# The purpose of this script is to
# download the PDO index data and
# package it neatly in a csv

# These methods borrow directly from those used for the SAFE report:
# https://github.com/pwoodworth-jefcoats/climate-indicators/blob/main/04-Pacific-Decadal-Oscillation.qmd

# Set up the environment
library(tidyverse, quietly = TRUE)
library(lubridate)
library(here)

# Load data
PDO_full <- read_csv("https://psl.noaa.gov/pdo/data/pdo.timeseries.ersstv5.csv", 
                     show_col_types = FALSE)

# Rename the PDO column to something reasonable
PDO_full <- rename(PDO_full, PDO = `PDO from ERSST V5 https://psl.noaa.gov/pdo/ Using EOF from 1920 to 2014 for N Pacific (see webpage)`)

# Remove the months with a missing value (-9999)
miss_idx <- which(PDO_full$PDO == -9999)
PDO_full <- slice(PDO_full, -miss_idx)

# Separate out the years and months
PDO_full <- PDO_full |>
  mutate(Month = month(Date)) |>
  mutate(Year = year(Date))

PDO <- PDO_full |>
  select(Month, Year, PDO)

# Save everything
write_csv(PDO, here("ClimateIndices", "PDO.csv"))
write_delim(PDO, here("ClimateIndices", "PDO.dat"))
