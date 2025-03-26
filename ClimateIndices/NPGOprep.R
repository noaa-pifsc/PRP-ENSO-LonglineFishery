# The purpose of this script is to
# download the NPGO index data and
# package it neatly in a csv

# Set up the environment
library(tidyverse, quietly = TRUE)
library(lubridate)
library(here)

# Donwload data
# Skip the rows that explain the data, which frustratingly includes the column
# names.
NPGO_full <-read_table("https://www.o3d.org/npgo/data/NPGO.txt",
                       skip = 26, 
                       col_names = FALSE)

# Add the column names in, which requires visiting the above url to get them
# We're going to align them with previous work, too.
colnames(NPGO_full) <- c("YEAR", "MONTH", "NPGO")

# And reorder the columns to align with previous work
NPGO <- NPGO_full |>
  relocate(MONTH)

# Save
write_csv(NPGO, here("ClimateIndices", "NPGO.csv"))
write_delim(NPGO, here("ClimateIndices", "NPGO.dat"))