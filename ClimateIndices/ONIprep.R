# The purpose of this script is to:
# download the ONI
# classify phases as El Niño, La Niña, Neutral, other
# package it neatly in a csv

# These methods borrow directly from those use for the SAFE report:
# https://github.com/pwoodworth-jefcoats/climate-indicators/blob/main/03-El-Nino-Southern-Oscillation.qmd

# Set up the environment
library(tidyverse, quietly = TRUE)
library(here)

# Load data
# While I very much would have liked to load the data directly from either
# https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt
# or
# https://psl.noaa.gov/data/correlation/oni.data
# both these site host data in a format that is not easily machine-readable.
# I spent over an hour on it and decided to move on.  
# Why is NOAA so BAD at providing accessible data?!
# If you have code to easily ingest data from either of the above sites,
# please do share it.
# Until then, I copied the data from the CPC site above, pasted it into Excel,
# and saved the file as ONI_full.csv.  (sigh.)
ONI_full <- read_csv(here("ClimateIndices", "ONI_full.csv"))

# Round data to the tenth of a degree
ONI_rounded <- round(ONI_full$ANOM, digits = 1)

# Identifying El Niños and La Niñas, based on the criterion that they persist for at least 5 consecutive seasons
ElNino_idx = which(ONI_rounded >= 0.5);
LaNina_idx = which(ONI_rounded <= -0.5);
Neutral_idx = which(ONI_rounded > -0.5 & ONI_rounded < 0.5);
ElNino <- array(0, dim = c(length(ONI_rounded), 3))
LaNina <- array(0, dim = c(length(ONI_rounded), 3))
Neutral <- array(0, dim = c(length(ONI_rounded), 3))

ElNino[ElNino_idx,1] = ONI_rounded[ElNino_idx]
LaNina[LaNina_idx,1] = ONI_rounded[LaNina_idx]
Neutral[,1] = NA # to avoid omitting months with ONI = 0
Neutral[Neutral_idx,1] = ONI_rounded[Neutral_idx]

if (ElNino[1,1] > 0) {
  ElNino[1,2] = 1
  ElNino[1,3] = 1
}

if (LaNina[1,1] < 0) {
  LaNina[1,2] = 1
  LaNina[1,3] = 1
}

if (!is.na(Neutral[1,1])) {
  Neutral[1,2] = 1
  Neutral[1,3] = 1
}

for (r in seq(2, length(ONI_rounded), 1)) {
  if (ElNino[r,1] > 0) {
    ElNino[r,2] = 1
    ElNino[r,3] = ElNino[r,2] + ElNino[r - 1,3]
  }
  
  if (LaNina[r,1] < 0) {
    LaNina[r,2] = 1
    LaNina[r,3] = LaNina[r,2] + LaNina[r - 1,3]
  }
  
  if (!is.na(Neutral[r,1])) {
    Neutral[r,2] = 1
    Neutral[r,3] = Neutral[r,2] + Neutral[r - 1,3]
  }
}


for (l in seq(4, 1, -1)) {
  for (r in seq(2, length(ONI_rounded), 1)) {
    if (ElNino[r,3] == 0 && ElNino[r - 1,3] <= l) {
      ElNino[r - 1,3] = 0
    }
    if (LaNina[r,3] == 0 && LaNina[r - 1,3] <= l) {
      LaNina[r - 1,3] = 0
    }
    if (Neutral[r,3] == 0 && Neutral[r - 1,3] <= l) {
      Neutral[r - 1,3] = 0
    }
  }
}

if (ElNino[dim(ONI_full)[1],2] == 1 & ElNino[dim(ONI_full)[1],3] == 1) {
  ElNino[dim(ONI_full)[1],3] = 0
}
if (LaNina[dim(ONI_full)[1],2] == 1 & LaNina[dim(ONI_full)[1],3] == 1) {
  LaNina[dim(ONI_full)[1],3] = 0
}
if (Neutral[dim(ONI_full)[1],2] == 1 & Neutral[dim(ONI_full)[1],3] == 1) {
  Neutral[dim(ONI_full)[1],3] = 0
}

pos_idx = which(ElNino[,3] == 0)
ElNino[pos_idx,1] = NA
neg_idx = which(LaNina[,3] == 0)
LaNina[neg_idx,1] = NA
neu_idx = which(Neutral[,3] == 0)
Neutral[neu_idx,1] = NA

# Add a column for CentralMonth (in keeping with prior format)
CentralMonth <- array(NA, dim = c(dim(ONI_full)[1], 1))
CentralMonth[which(ONI_full$SEAS == "DJF"), 1] <- 1
CentralMonth[which(ONI_full$SEAS == "JFM"), 1] <- 2
CentralMonth[which(ONI_full$SEAS == "FMA"), 1] <- 3
CentralMonth[which(ONI_full$SEAS == "MAM"), 1] <- 4
CentralMonth[which(ONI_full$SEAS == "AMJ"), 1] <- 5
CentralMonth[which(ONI_full$SEAS == "MJJ"), 1] <- 6
CentralMonth[which(ONI_full$SEAS == "JJA"), 1] <- 7
CentralMonth[which(ONI_full$SEAS == "JAS"), 1] <- 8
CentralMonth[which(ONI_full$SEAS == "ASO"), 1] <- 9
CentralMonth[which(ONI_full$SEAS == "SON"), 1] <- 10
CentralMonth[which(ONI_full$SEAS == "OND"), 1] <- 11
CentralMonth[which(ONI_full$SEAS == "NDJ"), 1] <- 12

# Put everything together
ONI_withPhases <- cbind(CentralMonth, ONI_full$YR, ONI_rounded, LaNina[,1], Neutral[,1], ElNino[,1])

# Rename the columns to match previous format
colnames(ONI_withPhases) <- c("CentralMonth", "YR", "ONI", "LaNina", "Neutral", "ElNino")

ONI_withPhases <- as_tibble(ONI_withPhases)
ONI <- ONI_withPhases |>
  select(CentralMonth, YR, ONI)

# Save everything
write_csv(ONI_withPhases, here("ClimateIndices", "ONI_withPhases.csv"))
write_csv(ONI, here("ClimateIndices", "ONI.csv"))

write_delim(ONI_withPhases, here("ClimateIndices", "ONI_withPhases.dat"))
write_delim(ONI, here("ClimateIndices", "ONI.dat"))