# Date: November 17, 2025
# Author: Johanna Wren
# Email: johanna.wren@noaa.gov
# Description: Load and run diagnostics for BRTs on data for the Climate Variability project.  


#------------------------------------------------
# Load Libraries
library(here)
library(dplyr)
library(ggplot2)
library(tidync)
library(tidyverse)

#------------------------------------------------
# Clear workspace
rm(list=ls())

# Set working directory
mainDir <- here()

#------------------------------------------------
# Load Data
#------------------------------------------------
## Fishery data
setwd(file.path(mainDir, 'FisheryData/BRT data'))
# list all the nc data files
ncfiles <- list.files(pattern='nc')
# read in all files into a list
dataList <- lapply(ncfiles, function(x) tidync(x) %>% hyper_tibble() %>% mutate_if(is.character, as.numeric))
# Fix colnames in the oxygen data
colnames(dataList[[2]]) <- c('Depth at Oxygen 2mLperL', colnames(dataList[[1]])[2:4])
# change month indexing in oxygen data to be consistent with all other datasets
# start with month=1 for Jan 1995. O2 starts with month=1 Jan 1993
dataList[[2]]$Month <- dataList[[2]]$Month-24
# put all files into one dataframe
dataDF <- dataList %>% 
  reduce(full_join, by=c('Longitude', 'Latitude', 'Month'))
#--note-- Catchability has five missing values in month 81 (Sept 2001)

#------------------------------------------------
## Climate indices
# Function to clean up and make the climate indicies consistent
ClimVarClean <- function(ClimIndex, MonName, YrName) {
  idx <- ClimIndex %>% 
    rename(Mon=MonName, Yr=YrName) %>% 
    filter(Yr >= 1995, Yr <= 2024) %>%
    arrange(Yr, Mon) 
  # add Month index (because I sorted the dataframe I'm just going to hardcode 1:360 here)
  idx$Month <- 1:nrow(idx)
  return(idx)
}

# Read in raw data
npgoRaw <- read.csv('NPGO.csv')
head(npgoRaw)
meiRaw <- read.csv('MEIv2.csv')
head(meiRaw)
oniRaw <- read.csv('ONI_withPhases.csv')
head(oniRaw)
pdoRaw <- read.csv('PDO.csv')
head(pdoRaw)
# Apply function to climate index data
npgo <- ClimVarClean(npgoRaw, MonName = 'MONTH', YrName = 'YEAR')
mei <- ClimVarClean(meiRaw, MonName = 'SecondMonth', YrName = 'Year')
oni <- ClimVarClean(oniRaw, MonName = 'CentralMonth', YrName = 'YR')
pdo <- ClimVarClean(pdoRaw, MonName = 'Month', YrName = 'Year')
# Combine all climate index data into one dataframe (wide format)
climIdx <- list(pdo, npgo, oni, mei) %>% 
  reduce(full_join, by=c('Month', 'Yr', 'Mon'))

#------------------------------------------------
# Combine fisheries data with climate index data and rearrange for a better flow
allDat <- dataDF %>% 
  left_join(climIdx, by='Month') %>% 
  rename(MonthIndex=Month, Year=Yr, Month=Mon) %>% 
  arrange(MonthIndex) %>% 
  select(Year, Month, MonthIndex, Longitude, Latitude, everything())

#------------------------------------------------
# Add bait and leader material
# Bait
# bait is sanma/sardine until Dec 2020, then a transition period Jan 2021 - Dec 2022, and milkfish starting Jan 2023
allDat <- allDat %>% 
  mutate(Bait=if_else(Year <= 2020, 'Sanma/Sardine', if_else(Year <= 2022, 'Transition', 'Milkfish')))
# Leader
# leader material is poorly reported prior to Jul 2000 (month 307), is wire Aug 2000 (308) - May 2021 (317), and mono starting in Jun 2021 (318)
allDat <- allDat %>% 
  mutate(Leader=if_else(MonthIndex <= 307, 'PoorlyReported', if_else(MonthIndex >= 318, 'Mono', 'Wire')))
