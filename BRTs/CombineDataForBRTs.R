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
# Change colnames in the oxygen data to be consistent with the other files
colnames(dataList[[2]]) <- c('Oxy_2mL', colnames(dataList[[1]])[2:4])
# change month indexing in oxygen data to be consistent with all other datasets
# start with month=1 for Jan 1995. O2 starts with month=1 Jan 1993
dataList[[2]]$Month <- dataList[[2]]$Month-24
# put all files into one dataframe
dataDF <- dataList %>% 
  reduce(full_join, by=c('Longitude', 'Latitude', 'Month'))
# make the column names single words for ease later. Selecting the second word in the string
names(dataDF)[6:13] <- str_split(names(dataDF)[6:13], boundary('word'), simplify = T)[,2]
#--note-- Catchability has five missing values in month 81 (Sept 2001)

# Calculate CPUE
dataDF <- dataDF %>% 
  mutate(BET_CPUE=Bigeye/(Effort/1000), Mahi_CPUE=Mahi/(Effort/1000), Pom_CPUE=Pomfret/(Effort/1000), SWO_CPUE=Swordfish/(Effort/1000), YFT_CPUE=Yellowfin/(Effort/1000))
# Then we add some presence/absence, and log transform the CPUE to make it more normally distributed
dataDF <- dataDF %>% mutate(BET_PA = if_else(Bigeye>0, 1, 0, missing=NA), YFT_PA=if_else(Yellowfin >0, 1, 0, missing=NA), Mahi_PA=if_else(Mahi>0, 1, 0, missing=NA), Pomfre_PA=if_else(Pomfret>0, 1, 0, missing=NA), SWO_PA=if_else(Swordfish>0, 1, 0, missing=NA),
       Log_BET_CPUE=log(BET_CPUE+0.5), Log_YFT_CPUE=log(YFT_CPUE+0.5), Log_Pom_CPUE=log(Pom_CPUE+0.5), Log_Mahi_CPUE=log(Mahi_CPUE+0.5), Log_SWO_CPUE=log(SWO_CPUE+0.5), 
       Random=runif(nrow(dataDF),0,10))

#------------------------------------------------
## Climate indices
# Function to clean up and make the climate indices consistent
ClimVarClean <- function(ClimIndexFile, MonName, YrName) {
  idx <- read.csv(ClimIndexFile) %>% 
    rename(Mon=MonName, Yr=YrName) %>% 
    filter(Yr >= 1995, Yr <= 2024) %>%
    arrange(Yr, Mon) 
  # add Month index (because I sorted the dataframe I'm just going to hardcode 1:360 here)
  idx$Month <- 1:nrow(idx)
  return(idx)
}

# Apply function to climate index data
npgo <- ClimVarClean('NPGO.csv', MonName = 'MONTH', YrName = 'YEAR')
mei <- ClimVarClean('MEIv2.csv', MonName = 'SecondMonth', YrName = 'Year')
oni <- ClimVarClean('ONI_withPhases.csv', MonName = 'CentralMonth', YrName = 'YR')
pdo <- ClimVarClean('PDO.csv', MonName = 'Month', YrName = 'Year')
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
  mutate(Bait=if_else(Year <= 2020, 'SanmaSardine', if_else(Year <= 2022, 'Transition', 'Milkfish')))
# Leader
# leader material is poorly reported prior to Jul 2000 (month 67), is wire Aug 2000 (308) - May 2021 (317), and mono starting in Jun 2021 (318)
allDat <- allDat %>% 
  mutate(Leader=if_else(MonthIndex <= 67, 'PoorlyReported', if_else(MonthIndex >= 318, 'Mono', 'Wire'))) %>% 
  mutate_if(is.character, as.factor) %>% 
  data.frame()   # the BRT functions need a df not a tibble!

#------------------------------------------------
# Save file
saveRDS(allDat, 'ENSO_BRTdata.rds')


#------------------------------------------------
#------------------------------------------------
# Spot check some of the variables I calculated
# Gear and bait
ggplot(allDat, aes(MonthIndex, Bait)) + geom_path()
ggplot(allDat, aes(MonthIndex, Leader)) + geom_path()

# Presence/Absence
allDat %>% filter(MonthIndex == 360) %>% 
  ggplot() + 
    geom_raster(aes(Longitude, Latitude, fill=as.factor(Mahi_PA))) + 
    scale_fill_brewer(palette = 'BuPu', direction = 1)

# CPUE
allDat %>% filter(MonthIndex == 360) %>% 
  ggplot() + 
    geom_raster(aes(Longitude, Latitude, fill=BET_CPUE)) + 
    scale_fill_distiller(palette = 'BuPu', direction = 1)
