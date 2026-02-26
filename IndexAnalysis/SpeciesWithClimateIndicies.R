# Diversity metrics and climate index correlations

# Correlations between climate indicies and % bigeye, yellowfin, swordfish, mahi, and pomfret, as well as CPUE for the same species
# Translating Phoebe's CatchDiversityWithClimateIndices.m matlab code for consistency.

#-------------------------------------------------------------------------------
# Load libraries and set up the environment
#-------------------------------------------------------------------------------
library(here)
library(dplyr)
library(readr)
library(terra)
library(tidyterra)
library(ggplot2)
library(rnaturalearth)

# set data directory
dataDir <- here('FisheryData', 'BRT data')

#-------------------------------------------------------------------------------
# Read in the data
#-------------------------------------------------------------------------------
# climate indices
setwd(dataDir)
oni <- read.csv('ONI_withPhases.csv') %>% 
  filter(YR >= 1995 & YR <= 2024) %>% 
  rename(Year=YR, Month=CentralMonth)
npgo <- read.csv('NPGO.csv') %>% 
  filter(YEAR >= 1995 & YEAR <= 2024) %>% 
  rename(Year=YEAR, Month=MONTH)
pdo <- read.csv('PDO.csv') %>% 
  filter(Year >= 1995 & Year <= 2024)
# join indices together
climIdx <- left_join(oni, pdo) %>% 
  left_join(npgo) %>% 
  mutate(MonthIndex=1:360)

# Fisheries data
# terra reads in the months as depth. Doesn't affect anything but just be aware
bet <- rast('TotalBigeyeCaught.nc')
swo <- rast('TotalSwordfishCaught.nc')
yft <- rast('TotalYellowfinCaught.nc')
pom <- rast('TotalPomfretCaught.nc')
mahi <- rast('TotalMahiCaught.nc')
pmus <- rast('TotalPMUS.nc')
effort <- rast('TotalEffort.nc')/1000
vessels <- rast('TotalVessels.nc')
fish <- list(bet, swo, yft, pom, mahi, pmus, effort)
names(fish) <- c('BET', 'SWO', 'YFT', 'POM', 'MAHI', 'PMUS', 'EFRT')

#-------------------------------------------------------------------------------
# Calculate fisheries metrics
#-------------------------------------------------------------------------------
# Calculate CPUE
cpue <- list()
for (i in 1:6) {
  cpue[[i]] <- fish[[i]]/fish$EFRT
  names(cpue)[i] <- names(fish)[i]
}

# Calcuclate % of catch
prop <- list()
for (i in 1:5) {
  prop[[i]] <- (fish[[i]]/fish$PMUS)*100
  names(prop)[i] <- names(fish)[i]
}

# Average over space each year-month
# CPUE
cpueTS <- list()
for (i in 1:6) {
  cpueTS[[i]] <- global(cpue[[i]], c('mean'), na.rm=T)
  cpueTS[[i]]$MonthIndex <- parse_number(rownames(cpueTS[[i]]))
  names(cpueTS)[i] <- names(cpue)[i]
}
# Percent catch by species
propTS <- list()
for (i in 1:5) {
  propTS[[i]] <- global(prop[[i]], c('mean'), na.rm=T)
  propTS[[i]]$MonthIndex <- parse_number(rownames(propTS[[i]]))
  names(propTS)[i] <- names(prop)[i]
}

# Climatologies (average over time)
# CPUE
cpueClim <- list()
for (i in 1:6) {
  cpueClim[[i]] <- mean(cpue[[i]], na.rm=T)
  names(cpueClim)[i] <- names(cpue)[i]
}
# Percent catch by species
propClim <- list()
for (i in 1:5) {
  propClim[[i]] <- mean(prop[[i]], na.rm=T)
  names(propClim)[i] <- names(prop)[i]
}
# And one for effort just in case
effortClim <- mean(effort, na.rm=T)
vesselClim <- sum(vessels, na.rm=T)
# Make a mask for confidentiality of maps
# any location with <3 vessels is set to NA, all others to 1
vesselClim[vesselClim < 3] <- NA
vesselClim[vesselClim >2] <- 1

#-------------------------------------------------------------------------------
# Correlations
#-------------------------------------------------------------------------------
# The R equivalent of matlabs 'corrcoef' is cor() but it doesn't give p-values, so I'm using cor.test() instead
# the p.value and estimate variables in the output holds the p-value and correlation respectively
# Set up empty dataframe
corVals <- data.frame(spp=NA, clim=NA, var=NA, r=NA, p=NA)
# Make vectors of species, climate variables 
species <- names(fish)[1:5]
climVars <- c('ONI', 'PDO', 'NPGO')
# make vectors longer since we are doing both CPUE and % of catch
evalVars <- c(rep('CPUE',(length(climVars)*length(species))), rep('Percent',,(length(climVars)*length(species))))
species <- c(species, species)
idx=1
# Loop through all species (twice), first for CPUE then % of catch
for (j in seq_along(species)) {
  # Select which metric dataset is used
  if (j <= 5) {
    dat <- cpueTS[[species[j]]] 
  } else { 
    dat <- propTS[[species[j]]]
  }
  # Loop through climate variables
  # Use Pearson correlation becauase that is what corrcoef in matlab uses to stay consistent
  for (i in seq_along(climVars)) {
    corTest <- cor.test(dat$mean, climIdx[,which(names(climIdx) == climVars[i])], method = 'pearson')
    corVals[idx,1] <- species[j]
    corVals[idx,2] <- climVars[i]
    corVals[idx,3] <- evalVars[idx]
    corVals[idx,4] <- corTest$estimate
    corVals[idx,5] <- corTest$p.value
    idx=idx+1
  }
}
# Make a smaller dataframe with just the significant (< 0.01) correlations
corValsSigOnly <- corVals[which(corVals$p < 0.01),]

#-------------------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------------------
setwd(here())
# Function to make a timeseries plot with CPUE and climate index
# Function is set up to do either CPUE or proportion of catch, not both in the same call. 
climVarTSplot <- function(TSData, climData, corData, fishUnit) {
  for (i in 1:nrow(corData)) {
    # Pull out the right data
    spp <- corData[i,'spp']
    clim <- corData[i,'clim']
    var <- corData[i,'var']
    posIdx <- climData[,clim]
    posIdx[which(climData[,clim] < 0)] <- 0
    negIdx <- climData[,clim]
    negIdx[which(climData[,clim] > 0)] <- 0
    xtext <- paste(climData$Year, climData$Month, sep='-')
    outName <- paste0('ENSO_sppWClimIndex_', spp, '_',var,'_',clim, '.png')
    yaxisLabel <- paste(spp, var, fishUnit)
    
    # Plot
    # Initiate saving file
    png(filename = outName, width = 9, height = 4, res = 150, units = 'in')
    # Set margins
    par(mar=c(5,4,4,5) + 0.1)
    # Clim index as bars in red and blue
    barplot(posIdx, ylim=c(-3,3), col='red',border='white', lwd=0.25)
    barplot(negIdx, ylim=c(-3,3), col='blue', border='white', lwd=0.25, add=T)
    text(x = 360, y=2.5, paste('r=',round(corData$r[i], 2)))
    text(x = 360, y=2, paste('p=',round(corData$p[i], 4)))
    mtext(clim, side=2, line=3)
    axis(2)
    # New plot on top of climate index
    par(new=T)
    # Time series of CPUE or percent of catch
    plot(TSData[[spp]]$mean, type='l', xlab='', ylab='', axes=F, lwd=1)
    mtext(yaxisLabel, side=4, line=3)
    axis(4)
    mtext('Date', side=1, line=3)
    # Add x-axis labels
    axis(1, at=seq(1,360,by=12), labels = xtext[seq(1,360,by=12)])
    box()
    # Save figure
    dev.off()
  }
}

# Proportion of catch
climVarTSplot(TSData = propTS, climData = climIdx, corData = corValsSigOnly[which(corValsSigOnly$var == 'Percent'),], fishUnit = 'of PMUS catch')
# CPUE
climVarTSplot(TSData = cpueTS, climData = climIdx, corData = corValsSigOnly[which(corValsSigOnly$var == 'CPUE'),], fishUnit = '(fish per 1000 hooks)')

#-------------------------------------------------------------------------------
# Make climatology plot
# Get land polygon. Can be done either with the 'borders' function during plotting, or beforehand with data from naturalearth
world <-ne_countries(country = 'united states of america', scale = 'medium', returnclass = 'sf') %>%
  select(1) #%>%
  #sf::st_transform(crs = crs(cpueClim))
worldcrop <- sf::st_crop(world, c(xmin=-180, xmax=-130, ymin=10,ymax=40))

# Loop through all species for both CPUE and % of catch
clims <- c(cpueClim, propClim)
for (i in 1:length(clims)) {
  # Apply confidentiality mask
  myDat <- clims[[i]]*vesselClim 
  sppNames <- names(clims)
  # name outfile
  outFile <- ifelse(i <= 6, paste0('ENSO_', sppNames[i], '_CPUE_meanClims.png'), paste0('ENSO_', sppNames[i], '_percCatch_meanClims.png'))
  print(outFile)
  # Make the plot
  p <- ggplot() +
    geom_spatraster(data=myDat) + 
    geom_sf(data=worldcrop, fill='gray30', color='white') +
    geom_segment(aes(x = c(-180), xend = c(-150), y = c(20), yend = 20)) +     # these are the lines dividing the area up into regions from Phoebe's paper
    geom_segment(aes(x = c(-180), xend = c(-150), y = c(26), yend = 26)) + 
    geom_vline(xintercept = -150) +    
    #borders('world', xlim=c(-180,-130), ylim=c(10,40), fill='darkgray', colour='white', linewidth=0.25) +
    coord_sf(expand = FALSE) +
    scale_fill_viridis_c(name = ifelse(i <= 6, paste(sppNames[i], 'CPUE'), paste(sppNames[i], '% of Catch')), na.value = 'transparent') +
    theme(panel.border = element_rect(color='black'),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "transparent"), 
          legend.key.height = unit(1, "null"),
          legend.key.width = unit(0.5, 'cm'),
          legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(size = 10)) +
    guides(fill=guide_colourbar(title.position='right', title.theme=element_text(angle=270, hjust=0.5, vjust=0.5)))
  ggsave(plot = p, filename = outFile, width = 8, height = 5, units = 'in')
}

       