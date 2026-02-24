# Diversity metrics and climate index correlations

# Correlations between climate indicies and % bigeye, yellowfin, swordfish, mahi, and pomfret, as well as CPUE for the same species
# Translating Phoebe's matlab code for consistency. Admittedly this is a lot more difficult than I thought. I don't remember matlab much any more

#-------------------------------------------------------------------------------
# Load libraries and set up the environment
#-------------------------------------------------------------------------------
library(here)
library(dplyr)
library(readr)
library(terra)
library(tidyterra)
library(ggplot2)

myDir <- here()

#-------------------------------------------------------------------------------
# Read in the data
#-------------------------------------------------------------------------------
# climate indices
oni <- read.csv('FisheryData/BRT data/ONI_withPhases.csv') %>% 
  filter(YR >= 1995 & YR <= 2024) %>% 
  rename(Year=YR, Month=CentralMonth)
npgo <- read.csv('FisheryData/BRT data/NPGO.csv') %>% 
  filter(YEAR >= 1995 & YEAR <= 2024) %>% 
  rename(Year=YEAR, Month=MONTH)
pdo <- read.csv('FisheryData/BRT data/PDO.csv') %>% 
  filter(Year >= 1995 & Year <= 2024)
# join indices together
climIdx <- left_join(oni, pdo) %>% 
  left_join(npgo) %>% 
  mutate(MonthIndex=1:360)

# Fisheries data
# terra reads in the months as depth. Doesn't affect anything but just be aware
bet <- rast('FisheryData/BRT data/TotalBigeyeCaught.nc')
swo <- rast('FisheryData/BRT data/TotalSwordfishCaught.nc')
yft <- rast('FisheryData/BRT data/TotalYellowfinCaught.nc')
pom <- rast('FisheryData/BRT data/TotalPomfretCaught.nc')
mahi <- rast('FisheryData/BRT data/TotalMahiCaught.nc')
pmus <- rast('FisheryData/BRT data/TotalPMUS.nc')
effort <- rast('FisheryData/BRT data/TotalEffort.nc')/1000
vessels <- rast('FisheryData/BRT data/TotalVessels.nc')
fish <- list(bet, swo, yft, pom, mahi, pmus, effort)
names(fish) <- c('BET', 'SWO', 'YFT', 'POM', 'MAHI', 'PMUS', 'EFRT')

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

# Sum over space each year-month
# CPUE
cpueTS <- list()
for (i in 1:6) {
  cpueTS[[i]] <- global(cpue[[i]], c('mean', 'sum'), na.rm=T)
  cpueTS[[i]]$MonthIndex <- parse_number(rownames(cpueTS[[i]]))
  names(cpueTS)[i] <- names(cpue)[i]
}
# Proportion fish by species
propTS <- list()
for (i in 1:5) {
  propTS[[i]] <- global(prop[[i]], c('mean', 'sum'), na.rm=T)
  propTS[[i]]$MonthIndex <- parse_number(rownames(propTS[[i]]))
  names(propTS)[i] <- names(prop)[i]
}

# Climatologies
cpueClim <- list()
for (i in 1:6) {
  cpueClim[[i]] <- mean(cpue[[i]], na.rm=T)
  names(cpueClim)[i] <- names(cpue)[i]
}

propClim <- list()
for (i in 1:5) {
  propClim[[i]] <- mean(prop[[i]], na.rm=T)
  names(propClim)[i] <- names(prop)[i]
}
# And one for effort just in case
effortClim <- mean(effort, na.rm=T)
vesselClim <- sum(vessels, na.rm=T)
# Make a mask for confidentiality of maps
vesselClim[vesselClim < 3] <- NA
vesselClim[vesselClim >2] <- 1

#-------------------------------------------------------------------------------
# Correlations
#-------------------------------------------------------------------------------
# The R equivalent of matlabs 'corrcoef' is cor() but it doesn't give p-values, so I'm using cor.test() instead
# the p.value and estimate variables in the output holds the p-value and correlation respectively

corVals <- data.frame(spp=NA, clim=NA, var=NA, r=NA, p=NA)
species <- names(fish)[1:5]
climVars <- c('ONI', 'PDO', 'NPGO')
evalVars <- c(rep('CPUE',15), rep('Percent',15))
species <- c(species, species)
idx=1
for (j in seq_along(species)) {
  if (j <= 5) {
    dat <- cpueTS[[species[j]]] 
  } else { 
    dat <- propTS[[species[j]]]
  }
  for (i in seq_along(climVars)) {
    corTest <- cor.test(dat$sum, climIdx[,which(names(climIdx) == climVars[i])], method = 'pearson')
    corVals[idx,1] <- species[j]
    corVals[idx,2] <- climVars[i]
    corVals[idx,3] <- evalVars[idx]
    corVals[idx,4] <- corTest$estimate
    corVals[idx,5] <- corTest$p.value
    idx=idx+1
  }
}
corValsSigOnly <- corVals[which(corVals$p < 0.01),]

## plots
# Make a timeseries plot with CPUE and climate index
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
climVarTSplot(TSData = propTS, climData = climIdx, corData = corValsSigOnly[6:11,], fishUnit = 'of PMUS catch')
# CPUE
climVarTSplot(TSData = cpueTS, climData = climIdx, corData = corValsSigOnly[1:5,], fishUnit = '(fish per 1000 hooks)')

# Make climatology plot
world <- rnaturalearth::ne_countries(country = 'united states of america', scale = 'medium', returnclass = 'sf') %>%
  dplyr::select(1) %>%
  sf::st_transform(crs = crs(cpueClim))
worldcrop <- sf::st_crop(world, c(xmin=-180, xmax=-130, ymin=10,ymax=40))

clims <- c(cpueClim, propClim)
for (i in 1:length(clims)) {
  # Apply confidentiality mask
  myDat <- clims[[i]]*vesselClim 
  sppNames <- names(clims)
  outFile <- ifelse(i <= 6, paste0('ENSO_', sppNames[i], '_CPUE_meanClims.png'), paste0('ENSO_', sppNames[i], '_percCatch_meanClims.png'))
  print(outFile)
  # Make the plot
  p <- ggplot() +
    geom_spatraster(data=myDat) + 
    geom_sf(data=worldcrop, fill='gray30', color='white') +
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

       