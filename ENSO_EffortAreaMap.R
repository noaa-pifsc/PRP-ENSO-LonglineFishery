library(ggplot2)
library(dplyr)
library(terra)
library(tidyterra)
library(here)
library(sf)

# Read in fisheries data
myDat <- read.csv('~/Documents/ClimateVariability/ENSO_LonglineFishery/FisheryData/Longline data/DeepSets.csv')
myDat <- myDat[-which(myDat$BH_LON == 0),]  
myDat <- myDat %>% 
  mutate(BH_LON = (BH_LON*-1))
head(myDat)

# Make a blank spatraster that covers the domain of the data with a 1 degree resolution
r <- rast(ext(-180,ceiling(max(myDat$BH_LON)),floor(min(myDat$BH_LAT)),ceiling(max(myDat$BH_LAT))), res=1, vals=NA)

# Turn the data into a spatraster
# Serial number data
# Make a spatial vector of the fishery data
maskVect <- vect(myDat[,c(12,11,2)], geom=c('BH_LON', 'BH_LAT'))
# Turn that into a raster using r (above) by counting the number of unique permitnumbers in each cell
# Have to use the 'field' flag here otherwise it won't work
maskRast <- rasterize(maskVect, r, field='PERMITNUM', fun='n_distinct', background=NA)
# Set all counts < 3 to NA, and everything else to 1
maskRast[maskRast$PERMITNUM_n_distinct < 3] <- NA
maskRast[!is.na(maskRast$PERMITNUM_n_distinct)] <- 1
plot(maskRast)

# Effort data
effortVect <- vect(myDat[,c(12,11,3)], geom=c('BH_LON', 'BH_LAT'))
effortRast <- rasterize(effortVect, r, field='HOOKSSET', fun='sum', background=NA)
# multiply by the mask to remove confidential data
effortRast <- (effortRast*maskRast)/1000000
plot(effortRast)

# Save as ncdf
writeCDF(effortRast, 'FisheryData/BRT data/ENSO_EffortRaster.nc')

### Can't figure out how to make a plot look nice in ggplot with a projection other than mercator. So I'm just saving the effort raster and plotting it in pyGMT
# # # change projection on the raster
# effortMoll <- project(effortRast, '+proj=moll +lon_0=210')
# 
# # get extent from lon lat to m
# lon_min <- 160; lon_max <- 260
# lat_min <- -2; lat_max <- 45
# # Create an extent object and project it to Mollweide
# e <- ext(lon_min, lon_max, lat_min, lat_max)
# e_sfc <- as.polygons(e, crs="+proj=longlat +datum=WGS84")
# moll_proj <- "+proj=moll +lon_0=210"
# e_moll <- project(e_sfc, moll_proj)
# new_ext <- ext(e_moll)
# 
# world <-ne_countries(scale = 'medium', returnclass = 'sf', continent = 'north america') %>%
#   select(1)
# worldcrop <- #st_crop(world, e) %>%
#   st_transform(world, '+proj=moll +lon_0=210')
# 
# ggplot() +
#   geom_spatraster(data = effortMoll) +
#   geom_sf(data=worldcrop) + 
#   scale_fill_viridis_c(na.value = 'transparent') +
#   coord_sf(crs = "+proj=moll +lon_0=220", datum = sf::st_crs(4326), xlim=new_ext[1:2], ylim=new_ext[3:4], expand = F) +
#   theme_minimal() + 
#   theme(panel.border = element_rect(color='black'))
# 
# # 
# # make basemap
# world <-ne_countries(scale = 'medium', returnclass = 'sf') %>%
#   select(1) #%>%
# #sf::st_transform(crs = crs(cpueClim))
# worldcrop <- sf::st_crop(world, c(xmin=160, xmax=260, ymin=-2,ymax=50))
# worldMoll <- st_transform(worldcrop, '+proj=moll +lon_0=200')
# 
# ggplot() +
#   geom_spatraster(data=effortMoll) + 
#   geom_sf(data=worldMoll, fill='gray30', color='white') +
#   geom_vline(xintercept = 180) +
#   coord_sf(expand = F, crs = '+proj=moll +lon_0=200') +
#   scale_fill_viridis_c(na.value = 'transparent') + 
#   #guides(fill=guide_colourbar(title.position='right', title.theme=element_text(angle=270, hjust=0.5, vjust=0.5))) +
#   theme(panel.border = element_rect(color='black'),
#         #panel.background = element_rect(fill = "white"),
#         legend.background = element_rect(fill = "transparent"), 
#         legend.key.height = unit(1, "null"),
#         legend.key.width = unit(0.5, 'cm'),
#         legend.margin = margin(0, 0, 0, 0),
#         legend.text = element_text(size = 10))

# # plotting R different projections (https://www.r-bloggers.com/2019/04/zooming-in-on-maps-with-sf-and-ggplot2/)
# 
# 
# # Source - https://stackoverflow.com/a/74566254
# # Posted by Robert Hijmans, modified by community. See post 'Timeline' for change history
# # Retrieved 2026-02-27, License - CC BY-SA 4.0
# 
# library(ggplot2)
# p <- project(effortRast, method="near", "+proj=moll +lon_0=-140", mask=TRUE)
# land <- ne_countries(scale='medium', returnclass='sf')
# ggplot() +
#   geom_spatraster(data=p) +
#   geom_sf(data=land, fill='transparent') +
#   coord_sf(xlim = c(-180,-110), ylim=c(-5,45)) +
#   scale_fill_viridis_c(na.value='transparent') +
#   scale_x_continuous(breaks=seq(-180,-110,by=10))
# 
# 
# # Define a bounding box in WGS84 (lon/lat)
# bbox_area <- st_bbox(c(xmin = -20, xmax = 45, ymin = 30, ymax = 73), crs = st_crs(4326))
# 
# # Crop the original WGS84 data
# world_cropped <- st_crop(ne_countries(scale = "medium", returnclass = "sf"), bbox_area)
# 
# # Plot the cropped data and use coord_sf to apply the Mollweide projection
# ggplot(data = world_cropped) +
#   geom_sf() +
#   coord_sf(crs = "+proj=moll") # Projection applied here
# 
# ggplot() +
#   geom_spatraster(data=p) +
#   scale_fill_viridis_c(na.value = 'transparent')
