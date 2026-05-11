# PRP-ENSO-LonglineFishery
This repository includes the code necessary to replicate the results published
in *paper information here*.  There's also some code from approaches we 
abandoned as our study took shape.  

The scripts ending in `.R`, `.Rmd`, and `.qmd` run in [R](https://www.r-project.org/), 
though note that the `.qmd` scripts use [Quarto](https://quarto.org/) in 
[R](https://www.r-project.org/).  The scripts ending in `.m` run in [Matlab](https://www.mathworks.com/products/matlab.html).
The scripts ending `.ipynb` run in [python](https://www.python.org/).  The
scripts ending in `.jnl` run in [PyFerret](https://ferret.pmel.noaa.gov/Ferret/documentation/pyferret)
(see also their [repository](https://github.com/NOAA-PMEL/PyFerret)).

## Data Access and Preparation
The scripts below were used to access and prepare data for use in our study:  

* [ONIprep.R](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/ClimateIndices/ONIprep.R), 
[PDOprep.R](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/ClimateIndices/PDOprep.R), 
and [NPGOprep.R](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/ClimateIndices/NPGOprep.R) 
were used to access and format the climate indices we used  
* [GODASaccess.jnl](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/OceanData/GODASaccess.jnl) 
was used to access the ocean temperature data  
* [oxygen_depth_of_isopleth_13Mar25.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/OceanData/oxygen_depth_of_isopleth_13Mar25.m) 
and [find_isopleth_depth.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/OceanData/find_isopleth_depth.m) 
were used to prepare ocean oxygen data after they were manually downloaded as 
described in the paper  
* [OceanDataRegrid.jnl](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/OceanData/OceanDataRegrid.jnl) 
was used to grid ocean temperature and oxygen data to a common horizontal grid 
* [ObserverAccess.Rmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/ObserverAccess.Rmd) 
and [LogbookAccess.Rmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/LogbookAccess.Rmd) 
were used to access observer and logbook data, respectively
* [ObserverCombine.Rmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/ObserverCombine.Rmd) 
was used to combine multiple observer data frames into a single data frame.
Likewise [LogbookCombine.Rmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/LogbookCombine.Rmd) 
for logbook data  
* [LogbookGridding.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/LogbookGridding.m) 
was used to grid fishery-dependent data to the same grid as was used for ocean data.  

## Data Analysis, Figures, and Tables
The scripts below were used for data analysis resulting in the figures and tables
as described below:  

* __Figure 1:__ [ENSO_EffortAreaMap.ipynb](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/ENSO_EffortAreaMap.ipynb) 
and [RegionalEffortTrends.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/RegionalEffortTrends.m)
* __Figure 2:__ [RegionalEffortWithClimateIndices.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/RegionalEffortWithClimateIndices.m) 
and [EffortClimateHistograms.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/EffortClimateHistograms.m)
* __Figure 3:__ [CPUEtimeseries_figs.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/CPUEtimeseries_figs.m)  
* __Figure 4:__ [CompositeOceanography.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/CompositeOceanography.m) 
and [redblueTecplot.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/redblueTecplot.m)  
* __Figure 5:__ [OceanographyClimateCatch_CorrelationMaps.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/OceanographyClimateCatch_CorrelationMaps.m) 
and [redblueTecplot.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/redblueTecplot.m)  
* __Figure 6:__ [OceanographyClimateCatch_CorrelationMaps.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/OceanographyClimateCatch_CorrelationMaps.m) 
and [redblueTecplot.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/redblueTecplot.m)  
* __Figure S1:__ [RegionalEffortWithClimateIndices.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/RegionalEffortWithClimateIndices.m) 
and [CPUEtimeseries_figs.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/CPUEtimeseries_figs.m)  
* __Figure S2:__ [CompositeOceanography.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/CompositeOceanography.m) 
and [redblueTecplot.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/redblueTecplot.m)  
* __Figure S3:__ [OceanographyClimateCatch_CorrelationMaps.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/OceanographyClimateCatch_CorrelationMaps.m) 
and [redblueTecplot.m](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/IndexAnalysis/redblueTecplot.m)  
* __Table 1:__ [IndexTable.qmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/Writing/IndexTable.qmd)  
* __Table 2:__ [CPUEtable.qmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/Writing/CPUEtable.qmd)  

The scripts [ExploreBaitVariation.Rmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/ExploreBaitVariation.Rmd) 
and [ExploreLeaderVariation.Rmd](https://github.com/noaa-pifsc/PRP-ENSO-LonglineFishery/blob/main/FisheryData/ExploreLeaderVariation.Rmd)
informed our discussion of bait and leader changes over our period of interest.

## Remaining files
The remaining files in this repository didn't contribute directly to our final
paper.  We've left them here because some were built into workflows we used. 
Revising the code into a stand-alone repository proved error-prone as we
started down that path.  
