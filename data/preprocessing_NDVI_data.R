###########################################################################################################################################
# PREPROCESSING FOR NDVI DATA
library(terra)
library(dplyr)
library(raster)

rm(list = ls())
# set the reproducibility-repository as working directory
setwd("D:/presi/Desktop/documenti/Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets/")

# load the downloaded data (from source available in README)
datadir <- file.path(getwd(), "data")
rawdata <- file.path(datadir, "MOD13C1.A2024113.061.2024134022931.hdf")
raster_data <- rast(rawdata[1])

# make it suitable for modeling
full_data <- raster::as.data.frame(raster_data, xy = TRUE)

# remove missing data
complete_indx <- complete.cases(full_data)

# wrangle outcomes the dataset of interest
NDVIdata <- full_data[complete_indx,] |> 
  dplyr::select(x, y, `"CMG 0.05 Deg 16 days NDVI"`, `"CMG 0.05 Deg 16 days red reflectance"`, `"CMG 0.05 Deg 16 days Avg sun zen angle"`) |> 
  transmute(lon = x,
            lat = y,
            lnNDVI = (log(1+((`"CMG 0.05 Deg 16 days NDVI"`*0.0001)+3000))), # +fill *scale factors (metadata)
            lnRedRefl = (log(1+((`"CMG 0.05 Deg 16 days red reflectance"`*0.0001)+3000))), # +fill *scale factors (metadata)
            lnZenith = (log(1+((`"CMG 0.05 Deg 16 days Avg sun zen angle"`*0.01)+10000))) ) # +fill *scale factors (metadata)

# reduce dimensions of data
set.seed(1997)
NDVIdata <- NDVIdata[sample(1:nrow(NDVIdata), 1.5e6, F),]

# save dataset for analysis in data foldere
save(NDVIdata, file = "data/NDVI_data_2024_05_12.RData")
