# Data Availability for Reproducibility 

"_Building Artificially Intelligent Geostatistical Systems Using Bayesian Predictive Stacking_" considers two sets of georeferenced real data, strictly related to critical problems in nowadays Climate Sciences. As a reminder, the analyses work on the Sea Surface Temperature data, and the Vegetation Index data, providing a full Bayesian inference over the entire global surface, including millions of units. Raw data is freely available for download from the following web portals: [NOAA Data Access](https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1sstd8dayR20190SQ_Lon0360.html), and [EarthData MOD13C1 v061](https://lpdaac.usgs.gov/products/mod13c1v061/), for whose metadata is available.

## Preprocessing

Since the original format and size of raw datasets, some preprocessing steps are mandatory before starting any analysis. Here two R scripts (`.R` format) are thus presented:
* `preprocessing_SST_data.R`: concerning the preprocessing of Sea Surface Temperature data, passing from raw data in `.nc` format (downloadable at the provided link) to the `.Rdata` object;
* `preprocessing_NDVI_data.R`: concerning the preprocessing of Vegetation Index data, passing from raw data in `.hdf` format (downloadable at the provided link) to the `.Rdata` object.

However, due to the massive dimensions, raw datasets are not loaded in this folder, but the preprocessed data are then available. In order to use this preprocessing script, first raw data must be downloaded from the link.

## Case Study Datasets

Introducing now `.Rdata` objects containing the preprocessed set of data used in Section 5. 
* `SST_data_2022_06_21.RData`: preprocessed data for univariate analysis in Section 5.1 (reference date of data 06/21/2022); 
* `NDVI_data_2024_05_12.RData`: preprocessed data for multivariate analysis in Section 5.2 (reference date of data 05/12/2024); 

The files here are analysis-ready, it suffices to load them (following the related scripts) to perform the analyses, as explained in the Workflow on main `README.md`.
