# Outputs

All the numbers, figures, and results of the experiment in "_Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach_", and its Supplementary appendix, are gathered here, as directly written by the execution of code in the `Script` folder. To this end, a short presentation of the `Output` folder contents follows.

## Result objects

The result objects here can be found in the reference work at:
* Section 4.1: `replications_results.RData`;
* Section 5: [`dataanalysis_multivariate.RData`](./dataanalysis_multivariate.RData), [`dataanalysis_multivariate250.RData`](./dataanalysis_multivariate250.RData), [`review_DataAppl_M.RData`](./review_DataAppl_M.RData);
* Supplement Section 2.3: [`upperbound.RData`](./upperbound.RData);
* Supplement Section 4.1: [`simulation_multivariate_5_500.RData`](./simulation_multivariate_5_500.RData), [`simulation_multivariate_5_1000.RData`](./simulation_multivariate_5_1000.RData), [`simulation_multivariate_10_500.RData`](./simulation_multivariate_10_500.RData), [`simulation_multivariate_10_1000.RData`](./simulation_multivariate_10_1000.RData), [`review_TimeComp_5k_M.RData`](./review_TimeComp_5k_M.RData), [`review_TimeComp_10k_M.RData`](./review_TimeComp_10k_M.RData);
* Supplement Section 6: [`dataanalysis_univariate.RData`](./dataanalysis_univariate.RData), [`dataanalysis_univariate250.RData`](./dataanalysis_univariate250.RData), [`review_DataAppl.RData`](./review_DataAppl.RData).

**Note:** The output file `replications_results.RData` is **not included in this repository** because its size exceeds GitHub's 100 MB limit (the file is approximately 210 MB). However, it is **fully reproducible** by running the script [`modifications_TL_M.R`](../script/modifications_TL_M.R).  
Please be aware that this script may take a **long time to execute**, depending on your system’s resources. If needed, the original `replications_results.RData` file can be provided upon request.

## Figures

The figures here can be found in the reference work at:
* Section 4.1: [`TL_pred.png`](./TL_pred.png), [`TL_post.png`](./TL_post.png);
* Section 4.2: [`heatmap-amortized.png`](./heatmap-amortized.png), [`parameters-amortized.png`](./parameters-amortized.png);
* Section 5: [`dataanalysis_multivariate_RR.png`](./dataanalysis_multivariate_RR.png), [`dataanalysis_multivariate_NDVI.png`](./dataanalysis_multivariate_NDVI.png), [`dataanalysis_multivariate_RR250.png`](./dataanalysis_multivariate_RR250.png), [`dataanalysis_multivariate_NDVI250.png`](./dataanalysis_multivariate_NDVI250.png), [`eda_multivariate.png`](./eda_multivariate.png);
* Supplement Section 2.3: [`upperbound_sim.png`](./upperbound_sim.png);
* Supplement Section 4.1: [`surface_M_5_500.png`](./surface_M_5_500.png), [`UC_M_5_500.png`](./UC_M_5_500.png), [`CIpost_M_5_500.png`](./CIpost_M_5_500.png), [`surface_M_5_1000.png`](./surface_M_5_1000.png), [`UC_M_5_1000.png`](./UC_M_5_1000.png), [`CIpost_M_5_1000.png`](./CIpost_M_5_1000.png), [`surface_M_10_500.png`](./surface_M_10_500.png), [`UC_M_10_500.png`](./UC_M_10_500.png), [`CIpost_M_10_500.png`](./CIpost_M_10_500.png), [`surface_M_10_1000.png`](./surface_M_10_1000.png), [`UC_M_10_1000.png`](./UC_M_10_1000.png), [`CIpost_M_10_1000.png`](./CIpost_M_10_1000.png);
* Supplement Section 4.2: [`subset_sens.png`](./subset_sens.png);
* Supplement Section 6: [`dataanalysis_univariate.png`](./dataanalysis_univariate.png), [`dataanalysis_univariate250.png`](./dataanalysis_univariate250.png), [`eda_univariate.png`](./eda_univariate.png).
