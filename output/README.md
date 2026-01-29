# Outputs

All the numbers, figures, and results of the experiment in "_Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach_", and its Supplementary appendix, are gathered here, as directly written by the execution of code in the `Script` folder. To this end, a short presentation of the `Output` folder contents follows.

## Result objects

The result objects here can be found in the reference work at:
* Section 5.1: `replicationsNNGPDBPS.RData`;
* Section 5.2: `replications_results.RData`;
* Section 6: [`dataanalysis_multivariate.RData`](./dataanalysis_multivariate.RData), [`dataanalysis_multivariate250.RData`](./dataanalysis_multivariate250.RData), [`review_DataAppl_M.RData`](./review_DataAppl_M.RData);
* Appendix Section C.1: [`simulation_multivariate_5_500.RData`](./simulation_multivariate_5_500.RData), [`simulation_multivariate_5_1000.RData`](./simulation_multivariate_5_1000.RData), [`simulation_multivariate_10_500.RData`](./simulation_multivariate_10_500.RData), [`simulation_multivariate_10_1000.RData`](./simulation_multivariate_10_1000.RData), [`review_TimeComp_5k_M.RData`](./review_TimeComp_5k_M.RData), [`review_TimeComp_10k_M.RData`](./review_TimeComp_10k_M.RData);
* Appendix Section C.2: [`upperbound.RData`](./upperbound.RData);

**Note:** The output files `replications_results.RData` and `replicationsNNGPDBPS.RData` are **not included in this repository** because their size exceeds GitHub's 100 MB limit. However, they are **fully reproducible** by running the scripts [`modifications_TL_M.R`](../script/modifications_TL_M.R), and [`replicationsNNGPDBPS.R`](../script/replicationsNNGPDBPS.R), respectively.
Please be aware that these scripts may take **long times to execute**, depending on your system’s resources. If needed, the original `replications_results.RData` or `replicationsNNGPDBPS.RData` file can be provided upon request.

## Figures

The figures here can be found in the reference work at:
* Section 5.2: [`TL_pred.png`](./TL_pred.png), [`TL_post.png`](./TL_post.png);
* Section 5.3: [`heatmap-amortized.png`](./heatmap-amortized.png), [`parameters-amortized.png`](./parameters-amortized.png);
* Section 6: [`dataanalysis_multivariate_RR.png`](./dataanalysis_multivariate_RR.png), [`dataanalysis_multivariate_NDVI.png`](./dataanalysis_multivariate_NDVI.png), [`dataanalysis_multivariate_RR250.png`](./dataanalysis_multivariate_RR250.png), [`dataanalysis_multivariate_NDVI250.png`](./dataanalysis_multivariate_NDVI250.png), [`eda_multivariate.png`](./eda_multivariate.png);
* Appendix Section C.1: [`surface_M_5_500.png`](./surface_M_5_500.png), [`UC_M_5_500.png`](./UC_M_5_500.png), [`CIpost_M_5_500.png`](./CIpost_M_5_500.png), [`surface_M_5_1000.png`](./surface_M_5_1000.png), [`UC_M_5_1000.png`](./UC_M_5_1000.png), [`CIpost_M_5_1000.png`](./CIpost_M_5_1000.png), [`surface_M_10_500.png`](./surface_M_10_500.png), [`UC_M_10_500.png`](./UC_M_10_500.png), [`CIpost_M_10_500.png`](./CIpost_M_10_500.png), [`surface_M_10_1000.png`](./surface_M_10_1000.png), [`UC_M_10_1000.png`](./UC_M_10_1000.png), [`CIpost_M_10_1000.png`](./CIpost_M_10_1000.png);
* Appendix Section C.2: [`upperbound_sim.png`](./upperbound_sim.png);
* Appendix Section C.3: [`subset_sens.png`](./subset_sens.png).
