# Bayesian-Transfer-Learning-for-GeoAI
This Repository contains the Reproducibility Material of "_Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee). The following includes a roadmap for this repository, which follows the Workflow to reproduce the analyses. Comprehensive descriptions, and suggestions, for performing the analyses are provided subsequently.
In addition, also all the functions implemented in the package `spBPS` ([**Luca Presicce**](https://lucapresicce.github.io/)) are available in the **code** folder, for further details see [spBPS R package repository](https://github.com/lucapresicce/spBPS).

--------------------------------------------------------------------------------
## Roadmap of the Repository

| Folder | Contents |
| :--- | :---: |
| **code** | `spBPS` package info & Rcpp sourcing functions in `.cpp` format |
| **data** | preprocessed dataset in `.Rdata` format & preprocessing scripts |
| **output** | data analyses/simulations results in `.Rdata` format & figures in paper/supplement  |
| **script** | data analyses/simulations working scripts in `.R` format |


---------------------------------------------------------------------------------------------------------------------------
## Workflow for Reproducible Results

This section provides an extensive Workflow to reproduce all the numbers, and figures displayed in "_Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach_". The Workflow is presented separately for each Section, and anticipated by a suggested setup to ease the execution of the analyses.

### Working directory

Since the structure of the R Scripts, the computations are organized considering the starting working directory of the entire repository. The scripts begin with:
```{r, echo = F, eval = F, collapse = TRUE}
setwd(".../Bayesian-Transfer-Learning-for-GeoAI")
```
where `".../"` represents the path on the user's machine, and then, the directory path where this repository should be placed before executing the scripts. The best option is to clone this repository on the local machine, by executing the following block of code in a `shell`. Once the location to clone this repository is chosen, open the command line and execute:
```{sh}
git clone https://github.com/lucapresicce/Bayesian-Transfer-Learning-for-GeoAI.git
```
If not possible, it is possible to execute the scripts by omitting the `setwd("...")` command, but it is mandatory to create two folders in the working directory:
* _code_: in which the `src` folder (from the `code` folder of this repository) must be copied, allowing the compilation of the `.cpp` file needed;
* _output_: this allows you to save the results and figures directly inside it.

### Package environments

The most important is the 'spBPS' package, for which installation of the `devtools` library is required:
```{r}
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```
Once devtools is available on the local machine, installation from the Github repository proceeds as follows:
```{r}
devtools::install_github("lucapresicce/spBPS")
```

### Section 4.1 - Transfer Learning in $\mathscr{M}$-closed & $\mathscr{M}$-open settings

Running `modifications_TL_M.R` produces the results, contained in the following objects: 
* _replication results_: `replication_results.RData`;
* _posterior metrics plot_: `TL_post.png`;
* _predictive metrics plot_: `TL_pred.png`.

In this section are displayed `TL_post.png`, `TL_pred.png` as Figures, and the contents of 50 replications, collected in `replication_results.Rdata`.

### Section 4.2 - Amortized Bayesian Inference

Running `ABI_matrix.R` produces the results, contained in the following objects: 
* _interpolation plots_: `heatmap-amortized.png`;
* _posterior credible interval plots_: `parameters-amortized.png`.

This section displayed `heatmap-amortized.png` and `parameters-amortized.png` as Figures.

### Section 5 - Application

Running `exec_analysis_multivariate.R`, and `exec_analysis_multivariate250.R`, produces the results, contained in the following objects: 
* _data analysis results_: `dataanalysis_multivariate.Rdata`, `dataanalysis_multivariate250.Rdata`;
* _interpolation & uncertainty quantification plots_: `dataanalysis_multivariate_RR.png`, `dataanalysis_multivariate_NDVI.png`, `dataanalysis_multivariate_RR250.png`,`dataanalysis_multivariate_NDVI250.png`;
* _exploratory spatial data analysis_: `eda_multivariate.png`.

Running `modifications_DataAppl_M.R`, produces the results, contained in the following objects: 
* _AI model competitor results_: `review_DataAppl_M.RData`;

In this section are displayed `dataanalysis_multivariate_RR.png`, and `dataanalysis_multivariate_NDVI.png` as Figures, while the results in `dataanalysis_multivariate.Rdata`, `dataanalysis_multivariate250.Rdata`, and `review_DataAppl_M.RData`, are described in the Section body along with Tables. While we present `eda_multivariate.png` in the Supplement Section 7.

### Supplement Section 2.3 - Monte Carlo approximation for upper bound simulations

Running `Asymp_MC_sim.R` produces the results, contained in the following object: 
* _KL divergence upper bound simulations_: `upperbound.RData`;
* _KL divergence upper bound plots_: `upperbound_sim.png`.

In this section is displayed `upperbound_sim.png`.

### Supplement Section 4.1 - Computational Performance

Running `exec_comparison_sim_M.R` produces the results, contained in the following objects: 
* _timing & RMSPE results_: `simulation_multivariate_5_500.Rdata`, `simulation_multivariate_5_1000.Rdata`, `simulation_multivariate_10_500.Rdata`, `simulation_multivariate_10_1000.Rdata`;
* _interpolation plots_: `surface_M_5_500.png`, `surface_M_5_1000.png`, `surface_M_10_500.png`, `surface_M_10_1000.png`;
* _uncertainty quantification plots_: `UC_M_5_500.png`, `UC_M_5_1000.png`, `UC_M_10_500.png`, `UC_M_10_1000.png`;
* _posterior credible interval plots_: `CIpost_M_5_500.png`, `CIpost_M_5_1000.png`, `CIpost_M_10_500.png`, `CIpost_M_10_1000.png`.

Running `modifications_TimeComp.R`, produces the results, contained in the following objects: 
* _Running time competitor_: `review_TimeComp_5k_M.RData`, `review_TimeComp_10k_M.RData`;

In this section are displayed `surface_M_5_500.png`, `UC_M_5_500.png`, `CIpost_M_5_500.png` as Figures, and the contents of `simulation_multivariate_5_500.Rdata`, `simulation_multivariate_5_1000.Rdata`, `simulation_multivariate_10_500.Rdata`, `simulation_multivariate_10_1000.Rdata`, `review_TimeComp_5k_M.RData`, and `review_TimeComp_10k_M.RData` within a Table.

Here the notation is the following: _type_setting_n_subsetsize_. For example, type = surface, setting = M (multivariate), n = 5 (thousand), and subset size = 500, lead to the surface plot interpolation of the $n=5000$ and $K=10$ dataset, for multivariate models, that is `surface_M_5_500.png`

### Supplement Section 4.2 - Subset size sensitivity

Running `exec_subset_sensitivity.R` produces the results, contained in the following object: 
* _subsets dimension sensitivity plot_: `subset_sens.png`.

In this section is displayed `subset_sens.png`.

### Supplement Section 6 - Application (Univariate)

Running `exec_analysis_univariate.R`, and `exec_analysis_univariate250.R`, produces the results, contained in the following objects: 
* _data analysis results_: `dataanalysis_univariate.Rdata`, `dataanalysis_univariate250.Rdata`;
* _interpolation & uncertainty quantification plots_: `dataanalysis_univariate.png`, `dataanalysis_univariate250.png`;
* _exploratory spatial data analysis_: `eda_univariate.png`.

Running `modifications_DataAppl.R`, produces the results, contained in the following objects: 
* _AI model competitor results_: `review_DataAppl.RData`;

In this section is displayed `dataanalysis_univariate.png` as a Figure, while the results in `dataanalysis_univariate.Rdata`, `dataanalysis_univariate250.Rdata`, and `review_DataAppl.RData`, are described in the Section body along with Tables. While we present `eda_univariate.png` in the Supplement Section 7.


--------------------------------------------------------------------------------
## Contacts

| **Author**|**Maintainer** |**Reference** |
| :--- | :--- | :--- |
| Luca Presicce (l.presicce@campus.unimib.it), Sudipto Banerjee (sudipto@ucla.edu) | Luca Presicce (l.presicce@campus.unimib.it) | "_Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee)  |



 .
