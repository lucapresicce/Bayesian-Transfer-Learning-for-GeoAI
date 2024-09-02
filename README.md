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
where `".../"` represents the path on the user's machine, and then, the directory path where this repository should be placed before executing the scripts. The best option is to clone this repository on the local machine, by executing the following block of code into a `shell`. Once the location to clone this repository is chosen, open the command line and execute:
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

### Section 4.1

Running `exec_comparison_sim_M.R` produces the results, contained in the following objects: 
* _timing & RMSPE results_: `simulation_multivariate_5_500.Rdata`, `simulation_multivariate_5_1000.Rdata`, `simulation_multivariate_10_500.Rdata`, `simulation_multivariate_10_1000.Rdata`;
* _interpolation plots_: `surface_M_5_500.png`, `surface_M_5_1000.png`, `surface_M_10_500.png`, `surface_M_10_1000.png`;
* _uncertainty quantification plots_: `UC_M_5_500.png`, `UC_M_5_1000.png`, `UC_M_10_500.png`, `UC_M_10_1000.png`;
* _posterior credible interval plots_: `CIpost_M_5_500.png`, `CIpost_M_5_1000.png`, `CIpost_M_10_500.png`, `CIpost_M_10_1000.png`.

In this section are displayed `surface_M_5_500.png`, `UC_M_5_500.png`, `CIpost_M_5_500.png` as Figures, and the contents of `simulation_multivariate_5_500.Rdata`, `simulation_multivariate_5_1000.Rdata`, `simulation_multivariate_10_500.Rdata`, `simulation_multivariate_10_1000.Rdata` within a Table. The remaining contents not shown in this Section are presented in the Supplement Section 3.1.

Here the notation is the following: _type_setting_n_subsetsize_. For example, type = surface, setting = M (multivariate), n = 5 (thousand), and subset size = 500, lead to the surface plot interpolation of the $n=5000$ and $K=10$ dataset, for multivariate models, that is `surface_M_5_500.png`

### Section 4.2

Running `exec_comparison_seq_M.R` produces the results, contained in the following objects: 
* _interpolation plots_: `surface_M_TLwell.png`, `surface_M_TLmis.png`, `surface_M_TLhms.png`, `surface_M_TLSbps.png`;
* _posterior credible interval plots_: `CIpost_M_TLwell.png`, `CIpost_M_TLmis.png`, `CIpost_M_TLhms.png`, `CIpost_M_TLbps.png`.

This section displayed `surface_M_TLbps.png`, and `CIpost_M_TLbps.png` as Figures, and the RMSPEs are placed within the interpolation plots `surface_M_TLwell.png`, `surface_M_TLmis.png`, `surface_M_TLhms.png`, `surface_M_TLSbps.png`. The remaining contents not shown in this Section are presented in the Supplement Section 3.2.

Here the notation is the following: _type_setting_TLspecification_. For example, type = surface, setting = M (multivariate), specification = hms (highly misspecified), lead to the surface plot interpolation for the multivariate highly misspecified model, that is `surface_M_TLhms.png`.

### Section 5.1

Running `exec_analysis_univariate.R`, and `exec_analysis_univariate250.R`, produces the results, contained in the following objects: 
* _data analysis results_: `dataanalysis_univariate.Rdata`, `dataanalysis_univariate250.Rdata`;
* _interpolation & uncertainty quantification plots_: `dataanalysis_univariate.png`, `dataanalysis_univariate250.png`;
* _exploratory spatial data analysis_: `eda_univariate.png`.

In this section is displayed `dataanalysis_univariate.png` as a Figure, while the results in `dataanalysis_univariate.Rdata`, and `dataanalysis_univariate250.Rdata`, are described in the Section body along with a Table. While we present `eda_univariate.png` in the Supplement Section 6.

### Section 5.2

Running `exec_analysis_multivariate.R`, and `exec_analysis_multivariate250.R`, produces the results, contained in the following objects: 
* _data analysis results_: `dataanalysis_multivariate.Rdata`, `dataanalysis_multivariate250.Rdata`;
* _interpolation & uncertainty quantification plots_: `dataanalysis_multivariate.png`, `dataanalysis_multivariate250.png`.
* _exploratory spatial data analysis_: `eda_multivariate.png`.

In this section is displayed `dataanalysis_multivariate.png` as a Figure, while the results in `dataanalysis_multivariate.Rdata`, and `dataanalysis_multivariate250.Rdata`, are described in the Section body along with a Table. While we present `eda_multivariate.png` in the Supplement Section 6.

### Supplement Section 3.3

Running `exec_subset_sensitivity.R` produces the results, contained in the following object: 
* _subsets dimension sensitivity plot_: `subset_sens.png`.

In this section is displayed `subset_sens.png`.

### Supplement Section 5.1

Running `exec_comparison_sim.R` produces the results, contained in the following objects: 
* _timing & RMSPE results_: `simulation_univariate_5_500.Rdata`, `simulation_univariate_5_1000.Rdata`, `simulation_univariate_10_500.Rdata`, `simulation_univariate_10_1000.Rdata`;
* _interpolation plots_: `surface_5_500.png`, `surface_5_1000.png`, `surface_10_500.png`, `surface_10_1000.png`;
* _uncertainty quantification plots_: `UC_5_500.png`, `UC_5_1000.png`, `UC_10_500.png`, `UC_10_1000.png`;
* _posterior credible interval plots_: `CIpost_5_500.png`, `CIpost_5_1000.png`, `CIpost_10_500.png`, `CIpost_10_1000.png`.

In this section are displayed `surface_5_500.png`, `UC_5_500.png`, `CIpost_5_500.png`, `surface_5_1000.png`, `UC_5_1000.png`, `CIpost_5_1000.png`, `surface_10_500.png`, `UC_10_500.png`, `CIpost_10_500.png`, `surface_10_1000.png`, `UC_10_1000.png`, `CIpost_10_1000.png` as Figures, and the contents of `simulation_univariate_5_500.Rdata`, `simulation_univariate_5_1000.Rdata`, `simulation_univariate_10_500.Rdata`, `simulation_univariate_10_1000.Rdata`  within a Table.

Here the notation is the following: _type_setting_n_subsetsize_. For example, type = surface, setting = (univariate), n = 5 (thousand), and subset size = 500, lead to the surface plot interpolation of the $n=5000$ and $K=10$ dataset, for univariate models, that is `surface_5_500.png`

### Supplement Section 5.2

Running `exec_comparison_seq.R` produces the results, contained in the following objects: 
* _interpolation plots_: `surface_TLwell.png`, `surface_TLmis.png`, `surface_TLhms.png`, `surface_TLSbps.png`;
* _posterior credible interval plots_: `CIpost_TLwell.png`, `CIpost_TLmis.png`, `CIpost_TLhms.png`, `CIpost_TLbps.png`.

In this section are displayed `surface_TLwell.png, `CIpost_TLwell.png`, `surface_TLmis.png`, `CIpost_TLmis.png`, `surface_TLhms.png`, `CIpost_TLhms.png`, `surface_TLbps.png`, `CIpost_TLbps.png`, and the RMSPEs are placed within the interpolation plots `surface_TLwell.png`, `surface_TLmis.png`, `surface_TLhms.png`, `surface_TLSbps.png`.

Here the notation is the following: _type_setting_TLspecification_. For example, type = surface, setting = (univariate), specification = hms (highly misspecified), lead to the surface plot interpolation for the univariate highly misspecified model, that is `surface_TLhms.png`.


--------------------------------------------------------------------------------
## Contacts

| **Author**|**Maintainer** |**Reference** |
| :--- | :--- | :--- |
| Luca Presicce (l.presicce@campus.unimib.it), Sudipto Banerjee (sudipto@ucla.edu) | Luca Presicce (l.presicce@campus.unimib.it) | "_Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee)  |



