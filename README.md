# Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets
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

This section provides an extensive Workflow to reproduce all the numbers, and figures displayed in "_Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach_" by [**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee. The Workflow is presented separately for each Section, and anticipated by a suggested setup to ease the execution of the analyses.

### Working directory

Since the structure of the R Scripts, the computations are organized considering the starting working directory of the entire repository. As a matter of fact, the scripts begin with:
```{r, echo = F, eval = F, collapse = TRUE}
setwd(".../Bayesian-Transfer-Learning-and-Divide-Conquer-Models")
```
where `".../"` represents the path on the author's laptop, and then, the directory path where this repository should be placed before executing the scripts. In particular, the best option is to clone this repository on the local machine, by executing the following block of code into a `shell`. Once the location to clone this repository, open the command line and execute:
```{sh}
git clone https://github.com/lucapresicce/Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets.git
```
If not possible, it is possible to execute the scripts by omitting the `setwd("...")` command, but it is mandatory to create two folders in the working directory:
* _code_: in which the `src` folder (from the `code` folder of this repository) must be copied, allowing the compilation of the `.cpp` file needed;
* _output_: this allows you to save the results and figures directly inside it.

<!--
### Package environments

In order to facilitate the reproducibility of the analyses, the `renv` R library was used to provide a working environment. The `renv.lock` file can be found in the `code` folder of this repository, providing all the packages needed to perform all the simulation experiments, and real data analyses.

Lastly, check the installation of the packages used in the scripts (within the `renv.lock` environment). The most important is the 'spBPS' package, for which installation the `devtools` library is required:
```{r}
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```
Once devtools is available on the local machine, installation from Github repository proceeds as follows:
```{r}
devtools::install_github("lucapresicce/spBPS")
```
-->

### Section 4.1

Running `exec_comparison_sim_M.R` produces the results, contained in the following objects: 
* _timing & RMSPE results_: `simulation_multivariate_5_500.Rdata`, `simulation_multivariate_5_1000.Rdata`, `simulation_multivariate_10_500.Rdata`, `simulation_multivariate_10_1000.Rdata`;
* _interpolation plots_: `surface_M_5_500.png`, `surface_M_5_1000.png`, `surface_M_10_500.png`, `surface_M_10_1000.png`;
* _uncertainty quantification plots_: `UC_M_5_500.png`, `UC_M_5_1000.png`, `UC_M_10_500.png`, `UC_M_10_1000.png`;
* _posterior credible interval plots_: `CIpost_M_5_500.png`, `CIpost_M_5_1000.png`, `CIpost_M_10_500.png`, `CIpost_M_10_1000.png`.

In this section are displayed `surface_M_5_500.png` in Figure 1, `UC_M_5_500.png` in Figure 2, `CIpost_M_5_500.png` in Figure 3, and the contents of `simulation_multivariate_5_500.Rdata`, `simulation_multivariate_5_1000.Rdata`, `simulation_multivariate_10_500.Rdata`, `simulation_multivariate_10_1000.Rdata` in Table 1.

Here the notation is the following: _type_setting_n_subsetsize_. For example, type = surface, setting = M (multivariate), n = 5 (thousand), and subsetsize = 500, lead to the surface plot interpolation of the $n=5000$ and $K=10$ dataset, for multivariate models, that is `surface_M_5_500.png`

### Section 4.2

Running `exec_comparison_seq_M.R` produces the results, contained in the following objects: 
* _interpolation plots_: `surface_M_TLwell.png`, `surface_M_TLmis.png`, `surface_M_TLhms.png`, `surface_M_TLSbps.png`;
* _posterior credible interval plots_: `CIpost_M_TLwell.png`, `CIpost_M_TLmis.png`, `CIpost_M_TLhms.png`, `CIpost_M_TLbps.png`.

In this section are displayed `surface_M_TLhms.png` in Figure 4, `CIpost_M_TLhms.png` in Figure 5, and the RMSPEs in Table 2 are placed within the interpolation plots `surface_M_TLwell.png`, `surface_M_TLmis.png`, `surface_M_TLhms.png`, `surface_M_TLSbps.png`.

Here the notation is the following: _type_setting_TLspecification_. For example, type = surface, setting = M (multivariate), specification = hms (highly misspecified), lead to the surface plot interpolation for the multivariate highly misspecified model, that is `surface_M_TLhms.png`.

### Section 5.1

Running `exec_analysis_univariate.R` produces the results, contained in the following objects: 
* _data analysis results_: `dataanalysis_univariate.Rdata`;
* _interpolation & uncertainty quantification plots_: `dataanalysis_univariate.png`.

In this section are displayed `dataanalysis_univariate.png` in Figure 6, while the results in `dataanalysis_univariate.Rdata` are described in the Section body.


### Section 5.2

Running `exec_analysis_multivariate.R` produces the results, contained in the following objects: 
* _data analysis results_: `dataanalysis_multivariate.Rdata`;
* _interpolation & uncertainty quantification plots_: `dataanalysis_multivariate.png`.

In this section are displayed `dataanalysis_multivariate.png` in Figure 7, while the results in `dataanalysis_multivariate.Rdata` are described in the Section body.

### Supplement Section 4.1

Running `exec_comparison_sim.R` produces the results, contained in the following objects: 
* _timing & RMSPE results_: `simulation_univariate_5_500.Rdata`, `simulation_univariate_5_1000.Rdata`, `simulation_univariate_10_500.Rdata`, `simulation_univariate_10_1000.Rdata`;
* _interpolation plots_: `surface_5_500.png`, `surface_5_1000.png`, `surface_10_500.png`, `surface_10_1000.png`;
* _uncertainty quantification plots_: `UC_5_500.png`, `UC_5_1000.png`, `UC_10_500.png`, `UC_10_1000.png`;
* _posterior credible interval plots_: `CIpost_5_500.png`, `CIpost_5_1000.png`, `CIpost_10_500.png`, `CIpost_10_1000.png`.

In this section are displayed `surface_5_500.png` in Figure 2, `UC_5_500.png` in Figure 3, `CIpost_5_500.png` in Figure 4, `surface_5_1000.png` in Figure 5, `UC_5_1000.png` in Figure 6, `CIpost_5_1000.png` in Figure 7, `surface_10_500.png` in Figure 8, `UC_10_500.png` in Figure 9, `CIpost_10_500.png` in Figure 10, `surface_10_1000.png` in Figure 11, `UC_10_1000.png` in Figure 12, `CIpost_10_1000.png` in Figure 13, and the contents of `simulation_univariate_5_500.Rdata`, `simulation_univariate_5_1000.Rdata`, `simulation_univariate_10_500.Rdata`, `simulation_univariate_10_1000.Rdata` in Table 1.

Here the notation is the following: _type_setting_n_subsetsize_. For example, type = surface, setting = (univariate), n = 5 (thousand), and subsetsize = 500, lead to the surface plot interpolation of the $n=5000$ and $K=10$ dataset, for univariate models, that is `surface_5_500.png`

### Supplement Section 4.2

Running `exec_comparison_seq_M.R` produces the results, contained in the following objects: 
* _interpolation plots_: `surface_TLwell.png`, `surface_TLmis.png`, `surface_TLhms.png`, `surface_TLSbps.png`;
* _posterior credible interval plots_: `CIpost_TLwell.png`, `CIpost_TLmis.png`, `CIpost_TLhms.png`, `CIpost_TLbps.png`.

In this section are displayed `surface_TLwell.png` in Figure 14, `CIpost_TLwell.png` in Figure 15, `surface_TLmis.png` in Figure 16, `CIpost_TLmis.png` in Figure 17, `surface_TLhms.png` in Figure 18, `CIpost_TLhms.png` in Figure 19, `surface_TLbps.png` in Figure 20, `CIpost_TLbps.png` in Figure 21, and the RMSPEs in Table 2 are placed within the interpolation plots `surface_TLwell.png`, `surface_TLmis.png`, `surface_TLhms.png`, `surface_TLSbps.png`.

Here the notation is the following: _type_setting_TLspecification_. For example, type = surface, setting = (univariate), specification = hms (highly misspecified), lead to the surface plot interpolation for the univariate highly misspecified model, that is `surface_TLhms.png`.

<!--

describes the folder structure for all the code used in [**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee "_Building Artificially Intelligent Geostatistical Systems Using Bayesian Predictive Stacking_". For any further questions, or problems/bugs concerning reproducibility, please contact the author/maintainer [**Luca Presicce**](https://lucapresicce.github.io/) (l.presicce@campus.unimib.it).
He will be happy to help you! :)

---------------------------------------------------------------------------------------------------------------------------
* reproducibility-code
  * analyses
    * univariate
      * exec_analysis_univariate.R **[contains R script for univariate data analysis]**
      * univariate.rar **[compressed archive - contains data for univariate data analysis]**
    * multivariate
      * exec_analysis_multivariate.R **[contains R script for multivariate data analysis]**
      * multivariate.rar **[compressed archive - contains data for multivariate data analysis]**
  * simulations
    * univariate
      * exec_comparison_sim.R **[contains R script for simulation: _ASMK vs SMK_]**
      * exe_comparison_seq.R **[contains R script for simulation: _Transfer learning_]**
    * multivariate
      * exec_comparison_sim_M.R **[contains R script for simulation: _ASMK vs SMK_]**
      
-->

--------------------------------------------------------------------------------
## Contacts

| **Author**|**Maintainer** |**Reference** |
| :--- | :--- | :--- |
| Luca Presicce (l.presicce@campus.unimib.it), Sudipto Banerjee (sudipto@ucla.edu) | Luca Presicce (l.presicce@campus.unimib.it) | "_Building Artificially Intelligent Geostatistical Systems Using Bayesian Predictive Stacking_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee)  |

<!--
Maintainer: l.presicce@campus.unimib.it
Reference: **Luca Presicce** and Sudipto Banerjee (2024+) *"Accelerated Meta-Kriging for massive Spatial dataset"* 
-->


