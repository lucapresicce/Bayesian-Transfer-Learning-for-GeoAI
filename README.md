# Bayesian-Transfer-Learning-and-Divide-Conquer-Models
This Repository contains the Reproducibility Material of "_Building Artificially Intelligent Geostatistical Systems Using Bayesian Predictive Stacking_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee). The following includes a roadmap for this repository, which follows the Workflow to reproduce the analyses. Comprehensive descriptions, and suggestions, for performing the analyses are provided subsequently.
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

This section provides an extensive Workflow to reproduce all the numbers, and figures displayed in "_Building Artificially Intelligent Geostatistical Systems Using Bayesian Predictive Stacking_" by [**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee. The Workflow is presented separately for each Section, and anticipated by a suggested setup to ease the execution of the analyses.

### Before starting

Since the structure of the R Scripts, the computations are organized considering the starting working directory of the entire repository. As a matter of fact, the scripts begin with:
```{r, echo = F, eval = F, collapse = TRUE}
setwd(".../Bayesian-Transfer-Learning-and-Divide-Conquer-Models")
```
where `".../"` represents the path on the author's laptop, and then, the directory path where this repository should be placed before executing the scripts. In particular, the best option is to clone this repository on the local machine, by executing the following block of code into a `shell`. Once the location to clone this repository, open the command line and execute:
```{sh}
git clone https://github.com/lucapresicce/Bayesian-Transfer-Learning-and-Divide-Conquer-Models.git
```
If not possible, it is possible to execute the scripts by omitting the `setwd("...")` command, but it is mandatory to create two folders in the working directory:
* _code_: in which the `src` folder (from the `code` folder of this repository) must be copied, allowing the compilation of the `.cpp` file needed;
* _output_: this allows you to save the results and figures directly inside it.

Lastly, check the installation of the packages used in the scripts. The most important is the 'spBPS' package, for which installation the `devtools` library is required:
```{r}
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```
Once devtools is available on the local machine, installation from Github repository proceeds as follows:
```{r}
devtools::install_github("lucapresicce/spBPS")
```
<!--

### Section 4.1

### Section 4.2

### Section 5.1

### Section 5.2

### Supplement Section 4.1

### Supplement Section 4.2


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


