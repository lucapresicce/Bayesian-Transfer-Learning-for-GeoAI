# Sourcing code for Reproducibility 

## spBPS R package
Most of the Data Analysis, and Simulations, used the code implemented in the R package `spBPS`, available at this repository [spBPS Github Repository](https://github.com/lucapresicce/spBPS). This is an optimized Rcpp-based package, which provides the main functions to perform the Double Bayesian Predictive Stacking approach presented in the manuscript related to this repository. Since the package is not already available on CRAN (working for submission, and hopefully soon available), we suggest the following procedure before starting the execution of the Script to reproduce the results. Firstly, must be installed the `devtools` R library, then, check for its presence on your device, otherwise install it:
```{r, echo = F, eval = F, collapse = TRUE}
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```
Once `devtools` is available on the local machine, installation proceeds as follows:
```{r, echo = F, eval = F, collapse = TRUE}
devtools::install_github("lucapresicce/spBPS")
```

## Rcpp source file
However, for the Scripts `exec_comparison_seq_M.R`, and `exec_comparison_seq.R` to reproduce the results in Section 4.2, and Supplement Section 5.2 respectively, it is mandatory to compile a `.cpp` file, which contains internal functions of package 'spBPS'. Even if the following procedure can be found within the scripts themself, we report here the execution code lines (already present in the scripts):
```{r, echo = F, eval = F, collapse = TRUE}
# checking for Rcpp library
if (!require(Rcpp)) {
  install.packages("Rcpp", dependencies = TRUE)
}

# sourcing Rcpp functions
Sys.setenv(PKG_CXXFLAGS = "-Ofast")
sourceCpp("code/src/code.cpp")
```
There are two specifications to point out here. Formerly, it is mandatory to set as a working directory the main folder of this repository (following the instructions in the Workflow to reproduce the analyses). Lastly, the optimizing flag is called before compiling by the command `Sys.setenv(PKG_CXXFLAGS = "-Ofast")`. That is not mandatory, however, computational performances may vary due to its effects on optimization matrix calculations.






