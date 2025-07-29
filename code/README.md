# Sourcing code for Reproducibility 

## spBPS R package
Most of the Data Analysis and Simulations used the code implemented in the R package `spBPS`, available on `CRAN `; otherwise, you can also find it at this repository [spBPS Github Repository](https://github.com/lucapresicce/spBPS). This is an optimized Rcpp-based package, which provides the main functions to perform the Double Bayesian Predictive Stacking approach presented in the manuscript related to this repository. We suggest installing the `spBPS` package directly from `CRAN` before starting the execution of the Script to reproduce the results. 
```{r, echo = F, eval = F, collapse = TRUE}
if (!require(spBPS)) {
  install.packages("spBPS", dependencies = TRUE)
}
```

Otherwise, you can follow the alternative procedure. Firstly, the `devtools` R library must be installed, then check for its presence on your device:
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
However, for some Scripts it is mandatory to compile a `.cpp` file, which contains internal functions of the package 'spBPS'. Even if the following procedure can be found within the scripts themself, we report here the execution code lines (already present in the scripts):
```{r, echo = F, eval = F, collapse = TRUE}
# checking for Rcpp library
if (!require(Rcpp)) {
  install.packages("Rcpp", dependencies = TRUE)
}

# sourcing Rcpp functions
Sys.setenv(PKG_CXXFLAGS = "-Ofast")
sourceCpp("code/src/code.cpp")
```
There are two specifications to point out here. Formerly, it was mandatory to set as a working directory the main folder of this repository (following the instructions in the Workflow to reproduce the analyses). Lastly, the optimizing flag is called before compiling by the command `Sys.setenv(PKG_CXXFLAGS = "-Ofast")`. That is not mandatory; however, computational performance may vary due to its effects on optimization matrix calculations.






