#####################################################################################################################################################
## Multivariate DATA ANALYSIS - ####################################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-for-GeoAI")

# Packages --------------------------------------------------------------------
library(spBPS)
library(Rcpp)
library(RcppArmadillo)
library(mniw)
library(ggplot2)
library(tictoc)
library(parallel)
library(doParallel)
library(foreach)
library(MBA)
library(classInt)
library(RColorBrewer)
library(sp)
library(fields)
library(sf)
library(geoR)
library(spBayes)
library(bayesplot)
library(raster)
library(corrplot)
library(mapproj)
library(rworldmap)
library(ggmap)
library(maps)
library(gridExtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(scoringRules)
library(h2o)
library(dplyr)

# Data loading ----------------------------------------------------------------

# load preprocessed RData
load("data/NDVI_data_2024_05_12.RData")
full_data <- NDVIdata
rm(list = c("NDVIdata"))
names(full_data)

# take a look to data structure and responses variables
head(full_data)
dim(full_data)

# train data and test data
set.seed(1997)
test_ind <- sample.int(nrow(full_data), floor(1/6 * nrow(full_data)) )
train_data <- full_data[-test_ind, ]
test_data<- full_data[test_ind, ]

# select responses and predictors sets
response_set <- c("lnNDVI", "lnRedRefl")
q <- length(response_set)
predictor_set <- c("lnZenith")
p <- length(predictor_set)+1

# define train dimensions
N <- nrow(train_data)
crd_S <- matrix(as.matrix(train_data[, c("lon","lat")]), ncol = 2)
Y_S   <- matrix(as.matrix(train_data[, response_set]), ncol = q)
X_S   <- cbind(1, matrix(as.matrix(train_data[, predictor_set]), ncol = (p-1)))

# define test dimensions
U <- nrow(test_data)
crd_U <- matrix(as.matrix(test_data[, c("lon","lat")]), ncol = 2)
Y_U   <- matrix(as.matrix(test_data[, response_set]), ncol = q)
X_U   <- cbind(1, matrix(as.matrix(test_data[, predictor_set]), ncol = (p-1)))

# remove full dataset and free memory
rm(list = c("full_data", "train_data", "test_data"))
gc()


# EDA -------------------------------------------------------------------------

# subsample for feasible EDA
set.seed(1997)
eda_ind <- sample.int(N, round(N*0.025))

# computing the maximum distance
d.max <- sqrt((max(crd_S[,1]) - min(crd_S[,1]))^2 +
                (max(crd_S[,2]) - min(crd_S[,2]))^2)
d.max # around 386.0908 (multiply by 111.139) ~> 42,909 KM

# check the variogram for the first response
v.res_y1 <- variog(coords = crd_S[eda_ind, ], data = Y_S[eda_ind, 1],
                   uvec = (seq(0, 0.8*d.max, length = 50)))

vario.fit_y1 <- variofit(v.res_y1, cov.model="exponential")
summary(vario.fit_y1)
print(vario.fit_y1)

# check the variogram for the second response
v.res_y2 <- variog(coords = crd_S[eda_ind, ], data = Y_S[eda_ind, 2],
                   uvec = (seq(0, 0.8*d.max, length = 50)))

vario.fit_y2 <- variofit(v.res_y2, cov.model="exponential")
summary(vario.fit_y2)
print(vario.fit_y2)

# Response correlation matrix
cor(Y_S)

# free memory
gc()

# SubSubsample for model fitting ----------------------------------------------

set.seed(1997)
n <- 1000000
subsample <- sample(1:N, n)
crd_s <- crd_S[subsample, ]
y     <- Y_S[subsample, ]
x     <- X_S[subsample, ]

set.seed(1997)
u <- 250000
subsampleu <- sample(1:U, u)
crd_u <- crd_U[subsampleu, ]
y_u   <- Y_U[subsampleu, ]
x_u   <- X_U[subsampleu, ]

set.seed(1997)
folds <- split(subsampleu, rep(1:5, each = 50000))

# Scoring functions -------------------------------------------------------

to_sample_mat <- function(pred, j) {
  if (is.array(pred) && length(dim(pred)) == 3) {
    t(pred[, , j])
  } else {
    matrix(pred[, j], ncol = 1)
  }
}

compute_crps <- function(y_u, pred) {
  q      <- ncol(y_u)
  result <- sapply(1:q, function(j) mean(crps_sample(y_u[, j], to_sample_mat(pred, j))))
  setNames(result, paste0("CRPS_Y", 1:q))
}

compute_es <- function(y_u, pred) {
  n <- nrow(y_u)
  if (is.array(pred) && length(dim(pred)) == 3) {
    es_vals <- sapply(1:n, function(i) es_sample(y_u[i, ], t(pred[, i, ])))
  } else {
    es_vals <- sqrt(rowSums((y_u - pred)^2))
  }
  mean(es_vals)
}

# Fit linear model----------------------------------------------------

# dimension
n <- nrow(x)
p <- ncol(x)

# Prior parameters
mu_B <- matrix(0, p, q)
V_B <- diag(10, p)
nu <- 3
Psi <- diag(q)

# Run the Gibbs sampler
n_iter <- 1000
burn_in <- 500
tic()
set.seed(1234)
samples <- spBPS::bayesMvLMconjugate(y, x, mu_B, V_B, nu, Psi, n_iter, burn_in)
fit_conj_time <- toc()

# Extract posterior samples
B_samples <- samples$B_samples
Sigma_samples <- samples$Sigma_samples

# Summary statistics
B_mean <- apply(B_samples, 1:2, mean)
apply(B_samples, 1:2, quantile, c(0.025, 0.5, 0.975))
Sigma_mean <- apply(Sigma_samples, 1:2, mean)
apply(Sigma_samples, 1:2, quantile, c(0.025, 0.5, 0.975))

pred_bayesMvLMconjugate_fast <- function(X_new, B_samples, Sigma_samples,
                                         parallel = FALSE, n_cores = NULL,
                                         store_samples = TRUE) {
  n_iter <- dim(B_samples)[3]
  n_new  <- nrow(X_new)
  q      <- dim(Sigma_samples)[1]
  
  mem_gb <- n_new * q * n_iter * 8 / 1e9
  if (store_samples && mem_gb > 1) {
    message(sprintf(
      "Attenzione: Y_pred_samples occuperà ~%.1f GB. ",
      "Usa store_samples=FALSE per risparmiare memoria.", mem_gb
    ))
  }
  
  one_iter <- function(iter) {
    B     <- B_samples[, , iter]
    Sigma <- Sigma_samples[, , iter]
    
    Mu <- X_new %*% B                          
    L  <- chol(Sigma)                          
    Z  <- matrix(rnorm(n_new * q), n_new, q)
    Mu + Z %*% L                               
  }
  
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
    message(sprintf("Esecuzione parallela su %d core...", n_cores))
    
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl, c("X_new", "B_samples", "Sigma_samples",
                                    "n_new", "q"),
                              envir = environment())
      samples_list <- parallel::parLapply(cl, seq_len(n_iter), one_iter)
    } else {
      samples_list <- parallel::mclapply(seq_len(n_iter), one_iter,
                                         mc.cores = n_cores)
    }
  } else {
    pb <- txtProgressBar(0, n_iter, style = 3)
    samples_list <- vector("list", n_iter)
    for (iter in seq_len(n_iter)) {
      samples_list[[iter]] <- one_iter(iter)
      setTxtProgressBar(pb, iter)
    }
    close(pb)
  }
  
  Y_pred <- Reduce("+", samples_list) / n_iter
  
  if (store_samples) {
    Y_pred_samples <- array(unlist(samples_list),
                            dim = c(n_new, q, n_iter))
    list(Y_pred_samples = Y_pred_samples, Y_pred = Y_pred)
  } else {
    list(Y_pred_samples = NULL, Y_pred = Y_pred)
  }
}

# model aveluation 
results_blm <- data.frame(fold     = 1:5,
                          rmspe_Y1 = NA, rmspe_Y2 = NA,
                          crps_Y1  = NA, crps_Y2  = NA,
                          es       = NA)

for (k in 1:5) {
  
  x_u_k <- x_u[folds[[k]], ]
  y_u_k <- y_u[folds[[k]], ]
  
  pred_k         <- pred_bayesMvLMconjugate_fast(x_u_k, B_samples, Sigma_samples)
  pred_samples_k <- aperm(pred_k$Y_pred_samples, c(3, 1, 2))
  
  results_blm[k, "rmspe_Y1"] <- sqrt(mean((y_u_k[, 1] - pred_k$Y_pred[, 1])^2))
  results_blm[k, "rmspe_Y2"] <- sqrt(mean((y_u_k[, 2] - pred_k$Y_pred[, 2])^2))
  crps_k <- compute_crps(y_u_k, pred_samples_k)
  results_blm[k, "crps_Y1"]  <- crps_k[1]
  results_blm[k, "crps_Y2"]  <- crps_k[2]
  results_blm[k, "es"]       <- compute_es(y_u_k, pred_samples_k)
  
  cat("BLM fold", k, "done\n")
}

# Summary
apply(results_blm[, -1], 2, mean)
apply(results_blm[, -1], 2, sd)

# AI competitors ----------------------------------------------------------

# initialize h2o cluster
h2o.init()

# structure data for methods
X_df <- as.data.frame(x)
Y_df <- as.data.frame(y)
names(Y_df) <- paste0("Y", 1:ncol(Y_df))
data_all <- cbind(X_df, Y_df)
data_h2o <- as.h2o(data_all)

# names for variables
predictors <- colnames(X_df)
response_vars <- colnames(Y_df)

# models container
models_list <- list()
timing_log <- data.frame(
  response = character(),
  model = character(),
  time_sec = numeric(),
  stringsAsFactors = FALSE
)

# training and timing (for each response)
for (response in response_vars) {
  
  models_list[[response]] <- list()
  
  # AutoML
  t_start <- Sys.time()
  aml <- h2o.automl(
    x = predictors,
    y = response,
    training_frame = data_h2o,
    max_models = 10,
    seed = 1234,
    nfolds = 10
  )
  t_end <- Sys.time()
  models_list[[response]]$automl <- aml
  timing_log <- rbind(timing_log, data.frame(
    response = response, model = "automl",
    time_sec = as.numeric(difftime(t_end, t_start, units = "secs"))
  ))
  
  # Deep Learning
  t_start <- Sys.time()
  dl <- h2o.deeplearning(
    x = predictors,
    y = response,
    training_frame = data_h2o,
    seed = 1234,
    nfolds = 10
  )
  t_end <- Sys.time()
  models_list[[response]]$deeplearning <- dl
  timing_log <- rbind(timing_log, data.frame(
    response = response, model = "deeplearning",
    time_sec = as.numeric(difftime(t_end, t_start, units = "secs"))
  ))
  
  # Distributed Random Forest
  t_start <- Sys.time()
  drf <- h2o.randomForest(
    x = predictors,
    y = response,
    training_frame = data_h2o,
    seed = 1234,
    nfolds = 10
  )
  t_end <- Sys.time()
  models_list[[response]]$drf <- drf
  timing_log <- rbind(timing_log, data.frame(
    response = response, model = "drf",
    time_sec = as.numeric(difftime(t_end, t_start, units = "secs"))
  ))
  
  # GB
  t_start <- Sys.time()
  xgb <- h2o.gbm(
    x = predictors,
    y = response,
    training_frame = data_h2o,
    seed = 1234,
    nfolds = 10
  )
  t_end <- Sys.time()
  models_list[[response]]$xgboost <- xgb
  timing_log <- rbind(timing_log, data.frame(
    response = response, model = "xgboost",
    time_sec = as.numeric(difftime(t_end, t_start, units = "secs"))
  ))
  
}

# see running times (divided for outcomes)
print(timing_log)

# folds predictions
results_aml <- data.frame(fold=1:5, rmspe_Y1=NA, rmspe_Y2=NA, crps_Y1=NA, crps_Y2=NA, es=NA)
results_dl  <- data.frame(fold=1:5, rmspe_Y1=NA, rmspe_Y2=NA, crps_Y1=NA, crps_Y2=NA, es=NA)
results_drf <- data.frame(fold=1:5, rmspe_Y1=NA, rmspe_Y2=NA, crps_Y1=NA, crps_Y2=NA, es=NA)
results_gbm <- data.frame(fold=1:5, rmspe_Y1=NA, rmspe_Y2=NA, crps_Y1=NA, crps_Y2=NA, es=NA)

for (k in 1:5) {
  
  x_u_k    <- X_U[folds[[k]], ]
  y_u_k    <- Y_U[folds[[k]], ]
  newdat_k <- as.h2o(as.data.frame(x_u_k))
  
  aml_pred <- cbind(as.vector(h2o.predict(models_list[["Y1"]]$automl@leader, newdat_k)),
                    as.vector(h2o.predict(models_list[["Y2"]]$automl@leader, newdat_k)))
  
  dl_pred  <- cbind(as.vector(h2o.predict(models_list[["Y1"]]$deeplearning, newdat_k)),
                    as.vector(h2o.predict(models_list[["Y2"]]$deeplearning,  newdat_k)))
  
  drf_pred <- cbind(as.vector(h2o.predict(models_list[["Y1"]]$drf, newdat_k)),
                    as.vector(h2o.predict(models_list[["Y2"]]$drf, newdat_k)))
  
  gbm_pred <- cbind(as.vector(h2o.predict(models_list[["Y1"]]$xgboost, newdat_k)),
                    as.vector(h2o.predict(models_list[["Y2"]]$xgboost, newdat_k)))
  
  # AutoML
  crps_k <- compute_crps(y_u_k, aml_pred)
  results_aml[k, "rmspe_Y1"] <- sqrt(mean((y_u_k[, 1] - aml_pred[, 1])^2))
  results_aml[k, "rmspe_Y2"] <- sqrt(mean((y_u_k[, 2] - aml_pred[, 2])^2))
  results_aml[k, "crps_Y1"]  <- crps_k[1]
  results_aml[k, "crps_Y2"]  <- crps_k[2]
  results_aml[k, "es"]       <- compute_es(y_u_k, aml_pred)
  
  # Deep Learning
  crps_k <- compute_crps(y_u_k, dl_pred)
  results_dl[k, "rmspe_Y1"] <- sqrt(mean((y_u_k[, 1] - dl_pred[, 1])^2))
  results_dl[k, "rmspe_Y2"] <- sqrt(mean((y_u_k[, 2] - dl_pred[, 2])^2))
  results_dl[k, "crps_Y1"]  <- crps_k[1]
  results_dl[k, "crps_Y2"]  <- crps_k[2]
  results_dl[k, "es"]       <- compute_es(y_u_k, dl_pred)
  
  # DRF
  crps_k <- compute_crps(y_u_k, drf_pred)
  results_drf[k, "rmspe_Y1"] <- sqrt(mean((y_u_k[, 1] - drf_pred[, 1])^2))
  results_drf[k, "rmspe_Y2"] <- sqrt(mean((y_u_k[, 2] - drf_pred[, 2])^2))
  results_drf[k, "crps_Y1"]  <- crps_k[1]
  results_drf[k, "crps_Y2"]  <- crps_k[2]
  results_drf[k, "es"]       <- compute_es(y_u_k, drf_pred)
  
  # GBM
  crps_k <- compute_crps(y_u_k, gbm_pred)
  results_gbm[k, "rmspe_Y1"] <- sqrt(mean((y_u_k[, 1] - gbm_pred[, 1])^2))
  results_gbm[k, "rmspe_Y2"] <- sqrt(mean((y_u_k[, 2] - gbm_pred[, 2])^2))
  results_gbm[k, "crps_Y1"]  <- crps_k[1]
  results_gbm[k, "crps_Y2"]  <- crps_k[2]
  results_gbm[k, "es"]       <- compute_es(y_u_k, gbm_pred)
  
  cat("h2o fold", k, "done\n")
}

# AutoML
apply(results_aml[, -1], 2, mean)
apply(results_aml[, -1], 2, sd)

# Deep Learning
apply(results_dl[, -1], 2, mean)
apply(results_dl[, -1], 2, sd)

# DRF
apply(results_drf[, -1], 2, mean)
apply(results_drf[, -1], 2, sd)

# GBM
apply(results_gbm[, -1], 2, mean)
apply(results_gbm[, -1], 2, sd)


# DBPS -------------------------------------------------------------------

# hyperparameters values by looking at variograms
variofitphi1 <- 1 / vario.fit_y1$cov.pars[2]; variofitphi1
variofitalpha1 <- vario.fit_y1$cov.pars[1] / (vario.fit_y1$nugget+vario.fit_y1$cov.pars[1]); variofitalpha1
variofitphi2 <- 1 / vario.fit_y2$cov.pars[2]; variofitphi2
variofitalpha2 <- vario.fit_y2$cov.pars[1] / (vario.fit_y2$nugget+vario.fit_y2$cov.pars[1]); variofitalpha2

# hyperparameters values
(alfa_seq <- sort(c(variofitalpha1, variofitalpha2)))
(phi_seq <- sort(2*c(variofitphi1, variofitphi2)))
hyperpar <- list(alpha = alfa_seq, phi = phi_seq)

# number of available computational cores
n.core <- parallel::detectCores(logical = F)-1

# fitting DBPS with subset_size = 500
tic("spBPS")
out500 <- spBPS(data = list(Y = y, X = x),
                 priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                               V_r = diag(10, p),
                               Psi = diag(1, q),
                               nu = 3),
                 coords = crd_s,
                 hyperpar = hyperpar,
                 subset_size = 500,
                 cv_folds = 5,
                 combine_method = "pseudoBMA",
                 draws = 200,
                 n_cores = n.core)
new_time <- toc()

# fitting DBPS with subset_size = 250
tic("spBPS")
out250 <- spBPS(data = list(Y = y, X = x),
                 priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                               V_r = diag(10, p),
                               Psi = diag(1, q),
                               nu = 3),
                 coords = crd_s,
                 hyperpar = hyperpar,
                 subset_size = 250,
                 cv_folds = 5,
                 combine_method = "pseudoBMA",
                 draws = 200,
                 n_cores = n.core)
new_time <- toc()

# folds prediction
results_dbps500 <- data.frame(fold = 1:5, rmspe_Y1=NA, rmspe_Y2=NA, crps_Y1=NA, crps_Y2=NA, es=NA)
results_dbps250 <- data.frame(fold = 1:5, rmspe_Y1=NA, rmspe_Y2=NA, crps_Y1=NA, crps_Y2=NA, es=NA)

for (k in 1:5) {
  
  crd_u_k <- crd_u[folds[[k]], ]
  x_u_k   <- x_u[folds[[k]], ]
  y_u_k   <- y_u[folds[[k]], ]
  
  pred_k <- predict(out500,
                    newdata        = list(X = x_u_k, coords = crd_u_k),
                    draws          = 200L,
                    n_cores        = n.core,
                    pred_batch_size = 25)
                          
  
  pred_samples_k <- aperm(pred_k$Yu, c(3, 1, 2))
  pred_mean_k    <- apply(pred_k$Yu, c(1, 2), mean)
  
  crps_k <- compute_crps(y_u_k, pred_samples_k)
  results_dbps500[k, "rmspe_Y1"] <- sqrt(mean((y_u_k[, 1] - pred_mean_k[, 1])^2))
  results_dbps500[k, "rmspe_Y2"] <- sqrt(mean((y_u_k[, 2] - pred_mean_k[, 2])^2))
  results_dbps500[k, "crps_Y1"]  <- crps_k[1]
  results_dbps500[k, "crps_Y2"]  <- crps_k[2]
  results_dbps500[k, "es"]       <- compute_es(y_u_k, pred_samples_k)
  
  pred_k <- predict(out250,
                    newdata        = list(X = x_u_k, coords = crd_u_k),
                    draws          = 200L,
                    n_cores        = n.core,
                    pred_batch_size = 25)
  
  pred_samples_k <- aperm(pred_k$Yu, c(3, 1, 2))  
  pred_mean_k    <- pred_mean_k    <- apply(pred_k$Yu, c(1, 2), mean)
  
  crps_k <- compute_crps(y_u_k, pred_samples_k)
  results_dbps250[k, "rmspe_Y1"] <- sqrt(mean((y_u_k[, 1] - pred_mean_k[, 1])^2))
  results_dbps250[k, "rmspe_Y2"] <- sqrt(mean((y_u_k[, 2] - pred_mean_k[, 2])^2))
  results_dbps250[k, "crps_Y1"]  <- crps_k[1]
  results_dbps250[k, "crps_Y2"]  <- crps_k[2]
  results_dbps250[k, "es"]       <- compute_es(y_u_k, pred_samples_k)
  
  cat("DBPS fold", k, "done\n")
}

apply(results_dbps500[, -1], 2, mean)
apply(results_dbps500[, -1], 2, sd)

apply(results_dbps250[, -1], 2, mean)
apply(results_dbps250[, -1], 2, sd)


# Results collection ----------------------------------------------------------

# prediction on first folds for plotting
crd_u_1 <- crd_u[folds[[1]], ]
x_u_1   <- x_u[folds[[1]], ]
y_u_1   <- y_u[folds[[1]], ]
pred_out <- predict(out500,
                    newdata        = list(X = x_u_1, coords = crd_u_1),
                    draws          = 200L,
                    n_cores        = n.core,
                    pred_batch_size = 25)

pred_mat_W <- pred_out$Wu
pred_mat_Y <- pred_out$Yu

# statistics computations W
post_mean_W <- apply(pred_mat_W, c(1,2), mean)
post_var_W <- apply(pred_mat_W, c(1,2), sd)
post_qnt_W <- apply(pred_mat_W, c(1,2), quantile, c(0.025, 0.975))

# statistics computations Y
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical correlation
post_cor_Y <- apply(pred_mat_Y, 3, cor)[2,]
post_cor_Y |> quantile(c(0.025, 0.5, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(y_u_1[,1] >= post_qnt_Y[1,,1] & y_u_1[,1] <= post_qnt_Y[2,,1]),
                mean(y_u_1[,2] >= post_qnt_Y[1,,2] & y_u_1[,2] <= post_qnt_Y[2,,2]))
# mean(Y_u >= post_qnt_Y[1,,] & Y_u <= post_qnt_Y[2,,])
cat("Empirical average coverage for Spatial process:", round(mean(coverage_Y), 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
(rmspe_Y <- sqrt( colMeans( (y_u_1 - post_mean_Y)^2 ) )); mean(rmspe_Y)
(mape_Y <- colMeans( abs(y_u_1 - post_mean_Y) ) ); mean(mape_Y)

# Posterior inference -----------------------------------------------------

beta_smp  <- out500$posterior$beta
sigma_smp <- out500$posterior$sigma
K <- length(out500$weights_global)
R <- 200

post_mean_beta <-  apply(beta_smp, c(1,2), mean)
post_var_beta <- apply(beta_smp, c(1,2), sd)
post_qnt_beta <- apply(beta_smp, c(1,2), quantile, c(0.025, 0.5, 0.975))

post_mean_sigma <- apply(sigma_smp, c(1,2), mean)
post_var_sigma <- apply(sigma_smp, c(1,2), sd)
post_qnt_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025, 0.5, 0.975))

Wbps <- matrix(out500$weights_global)
(post_mean_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% out500$weights_local[[k]]) %*% Wbps)
(post_var_hyp  <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% out500$weights_local[[k]])^2 %*% Wbps - (post_mean_hyp^2))

posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                       t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3])
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]")

# Comparison results ----------------------------------------------------------

summary_table <- rbind(
  data.frame(model = "DBPS500",  t(apply(results_dbps500[, -1], 2, mean))),
  data.frame(model = "DBPS250",  t(apply(results_dbps250[, -1], 2, mean))),
  data.frame(model = "BLM",      t(apply(results_blm[,  -1], 2, mean))),
  data.frame(model = "AutoML",   t(apply(results_aml[,  -1], 2, mean))),
  data.frame(model = "DL",       t(apply(results_dl[,   -1], 2, mean))),
  data.frame(model = "DRF",      t(apply(results_drf[,  -1], 2, mean))),
  data.frame(model = "GBM",      t(apply(results_gbm[,  -1], 2, mean)))
)
print(summary_table)

# Plotting result -------------------------------------------------------------

common_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(colour = "black"),
    axis.text       = element_text(size = 16, face = "bold"),
    axis.title      = element_text(size = 16, face = "bold"),
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 16),
    legend.key.height = unit(1, "cm")
  )

x_breaks <- c(-120, -60, 0, 60, 120)
x_labels <- c("120°W", "60°W", "0°", "60°E", "120°E")
y_breaks <- c(-40, -20, 0, 20, 40, 60, 80)
y_labels <- c("40°S", "20°S", "0°", "20°N", "40°N", "60°N", "80°N")

interp_surface <- function(lon, lat, z, no.X = 300, no.Y = 300) {
  surf <- MBA::mba.surf(cbind(lon, lat, z), no.X = no.X, no.Y = no.Y,
                        exten = FALSE, sp = TRUE, h = 12)$xyz.est
  data.frame(lon = surf$x, lat = surf$y, z = as.vector(surf$z))
}

world        <- ne_countries(scale = "medium", returnclass = "sf")
world_raster <- raster::raster(raster::extent(-180, 180, -90, 90), resolution = 1)
world_raster <- raster::rasterize(world, world_raster, field = 1, background = NA)

mask_ocean <- function(df) {
  df[!is.na(raster::extract(world_raster, df[, c("lon", "lat")])), ]
}

interp_and_mask <- function(lon, lat, z) {
  mask_ocean(interp_surface(lon, lat, z))
}

plot_map <- function(df, fill_label, pal = NULL, limits = NULL, residual = FALSE) {
  if (residual) {
    sym_lim <- if (!is.null(limits)) limits else {
      lim <- max(abs(df$z), na.rm = TRUE); c(-lim, lim)
    }
    col_pal <- colorRampPalette(brewer.pal(11, "RdGy"))(101)
    p_scale <- scale_fill_gradientn(
      colors   = col_pal,
      limits   = sym_lim,
      breaks   = pretty(sym_lim, n = 5),
      name     = fill_label,
      na.value = NA
    )
  } else {
    data_range <- if (!is.null(limits)) limits else range(df$z, na.rm = TRUE)
    breaks     <- pretty(data_range, n = 5)
    col_pal    <- if (pal == "BrBG") {
      colorRampPalette(brewer.pal(9, pal))(100)
    } else {
      rev(colorRampPalette(brewer.pal(9, pal))(100))
    }
    p_scale <- scale_fill_gradientn(colors   = col_pal, breaks = breaks,
                                    limits   = limits,  name   = fill_label,
                                    na.value = NA)
  }
  
  ggplot(df, aes(x = lon, y = lat, fill = z)) +
    geom_raster(interpolate = TRUE) +
    borders("world", colour = "black", fill = NA, linewidth = 0.3) +
    p_scale +
    coord_fixed(ratio = 1.3, ylim = c(-55, 80), xlim = c(-180, 180)) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels) +
    labs(x = NULL, y = NULL) +
    common_theme
}

plot_scatter <- function(y_true, y_pred, y_lower = NULL, y_upper = NULL,
                         xlab = "True", ylab = "Predicted",
                         show_errorbars = TRUE, n_sample = 1000,
                         legend_style = "hline") {  # "hline" o "square"
  set.seed(1997)
  idx  <- sample(length(y_true), min(n_sample, length(y_true)))
  df_s <- data.frame(
    true  = y_true[idx],
    pred  = y_pred[idx],
    lower = if (!is.null(y_lower)) y_lower[idx] else NA_real_,
    upper = if (!is.null(y_upper)) y_upper[idx] else NA_real_
  )
  
  ols_slope <- coef(lm(pred ~ true - 1, data = df_s))
  ax_range  <- range(c(df_s$true, df_s$pred, df_s$lower, df_s$upper), na.rm = TRUE)
  
  df_leg <- data.frame(true = NA_real_, pred = NA_real_,
                       Line = factor(c("Identity", "OLS"), levels = c("Identity", "OLS")))
  
  p <- ggplot(df_s, aes(x = true, y = pred))
  
  if (show_errorbars && !is.null(y_lower)) {
    p <- p + geom_errorbar(aes(ymin = lower, ymax = upper),
                           width = 0, linewidth = 0.15, alpha = 0.25, color = "#D41159")
  }
  
  p <- p +
    geom_point(size = 0.8, alpha = 0.5, color = "#1A85FF") +
    # linee reali senza legenda
    geom_abline(slope = 1,         intercept = 0, linetype = "dashed",
                color = "black",   linewidth = 0.8) +
    geom_abline(slope = ols_slope, intercept = 0, linetype = "solid",
                color = "#D41159", linewidth = 0.8)
  
  if (legend_style == "hline") {
    p <- p +
      geom_line(data = df_leg,
                aes(x = true, y = pred, color = Line, linetype = Line),
                linewidth = 0.8, key_glyph = "path") +
      scale_color_manual(name = NULL,
                         values = c("Identity" = "black", "OLS" = "#D41159")) +
      scale_linetype_manual(name = NULL,
                            values = c("Identity" = "dashed", "OLS" = "solid"))
    
  } else {
    p <- p +
      geom_point(data = df_leg,
                 aes(x = true, y = pred, fill = Line),
                 shape = 22, size = 5, color = NA) +
      scale_fill_manual(name = NULL,
                        values = c("Identity" = "black", "OLS" = "#D41159"))
  }
  
  p +
    labs(x = xlab, y = ylab) +
    common_theme +
    theme(legend.position = "right",
          legend.title    = element_blank())
}

lon_u <- crd_u_1[, 1]; lat_u <- crd_u_1[, 2]

NDVI_test_df  <- interp_and_mask(lon_u, lat_u, y_u_1[, 1])
NDVI_hat_df   <- interp_and_mask(lon_u, lat_u, post_mean_Y[, 1])
NDVI_resid_df <- interp_and_mask(lon_u, lat_u, y_u_1[, 1] - post_mean_Y[, 1])
RR_test_df    <- interp_and_mask(lon_u, lat_u, y_u_1[, 2])
RR_hat_df     <- interp_and_mask(lon_u, lat_u, post_mean_Y[, 2])
RR_resid_df   <- interp_and_mask(lon_u, lat_u, y_u_1[, 2] - post_mean_Y[, 2])

lims_ndvi <- range(c(NDVI_test_df$z, NDVI_hat_df$z),  na.rm = TRUE)
lims_rr   <- range(c(RR_test_df$z,   RR_hat_df$z),    na.rm = TRUE)
lim_resid <- max(abs(c(NDVI_resid_df$z, RR_resid_df$z)), na.rm = TRUE)
resid_limits <- c(-lim_resid, lim_resid)

p_ndvi_test  <- plot_map(NDVI_test_df,  "NDVI",     pal = "BrBG", limits = lims_ndvi)
p_ndvi_hat   <- plot_map(NDVI_hat_df,   "NDVI",     pal = "BrBG", limits = lims_ndvi)
p_ndvi_resid <- plot_map(NDVI_resid_df, "Residual", residual = TRUE,  limits = resid_limits)
p_rr_test    <- plot_map(RR_test_df,    "RR",  pal = "RdYlBu", limits = lims_rr)
p_rr_hat     <- plot_map(RR_hat_df,     "RR",  pal = "RdYlBu", limits = lims_rr)
p_rr_resid   <- plot_map(RR_resid_df,   "Residual", residual = TRUE,  limits = resid_limits)

sc_ndvi <- plot_scatter(
  y_true  = y_u_1[, 1], y_pred  = post_mean_Y[, 1],
  y_lower = post_qnt_Y[1, , 1], y_upper = post_qnt_Y[2, , 1],
  xlab = "True", ylab = "Predicted",
  show_errorbars = TRUE,
  legend_style = "hline",
  n_sample = 5000
)
sc_rr <- plot_scatter(
  y_true  = y_u_1[, 2], y_pred  = post_mean_Y[, 2],
  y_lower = post_qnt_Y[1, , 2], y_upper = post_qnt_Y[2, , 2],
  xlab = "True", ylab = "Predicted",
  show_errorbars = TRUE,
  legend_style = "hline",
  n_sample = 5000
)

final_ndvi <- (p_ndvi_test | p_ndvi_hat) / (p_ndvi_resid | sc_ndvi)
final_rr   <- (p_rr_test   | p_rr_hat)  / (p_rr_resid   | sc_rr)

png("output/NDVI500_map.png",
    width = 20, height = 10, units = "in", res = 300, family = "sans")
print(final_ndvi)
dev.off()

png("outbput/RR500_map.png",
    width = 20, height = 10, units = "in", res = 300, family = "sans")
print(final_rr)
dev.off()
