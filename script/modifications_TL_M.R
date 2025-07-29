## MULTIVARIATE TRANSFER LEARNING SIMULATION ########################################################################################################

rm(list = ls())
gc()
# setwd(".../Bayesian-Transfer-Learning-for-GeoAI")

# Packages ----------------------------------------------------------------

library(spBPS)
library(Rcpp)
library(RcppArmadillo)
library(mniw)
library(MCMCpack)
library(ggplot2)
library(tictoc)
library(Mposterior)
library(parallel)
library(doParallel)
library(foreach)
library(MBA)
library(classInt)
library(RColorBrewer)
library(sp)
library(fields)
library(microbenchmark)
library(RcppClock)
library(bayesplot)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(patchwork)


# loading fitting function for Bayesian transfer learning
Sys.setenv(PKG_CXXFLAGS = "-Ofast"); sourceCpp("code/src/code.cpp")

# BPS fixed features
# function for the fit loop
fit_loop <- function(i) {
  
  Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
  p <- ncol(Xi); q <- ncol(Yi)
  bps <- spBPS::BPS_weights_MvT(data = list(Y = Yi, X = Xi),
                                priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                              V_r = diag(10, p),
                                              Psi = diag(1, q),
                                              nu = 3), coords = crd_i,
                                hyperpar = list(alpha = alfa_seq, phi = phi_seq), K = 5)
  w_hat <- bps$W
  epd <- bps$epd
  
  result <- list(epd, w_hat)
  return(result)
  
}

# function for the pred loop
pred_loop <- function(r) {
  
  ind_s <- subset_ind[r]
  Ys <- data_part$Y_list[[ind_s]]; Xs <- data_part$X_list[[ind_s]]; crds <- data_part$crd_list[[ind_s]]; Ws <- W_list[[ind_s]]
  result <- spBPS::BPS_post_MvT(data = list(Y = Ys, X = Xs), coords = crds,
                                X_u = X_u, crd_u = crd_u,
                                priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                              V_r = diag(10, p),
                                              Psi = diag(1, q),
                                              nu = 3),
                                hyperpar = list(alpha = alfa_seq, phi = phi_seq),
                                W = Ws, R = 1)
  
  return(result)
}

# number of clusters for parallel implementation
n.core <- parallel::detectCores(logical = F)-1

# list of function
funs_fit <- lsf.str()[which(lsf.str() != "fit_loop")]

# list of function
funs_pred <- lsf.str()[which(lsf.str() != "pred_loop")]

# Replications ------------------------------------------------------------

# fixed quantities
n <- 5000
u <- 250
p <- 2
q <- 3

# parameters
B <- matrix(c(-0.75, 2.20, 1.05, -1.1, -0.35, 0.45), p, q)
sigma2 <- matrix(c(2, 0.8, 0.2, 0.8, 2, -0.45, 0.2, -0.45, 2), q, q)
alfa <- 0.8
phi <- 4

set.seed(97)
# generate constant data structure
crd <- matrix(runif((n+u) * 2), ncol = 2)
D <- arma_dist(crd); gc()
Rphi <- exp(-phi * D)
rm("D"); gc()

# replications loop
BB <- 50
out_list <- vector(mode = "list", length = BB)
pb <- txtProgressBar(min = 0, max = BB, style = 3)
for (b in 1:BB) {
  
  cat("Computation started ... \n")
  
  # define models sub-list
  replication <- vector(mode = "list", length = 6)
  
  # Fix the seed for reproducibility
  set.seed(b)
  
  # Data generation -------------------------------------------------------
  X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
  W_or <- matrix(0, n+u, q) + mniw::rMNorm(1, Lambda = matrix(0, n+u, q), SigmaR = Rphi, SigmaC = sigma2)
  Y_or <- X_or %*% B + W_or + mniw::rMNorm(1, Lambda = matrix(0, n+u, q), SigmaR = diag((1/alfa)-1, n+u), SigmaC = sigma2)
  gc()
  
  # sample data
  crd_s <- crd[1:n, ]
  X <- X_or[1:n, ]
  W <- W_or[1:n, ]
  Y <- Y_or[1:n, ]
  
  # prediction data
  crd_u <- crd[-(1:n), ]
  X_u <- X_or[-(1:n), ]
  W_u <- W_or[-(1:n), ]
  Y_u <- Y_or[-(1:n), ]
  
  # save test data
  replication[[6]] <- Y_u
  
  # free memory
  rm(list = c("Y_or","W_or","X_or"))
  
  # data subsetting
  subset_size <- 500
  K <- n/subset_size
  data_part <- spBPS::subset_data(data = list(Y = Y, X = X, crd = crd_s), K = K)
  
  cat("Data generated ... \n")
  
  # well specified - Transfer Learning ------------------------------------
  
  # prior inizialization
  p <- ncol(X)
  priors_list <- vector("list", K+1)
  priors_list[[1]] <- list(mu_B = matrix(0, nrow = p, ncol = q),
                           V_r = diag(10, p),
                           Psi = diag(1, q),
                           nu = 3)
  R <- 250
  W_smp <-      vector("list", K)
  beta_smp <-   vector("list", K)
  sigma2_smp <- vector("list", K)
  
  # setting the true hyperparameters
  alfa_hat <- alfa
  phi_hat <- phi
  hpar <- list(alpha = alfa_hat, phi = phi_hat)
  
  tic()
  for (i in 1:K) {
    
    out <- fit_cpp_MvT(data = list(Y = data_part$Y_list[[i]], X = data_part$X_list[[i]]),
                       priors = priors_list[[i]],
                       coords = data_part$crd_list[[i]],
                       hyperpar = hpar)
    
    priors_list[[i+1]] <- list(mu_B = out$mu_star[1:p,],
                               V_r = out$V_star[1:p, 1:p],
                               Psi = out$Psi_star,
                               nu = out$nu_star)
    
    smp <- post_draws_MvT(out, R, F, p)
    
    W_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$beta[-(1:p),], simplify = "array")
    beta_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$beta[1:p,], simplify = "array")
    sigma2_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$sigma, simplify = "array")
    
  }
  fit_time_ws <- toc()
  
  # collect posteroir inference
  sigma2_mat <- matrix(sigma2_smp[[K]], nrow = R, byrow = T)
  beta_mat <- matrix(beta_smp[[K]], nrow = R, byrow = T)
  
  posterior_ws <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
  colnames(posterior_ws)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                             "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")
  
  # collect predictions
  YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
  predictions_ws <- r_pred_marg_MvT(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                                    d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                                    hyperpar = hpar, poster = out)$Yu
  
  replication[[1]] <- list("posterior"  = posterior_ws,
                           "prediction" = predictions_ws,
                           "time" = fit_time_ws)
  
  
  cat("ws - completed ... \n")
  
  # misspecified - Transfer Learning ------------------------------------
  
  # prior inizialization
  p <- ncol(X)
  priors_list <- vector("list", K+1)
  priors_list[[1]] <- list(mu_B = matrix(0, nrow = p, ncol = q),
                           V_r = diag(10, p),
                           Psi = diag(1, q),
                           nu = 3)
  R <- 250
  W_smp <-      vector("list", K)
  beta_smp <-   vector("list", K)
  sigma2_smp <- vector("list", K)
  
  # setting the true hyperparameters
  set.seed(1997)
  alfa_hat <- alfa + runif(1, -0.5, 0.2)
  phi_hat <- phi + runif(1, -10, 10)
  hpar <- list(alpha = alfa_hat, phi = phi_hat)
  
  tic()
  for (i in 1:K) {
    
    out <- fit_cpp_MvT(data = list(Y = data_part$Y_list[[i]], X = data_part$X_list[[i]]),
                       priors = priors_list[[i]],
                       coords = data_part$crd_list[[i]],
                       hyperpar = hpar)
    
    priors_list[[i+1]] <- list(mu_B = out$mu_star[1:p,],
                               V_r = out$V_star[1:p, 1:p],
                               Psi = out$Psi_star,
                               nu = out$nu_star)
    
    smp <- post_draws_MvT(out, R, F, p)
    
    W_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$beta[-(1:p),], simplify = "array")
    beta_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$beta[1:p,], simplify = "array")
    sigma2_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$sigma, simplify = "array")

  }
  fit_time_ms <- toc()
  
  # collect posteroir inference
  sigma2_mat <- matrix(sigma2_smp[[K]], nrow = R, byrow = T)
  beta_mat <- matrix(beta_smp[[K]], nrow = R, byrow = T)
  
  posterior_ms <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
  colnames(posterior_ms)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                             "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")
  
  # collect predictions
  YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
  predictions_ms <- r_pred_marg_MvT(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                                    d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                                    hyperpar = hpar, poster = out)$Yu
  
  replication[[2]] <- list("posterior"  = posterior_ms,
                           "prediction" = predictions_ms,
                           "time" = fit_time_ms)
  
  
  cat("ms - completed ... \n")
  
  # highly misspecified - Transfer Learning ------------------------------------
  
  # prior inizialization
  p <- ncol(X)
  priors_list <- vector("list", K+1)
  priors_list[[1]] <- list(mu_B = matrix(0, nrow = p, ncol = q),
                           V_r = diag(10, p),
                           Psi = diag(1, q),
                           nu = 3)
  R <- 250
  W_smp <-      vector("list", K)
  beta_smp <-   vector("list", K)
  sigma2_smp <- vector("list", K)
  
  # setting the true hyperparameters
  alfa_hat <- 0.25
  phi_hat <- 50
  hpar <- list(alpha = alfa_hat, phi = phi_hat)
  
  tic()
  for (i in 1:K) {
    
    out <- fit_cpp_MvT(data = list(Y = data_part$Y_list[[i]], X = data_part$X_list[[i]]),
                       priors = priors_list[[i]],
                       coords = data_part$crd_list[[i]],
                       hyperpar = hpar)
    
    priors_list[[i+1]] <- list(mu_B = out$mu_star[1:p,],
                               V_r = out$V_star[1:p, 1:p],
                               Psi = out$Psi_star,
                               nu = out$nu_star)
    
    smp <- post_draws_MvT(out, R, F, p)
    
    W_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$beta[-(1:p),], simplify = "array")
    beta_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$beta[1:p,], simplify = "array")
    sigma2_smp[[i]] <- sapply(1:R, function(r)smp[[r]]$sigma, simplify = "array")

  }
  fit_time_hms <- toc()
  
  # collect posteroir inference
  sigma2_mat <- matrix(sigma2_smp[[K]], nrow = R, byrow = T)
  beta_mat <- matrix(beta_smp[[K]], nrow = R, byrow = T)
  
  posterior_hms <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
  colnames(posterior_hms)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                              "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")
  
  # collect predictions
  YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
  predictions_hms <- r_pred_marg_MvT(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                                     d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                                     hyperpar = hpar, poster = out)$Yu
  
  replication[[3]] <- list("posterior"  = posterior_hms,
                           "prediction" = predictions_hms,
                           "time" = fit_time_hms)
  
  
  cat("hms - completed ... \n ")
  
  # BPS closed - Parallel Learning ------------------------------------------
  
  alfa_seq <- seq(0.75, 0.85, 0.05)
  phi_seq <- seq(8, 12, 2)

  # starting cluster
  cl <- makeCluster(n.core)  
  registerDoParallel(cl)
  
  # parallelized subset computation of GP in different cores
  tic("fit_time_bps")
  obj_fit <- foreach(i = 1:K, .noexport = funs_fit) %dopar% { fit_loop(i) }
  
  gc()
  # Combination using double BPS
  comb_bps <- BPS_combine(obj_fit, K, 1)
  fit_time_bps <- toc()
  
  Wbps <- comb_bps$W/sum(comb_bps$W)
  W_list <- comb_bps$W_list
  
  gc()
  # parallelized subset computation of GP in different cores
  R <- 250
  subset_ind <- sample(1:K, R, T, Wbps)
  predictions <- foreach(r = 1:R, .noexport = funs_pred) %dopar% { pred_loop(r) }
  
  # closing cluster
  stopCluster(cl)
  gc()
  
  # collect posteroir inference
  sigma2_mat <- matrix(sapply(1:R, function(r){predictions[[r]]$Post[[1]]$sigma}, simplify = "array"), nrow = R, byrow = T)
  beta_mat <- matrix(sapply(1:R, function(r){predictions[[r]]$Post[[1]]$beta[1:p,]}, simplify = "array"), nrow = R, byrow = T)
  
  posterior_bps <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
  colnames(posterior_bps)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                              "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")
  
  # collect predictions
  predictions_bps <- sapply(1:R, function(r){predictions[[r]]$Pred[[1]]$Yu}, simplify = "array")
  
  replication[[4]] <- list("posterior"  = posterior_bps,
                           "prediction" = predictions_bps,
                           "time" = fit_time_bps,
                           "weights" = Wbps)
  
  
  cat("bps closed - completed ... \n")
  
  # BPS open - Parallel Learning --------------------------------------------
  
  set.seed(97)
  (alfa_seq <- runif(3, 0.5, 1))
  set.seed(97)
  (phi_seq <- sample(1:50, 3, F))
  
  # starting cluster
  cl <- makeCluster(n.core)  
  registerDoParallel(cl)
  
  # parallelized subset computation of GP in different cores
  tic("fit_time_bps")
  obj_fit <- foreach(i = 1:K, .noexport = funs_fit) %dopar% { fit_loop(i) }
  
  gc()
  # Combination using double BPS
  comb_bps <- BPS_combine(obj_fit, K, 1)
  fit_time_open <- toc()
  
  Wbps <- comb_bps$W/sum(comb_bps$W)
  W_list <- comb_bps$W_list
  
  gc()
  # parallelized subset computation of GP in different cores
  R <- 250
  subset_ind <- sample(1:K, R, T, Wbps)
  predictions <- foreach(r = 1:R, .noexport = funs_pred) %dopar% { pred_loop(r) }
  
  # closing cluster
  stopCluster(cl)
  gc()
  
  # collect posteroir inference
  sigma2_mat <- matrix(sapply(1:R, function(r){predictions[[r]]$Post[[1]]$sigma}, simplify = "array"), nrow = R, byrow = T)
  beta_mat <- matrix(sapply(1:R, function(r){predictions[[r]]$Post[[1]]$beta[1:p,]}, simplify = "array"), nrow = R, byrow = T)
  
  posterior_open <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
  colnames(posterior_open)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                               "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")
  
  # collect predictions
  predictions_open <- sapply(1:R, function(r){predictions[[r]]$Pred[[1]]$Yu}, simplify = "array")
  
  replication[[5]] <- list("posterior"  = posterior_open,
                           "prediction" = predictions_open,
                           "time" = fit_time_open,
                           "weights" = Wbps)
  
  
  cat("bps open - completed ... \n")
  
  # assign models sub-list ------------------------------------------------
  out_list[[b]] <- replication
  
  # free memory
  rm(list = c("replication"))
  
  # pb
  setTxtProgressBar(pb, b)
  
  
  }

# free the environment
rm("Rphi"); gc()

# save the replications
save(out_list, file = "replications_results.RData")
# load("replications_results.RData")


# Results collection ------------------------------------------------------

# extract test data
test <- sapply(out_list, function(b){ b[[6]] }, simplify = "array")

# result collections
metrics_list <- vector(mode = "list", 5)
for (i in 1:5) {
  
  # extract results
  post <- sapply(out_list, function(b){ b[[i]]$posterior }, simplify = "array")
  pred <- sapply(out_list, function(b){ b[[i]]$prediction }, simplify = "array")

  # set the true parameters
  par_names <- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                 "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")
  true_par <- c(sigma2[upper.tri(sigma2,T)], matrix(B))
  true_mat <- matrix(true_par, nrow = length(true_par), ncol = BB)
  
  # posterior
  post_med_pars <- apply(post, c(2,3), mean)
  post_low_pars <- apply(post, c(2,3), quantile, c(0.025))
  post_upp_pars <- apply(post, c(2,3), quantile, c(0.975))
  
  # posterior metrics
  post_embias <- matrix(rowMeans(post_med_pars - true_mat), dimnames = list(par_names, "emp bias"))
  post_coverg <- matrix(rowMeans((post_low_pars < true_mat) &  (post_upp_pars > true_mat)), , dimnames = list(par_names, "coverage"))
  post_stddev <- matrix(rowMeans(apply(post, c(2,3), sd)), dimnames = list(par_names, "std dev"))
  
  # predictions
  if (i <= 3) {
    
    post_mean_Y <- apply(pred, c(1,2,4), mean)
    post_low_Y  <- apply(pred, c(1,2,4), quantile, c(0.025))
    post_upp_Y  <- apply(pred, c(1,2,4), quantile, c(0.975))
    post_var_Y  <- apply(pred, c(1,2,4), sd)
    
  } else {
    
    post_mean_Y <- apply(pred, c(1,2,5), mean)
    post_low_Y  <- apply(pred, c(1,2,5), quantile, c(0.025))
    post_upp_Y  <- apply(pred, c(1,2,5), quantile, c(0.975))
    post_var_Y  <- apply(pred, c(1,2,5), sd)
    
  }
  
  # coverage
  covrg_mat_Y <- sapply(1:BB, function(b) colMeans(test[,,b] >= post_low_Y[,,b] & test[,,b] <= post_upp_Y[,,b]))
  
  # pred int length
  leng_mat_Y <- sapply(1:BB, function(b) colMeans(post_upp_Y[,,b] - post_low_Y[,,b]))
  
  # Root Mean Square Prediction Error
  rmspe_mat <- sapply(1:BB, function(b) sqrt(colMeans( (test[,,b] - post_mean_Y[,,b])^2 )) )
  
  # Bias
  errors_mat  <- sapply(1:BB, function(b) (test[,,b] - post_mean_Y[,,b]), simplify = "array")
  bias_mat  <- sapply(1:BB, function(b) colMeans(test[,,b] - post_mean_Y[,,b]) )
  
  # Variance of errors around the bias (i.e., residual variance)
  prd_std_mat <- sapply(1:BB, function(b) colMeans(sqrt((errors_mat[,,b] - matrix(bias_mat[,b], 
                                                                                  nrow = nrow(errors_mat[,,b]), 
                                                                                  ncol = ncol(errors_mat[,,b]), 
                                                                                  byrow = TRUE))^2)) )
  
  # collect data
  metrics_list[[i]] <- list("pred_cover"  = covrg_mat_Y,
                            "pred_length" = leng_mat_Y,
                            "pred_rmpse"  = rmspe_mat,
                            "pred_bias"   = bias_mat,
                            "pred_std"    = prd_std_mat,
                            "post_bias"   = post_embias,
                            "post_cover"  = post_coverg,
                            "post_std"    = post_stddev)
  
}

post_bias_array   <- sapply(metrics_list, function(m) m$post_bias )#, simplify = "array")
post_cover_array  <- sapply(metrics_list, function(m) m$post_cover)#, simplify = "array")
post_std_array    <- sapply(metrics_list, function(m) m$post_std  )#, simplify = "array") 

pred_cover_array  <- sapply(metrics_list, function(m) m$pred_cover, simplify = "array")
pred_length_array <- sapply(metrics_list, function(m) m$pred_length, simplify = "array")
pred_rmpse_array  <- sapply(metrics_list, function(m) m$pred_rmpse, simplify = "array")
pred_bias_array   <- sapply(metrics_list, function(m) m$pred_bias, simplify = "array")
pred_std_array    <- sapply(metrics_list, function(m) m$pred_std, simplify = "array")


# plot results ------------------------------------------------------------

# Posterior

# parameter names 
param_names <- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]",
                 "Sigma[2,3]", "Sigma[3,3]", "beta[1,1]", "beta[2,1]",
                 "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")
param_labels <- c(
  expression(Sigma[1*","*1]), expression(Sigma[1*","*2]), expression(Sigma[2*","*2]),
  expression(Sigma[1*","*3]), expression(Sigma[2*","*3]), expression(Sigma[3*","*3]),
  expression(beta[1*","*1]), expression(beta[2*","*1]), expression(beta[1*","*2]),
  expression(beta[2*","*2]), expression(beta[1*","*3]), expression(beta[2*","*3])
)
names(param_labels) <- param_names

# setting names
setting_names <- c("WS", "MS", "HMS", "BPS-C", "BPS-O")

# setting colors
setting_colors <- c(
  "WS"    = "#fdd835",   # Amber/yellow — well specified
  "MS"    = "#fdae61",   # Soft orange — misspecified
  "HMS"   = "#d7191c",   # Deep red — highly misspecified
  "BPS-C" = "#1976d2",   # Strong blue — BPS closed
  "BPS-O" = "#2ecc71"    # Emerald — BPS open
)

# metric names
metric_names <- c("Avg. Empirical Bias", "Avg. Coverage", "Avg. Std. Dev.")
metric_list <- list(
  "Bias"     = post_bias_array,
  "Coverage" = post_cover_array,
  "Std"      = post_std_array)

# collect data into long format
long_df <- do.call(rbind, lapply(1:3, function(metric_index) {
  mat <- metric_list[[metric_index]]
  df <- melt(mat)
  colnames(df) <- c("Parameter", "Setting", "Value")
  df$Setting <- factor(df$Setting, levels = 1:5, labels = setting_names)
  df$Parameter <- factor(param_names[df$Parameter], levels = param_names)
  df$Metric <- metric_names[metric_index]
  df
}))

# y-axis limits per metric
y_lims <- list(
  "Bias" = c(-2, 2),
  "Coverage" = c(0, 1),
  "Std" = c(0, max(long_df$Value[long_df$Metric == "Avg. Std. Dev."]) * 1.1)
)

# plot
plot_list <- lapply(metric_names, function(metric) {
  df_sub <- subset(long_df, Metric == metric)
  
  # Define y limits
  ylim_val <- switch(metric,
                     "Avg. Empirical Bias" = c(-2, 1),
                     "Avg. Coverage"       = c(0, 1),
                     "Avg. Std. Dev."      = c(0, max(df_sub$Value) * 1.1))
  
  ggplot(df_sub, aes(x = Parameter, y = Value, fill = Setting)) +
    geom_col(position = position_dodge2(width = 0.1), width = 0.5) +  # grouped bars
    scale_fill_manual(values = setting_colors) +
    scale_x_discrete(labels = param_labels) +
    labs(x = " ", y = metric, fill = "Setting") +
    coord_cartesian(ylim = ylim_val) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "right"
    )
})

# Combine plots vertically
final_plot <- plot_list[[1]] / plot_list[[2]] / plot_list[[3]] + 
  plot_layout(heights = c(1,1,1))

width <- 360*3
height <- 360*2
pointsize <- 12
png("TL_post.png", width = width, height = height, pointsize = pointsize, family = "sans")
print(final_plot)
dev.off()


# PREDICTIVE ----

# define metrics
metric_names <- c("MSPE", "Pred Int Width", "Absolute Bias", "Variance")

# define list of arrays
array_list <- list(pred_rmpse_array^2, pred_length_array, abs(pred_bias_array), pred_std_array^2)

# setting labels
setting_labels <- c("WS", "MS", "HMS", "BPS-C", "BPS-O")

# setting colors
setting_colors <- c(
  "WS"    = "#fdd835",   # Amber/yellow — well specified
  "MS"    = "#fdae61",   # Soft orange — misspecified
  "HMS"   = "#d7191c",   # Deep red — highly misspecified
  "BPS-C" = "#1976d2",   # Strong blue — BPS closed
  "BPS-O" = "#2ecc71"    # Emerald — BPS open
)

# convert all arrays into one long data frame
long_df <- do.call(rbind, lapply(1:4, function(metric_index) {
  arr <- array_list[[metric_index]]
  df <- melt(arr)
  colnames(df) <- c("Response", "ValueID", "Setting", "Value")
  df$Setting <- factor(df$Setting, levels = 1:5, labels = setting_labels)
  df$Metric <- metric_names[metric_index]
  df
}))

# labeling for plotting
long_df <- long_df %>%
  mutate(
    Response = paste0("Y", Response),
    Group = interaction(Response, Setting, sep = " × ")
  )
long_df$Metric <- factor(long_df$Metric, levels = metric_names)

# plot
final_plot2 <- 
  ggplot(long_df, aes(x = Group, y = Value, fill = Setting)) +
  
  geom_boxplot(width = 0.35, outlier.shape = NA, color = "black", alpha = 0.75, linewidth = 0.6) +
  
  facet_wrap(~ Metric, scales = "free_y") +
  coord_cartesian(ylim = c(0, NA)) +  # Y-axis starts at 0
  scale_fill_manual(values = setting_colors, name = "Setting") +
  
  labs(x = "Response × Setting", y = NULL) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

width <- 360*3
height <- 360*2
pointsize <- 12
png("TL_pred.png", width = width, height = height, pointsize = pointsize, family = "sans")
print(final_plot2)
dev.off()

