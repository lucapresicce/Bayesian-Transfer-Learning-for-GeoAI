# REPLICATION FOR SIMULATION - NNGP VS DBPS 

rm(list = ls())
gc()

# Packages --------------------------------------------------------------------

library(spBPS)
library(spBayes)
library(spNNGP)
library(Rcpp)
library(RcppArmadillo)
library(mniw)
library(MCMCpack)
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


# Replications ------------------------------------------------------------

B <- 50
results <- array(0, c(10, 5, B)) # models - metrics - replications
pb <- txtProgressBar(0, B ,0, style = 3)
for (b in 1:B) {
  # Data generation ---------------------------------------
  
  set.seed(4-8-15-16-23-42+b)
  N <- 2250
  p <- 2
  q <- 1
  coords <- cbind(runif(N), runif(N))
  X <- cbind(1, rnorm(N))   # intercept + covariate
  beta_true <- c(1.0, 0.5)
  sigma2 <- 1.0; phi <- 4; tau2 <- 0.25
  
  # distance matrix (might be heavy for 10k; fields::rdist is efficient)
  D <- rdist(coords)
  Sigma <- sigma2 * exp(-D/phi)   # exponential
  w <- as.numeric(t(chol(Sigma)) %*% rnorm(N))  # simulate GP
  eps <- rnorm(N, 0, sqrt(tau2))
  y <- as.numeric(X %*% beta_true + w + eps)
  
  # train test data
  n <- 2000
  u <- N - n
  
  crd <- coords[1:n,]
  Y <- y[1:n]
  x <- X[1:n,]
  data <- data.frame(y=y, x=X[,2], lon=coords[,1], lat=coords[,2])
  
  crdu <- coords[-(1:n),]
  Yu <- y[-(1:n)]
  xu <- X[-(1:n),]
  
  train <- data[1:n,]
  test <- data[-(1:n),]
  
  
  # NNGP (m = 5) ------------------------------------------
  n.samples <- 5000
  starting <- list("phi"=1, "sigma.sq"=1, "tau.sq"=0.5)
  phi.min <- 0.1
  phi.max <- 50
  priors <- list("phi.Unif" = c(phi.min, phi.max),
                 "sigma.sq.IG" = c(2, 1),
                 "tau.sq.IG"   = c(2, 0.5))
  tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
  cov.model <- "exponential"
  
  tictoc::tic()
  m.s <- spNNGP(Y~x-1, coords=crd, starting=starting, method="latent", n.neighbors=5,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads = 5)
  fit_time <- tictoc::toc()
  
  tictoc::tic()
  s.pred <- predict(m.s, X.0 = xu, coords.0 = crdu,
                    n.omp.threads = 5, n.report = 50)
  pred_time <- tictoc::toc()
  gc()
  
  y.hat <- apply(s.pred$p.y.0, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[1,,b] <- c(avg_width, var_width, mse, cover , as.numeric(fit_time$toc-fit_time$tic) + as.numeric(pred_time$toc-pred_time$tic))
  
  
  # NNGP (m = 10) -----------------------------------------
  n.samples <- 5000
  starting <- list("phi"=1, "sigma.sq"=1, "tau.sq"=0.5)
  phi.min <- 0.1
  phi.max <- 50
  priors <- list("phi.Unif" = c(phi.min, phi.max),
                 "sigma.sq.IG" = c(2, 1),
                 "tau.sq.IG"   = c(2, 0.5))
  tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
  cov.model <- "exponential"
  
  tictoc::tic()
  m.s <- spNNGP(Y~x-1, coords=crd, starting=starting, method="latent", n.neighbors=10,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads = 5)
  fit_time <- tictoc::toc()
  
  tictoc::tic()
  s.pred <- predict(m.s, X.0 = xu, coords.0 = crdu,
                    n.omp.threads = 5, n.report = 50)
  pred_time <- tictoc::toc()
  gc()
  
  y.hat <- apply(s.pred$p.y.0, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[2,,b] <- c(avg_width, var_width, mse, cover , as.numeric(fit_time$toc-fit_time$tic) + as.numeric(pred_time$toc-pred_time$tic))
  
  
  # NNGP (m = 20) -----------------------------------------
  n.samples <- 5000
  starting <- list("phi"=1, "sigma.sq"=1, "tau.sq"=0.5)
  phi.min <- 0.1
  phi.max <- 50
  priors <- list("phi.Unif" = c(phi.min, phi.max),
                 "sigma.sq.IG" = c(2, 1),
                 "tau.sq.IG"   = c(2, 0.5))
  tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
  cov.model <- "exponential"
  
  tictoc::tic()
  m.s <- spNNGP(Y~x-1, coords=crd, starting=starting, method="latent", n.neighbors=20,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads = 5)
  fit_time <- tictoc::toc()
  
  tictoc::tic()
  s.pred <- predict(m.s, X.0 = xu, coords.0 = crdu,
                    n.omp.threads = 5, n.report = 50)
  pred_time <- tictoc::toc()
  gc()
  
  y.hat <- apply(s.pred$p.y.0, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[3,,b] <- c(avg_width, var_width, mse, cover , as.numeric(fit_time$toc-fit_time$tic) + as.numeric(pred_time$toc-pred_time$tic))
  
  # DBPS (K = 5) ------------------------------------------
  
  # hyperparameters values
  alfa_seq <- c(0.7, 0.8, 0.9)
  phi_seq <- c(3, 4, 5)
  hyperpar <- list(alpha = alfa_seq, phi = phi_seq)
  
  # fitting DBPS
  n.core <- parallel::detectCores(logical = F)-1
  tic("spBPS")
  out <- spBPS::spBPS(data = list(Y = matrix(Y), X = x),
                      priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                    V_r = diag(10, p),
                                    Psi = diag(1, q),
                                    nu = 3),
                      coords = crd,
                      hyperpar = hyperpar,
                      subset_size = n/5,
                      cv_folds = 5,
                      rp = 1,
                      combine_method = "pseudoBMA",
                      draws = 200,
                      newdata = list(X = xu, coords = crdu),
                      include_latent = F,
                      cores = n.core)
  new_time <- toc()
  gc()
  
  pred_mat_Y <- sapply(1:200, function(r){out$predictive[[r]]$Yu})

  y.hat <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[4,,b] <- c(avg_width, var_width, mse, cover , as.numeric(new_time$toc-new_time$tic))
  
  # correcting fro disagreement
  pred_mat_Y <- sapply(1:200, function(r){out$predictive[[r]]$Yu})
  disagree_array <- sapply(1:200, function(r){out$predictive[[r]]$MY})
  pred_mat_Y <- pred_mat_Y - disagree_array
  agreement <- matrix(rowMeans(disagree_array))
  agree_array <- replicate(200, agreement, simplify = "matrix")
  pred_mat_Y <- pred_mat_Y + agree_array
  y.hat <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[5,,b] <- c(avg_width, var_width, mse, cover , as.numeric(new_time$toc-new_time$tic))
  
  
  # DBPS (K = 10) ------------------------------------------
  
  # hyperparameters values
  alfa_seq <- c(0.7, 0.8, 0.9)
  phi_seq <- c(3, 4, 5)
  hyperpar <- list(alpha = alfa_seq, phi = phi_seq)
  
  # fitting DBPS
  n.core <- parallel::detectCores(logical = F)-1
  tic("spBPS")
  out <- spBPS::spBPS(data = list(Y = matrix(Y), X = x),
                      priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                    V_r = diag(10, p),
                                    Psi = diag(1, q),
                                    nu = 3),
                      coords = crd,
                      hyperpar = hyperpar,
                      subset_size = n/10,
                      cv_folds = 5,
                      rp = 1,
                      combine_method = "pseudoBMA",
                      draws = 200,
                      newdata = list(X = xu, coords = crdu),
                      include_latent = F,
                      cores = n.core)
  new_time <- toc()
  gc()
  
  pred_mat_Y <- sapply(1:200, function(r){out$predictive[[r]]$Yu})
  
  y.hat <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[6,,b] <- c(avg_width, var_width, mse, cover , as.numeric(new_time$toc-new_time$tic))
  
  # correcting fro disagreement
  pred_mat_Y <- sapply(1:200, function(r){out$predictive[[r]]$Yu})
  disagree_array <- sapply(1:200, function(r){out$predictive[[r]]$MY})
  pred_mat_Y <- pred_mat_Y - disagree_array
  agreement <- matrix(rowMeans(disagree_array))
  agree_array <- replicate(200, agreement, simplify = "matrix")
  pred_mat_Y <- pred_mat_Y + agree_array
  y.hat <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[7,,b] <- c(avg_width, var_width, mse, cover , as.numeric(new_time$toc-new_time$tic))
  
  # DBPS (K = 20) ------------------------------------------
  
  # hyperparameters values
  alfa_seq <- c(0.70, 0.8, 0.9)
  phi_seq <- c(3, 4, 5)
  hyperpar <- list(alpha = alfa_seq, phi = phi_seq)
  
  # fitting DBPS
  n.core <- parallel::detectCores(logical = F)-1
  tic("spBPS")
  out <- spBPS::spBPS(data = list(Y = matrix(Y), X = x),
                      priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                    V_r = diag(10, p),
                                    Psi = diag(1, q),
                                    nu = 3),
                      coords = crd,
                      hyperpar = hyperpar,
                      subset_size = n/20,
                      cv_folds = 5,
                      rp = 1,
                      combine_method = "pseudoBMA",
                      draws = 200,
                      newdata = list(X = xu, coords = crdu),
                      include_latent = F,
                      cores = n.core)
  new_time <- toc()
  
  pred_mat_Y <- sapply(1:200, function(r){out$predictive[[r]]$Yu})
  
  y.hat <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[8,,b] <- c(avg_width, var_width, mse, cover , as.numeric(new_time$toc-new_time$tic))
  
  # correcting fro disagreement
  pred_mat_Y <- sapply(1:200, function(r){out$predictive[[r]]$Yu})
  disagree_array <- sapply(1:200, function(r){out$predictive[[r]]$MY})
  pred_mat_Y <- pred_mat_Y - disagree_array
  agreement <- matrix(rowMeans(disagree_array))
  agree_array <- replicate(200, agreement, simplify = "matrix")
  pred_mat_Y <- pred_mat_Y + agree_array
  y.hat <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[9,,b] <- c(avg_width, var_width, mse, cover , as.numeric(new_time$toc-new_time$tic))
  
  gc()
  # Full GP -----------------------------------------------
  
  n.samples <- 500
  starting <- list("phi"=1, "sigma.sq"=1, "tau.sq"=0.5)
  phi.min <- 0.1
  phi.max <- 50
  priors <- list("phi.Unif" = c(phi.min, phi.max),
                 "sigma.sq.IG" = c(2, 1),
                 "tau.sq.IG"   = c(2, 0.5))
  tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
  cov.model <- "exponential"
  
  tictoc::tic()
  m.s <- spLM(Y~x-1, coords=crd, starting=starting, 
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, verbose=T)
  fit_time <- tictoc::toc()
  
  tictoc::tic()
  s.pred <- spBayes::spPredict(m.s, pred.covars = xu, pred.coords = crdu,
                               n.omp.threads = 5)
  pred_time <- tictoc::toc()
  gc()
  
  y.hat <- apply(s.pred$p.y.predictive.samples, 1, quantile, c(0.025, 0.5, 0.975))
  avg_width <- mean(y.hat[3,]-y.hat[1,])
  var_width <- var(y.hat[3,]-y.hat[1,])
  mse <- mean((Yu- y.hat[2,])^2)
  cover <- mean(Yu >= y.hat[1,] & Yu <= y.hat[3,])
  
  results[10,,b] <- c(avg_width, var_width, mse, cover , as.numeric(fit_time$toc-fit_time$tic) + as.numeric(pred_time$toc-pred_time$tic))
  

  # check lap ---------------------------------------------------------------
  setTxtProgressBar(pb, b)
  
}


# Results -----------------------------------------------------------------

raw_tab <- apply(results, c(1,2), mean)
rownames(raw_tab) <- c("NNGP (m = 5)", "NNGP (m = 10)", "NNGP (m = 20)",
                       "DBPS (K = 5)", "DBPS (K = 5) DT", "DBPS (K = 10)", "DBPS (K = 10) DT", "DBPS (K = 20)", "DBPS (K = 20) DT",
                       "Full GP")
colnames(raw_tab) <- c("Avg. Pred. int. width", "Avg. Pred. int. widht var", "RMSPE", "Avg. emp. coverage", "Time (sec)")
raw_tab

# save.image("replicationsNNGPDBPS.RData")
# load("~/replicationsNNGPDBPS.RData")