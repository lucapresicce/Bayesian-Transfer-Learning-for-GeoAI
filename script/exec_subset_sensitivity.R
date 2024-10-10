## SUBSET SENSITIVITY ANALYSIS ######################################################################################################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-for-GeoAI")

# Packages --------------------------------------------------------------------
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
library(bayesplot)

# Data generation -------------------------------------------------------------

# dimensions
n <- 5000
m <- 500
u <- 250
p <- 2
q <- 2

# parameters
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- matrix(c(1, 0.40, -0.3, 1), q, q)
alfa <- 0.8
phi <- 4
set.seed(95)
# generate sintethic data
crd <- matrix(runif((n+u) * 2), ncol = 2)
X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
D <- arma_dist(crd)
gc()
Rphi <- exp(-phi * D)
rm("D"); gc()
W_or <- matrix(0, n+u, q) + mniw::rMNorm(1, Lambda = matrix(0, n+u, q), SigmaR = Rphi, SigmaC = sigma2)
rm("Rphi"); gc()
Y_or <- X_or %*% B + W_or + mniw::rMNorm(1, Lambda = matrix(0, n+u, q), SigmaR = diag((1/alfa)-1, n+u), SigmaC = sigma2)
gc()

# sample data
crd_s <- crd[1:n, ]
X <- X_or[1:n, ]
W <- W_or[1:n, ]
Y <- Y_or[1:n, ]

# prediction data
crd_u <- crd[-(1:(n-m)), ]
X_u <- X_or[-(1:(n-m)), ]
W_u <- W_or[-(1:(n-m)), ]
Y_u <- Y_or[-(1:(n-m)), ]

# Subset posterior models -----------------------------------------------------

# hyperparameters values
alfa_seq <- c(0.7, 0.8, 0.9)
phi_seq <- c(3, 4, 5)

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


# Subset sensitivity ------------------------------------------------------

res0 <- c(25, 50, 100, 250, 500, 1000, 1250) # subset sizes
res1 <- NULL                                 # result container for rmspe
res2 <- NULL                                 # result container for fit time
res3 <- NULL                                 # result container for sample iter time

pb <- txtProgressBar(0, length(res0), style = 3)
for (i in 1:length(res0)) {
  
  subset_size <- res0[i]
  K <- n/subset_size
  data_part <- subset_data(data = list(Y = Y, X = X, crd = crd_s), K = K)
  
  # BPS parallel fit -------------------------------------------------------
  
  # number of clusters for parallel implementation
  n.core <- parallel::detectCores(logical=F)-1
  
  # list of function
  funs_fit <- lsf.str()[which(lsf.str() != "fit_loop")]
  
  # list of function
  funs_pred <- lsf.str()[which(lsf.str() != "pred_loop")]
  
  # starting cluster
  cl <- makeCluster(n.core)
  registerDoParallel(cl)
  
  # timing
  tic("total")
  
  # parallelized subset computation of GP in different cores
  tic("fit")
  obj_fit <- foreach(i = 1:K, .noexport = funs_fit) %dopar% { fit_loop(i) }
  fit_time <- toc()
  
  gc()
  # Combination using double BPS
  tic("comb")
  comb_bps <- BPS_combine(obj_fit, K, 1)
  comb_time <- toc()
  Wbps <- comb_bps$W
  W_list <- comb_bps$W_list
  
  gc()
  # parallelized subset computation of GP in different cores
  R <- 250
  subset_ind <- sample(1:K, R, T, Wbps)
  tic("prediction")
  predictions <- foreach(r = 1:R, .noexport = funs_pred) %dopar% { pred_loop(r) }
  prd_time <- toc()
  
  # timing
  tot_time <- toc()
  
  # closing cluster
  stopCluster(cl)
  gc()
  
  # Results collection ----------------------------------------------------------
  
  # statistics computations Y
  pred_mat_Y <- sapply(1:R, function(r){predictions[[r]]$Pred[[1]]$Yu}, simplify = "array")
  post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
  
  # Root Mean Square Prediction Error
  rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) )
  
  # Save result -----------------------------------------------------------------
  
  res1 <- c(res1, mean(rmspe_Y))
  res2 <- c(res2, as.numeric(fit_time$toc-fit_time$tic))
  tic()
  tmp <- fit_loop(sample(1:K, 1))
  iter_time <- toc()
  res3 <- c(res3, as.numeric(iter_time$toc - iter_time$tic))
  
  setTxtProgressBar(pb = pb, value = i)
  
}

# Min-Max normalization (to better comparison)
res1_norm <- (res1-min(res1))/(max(res1)-min(res1))
res2_norm <- (res2-min(res2))/(max(res2)-min(res2))

# plotting results
width <- 360*1.5
height <- 360
pointsize <- 12
png("output/subset_sens.png", width = width, height = height, pointsize = pointsize, family = "sans")

par(mar = c(5, 4, 4, 6) + 0.1)
plot(res0, res1_norm, type = "l", col = "cornflowerblue", lwd = 2, lty = 1,
     xlab = "Subset size", ylab = "Average RMSPE", main = "", frame.plot = FALSE, axes = F,)
box(bty = "l")
axis(2)
axis(1) 
par(new = TRUE)
plot(res0, res2_norm, type = "l", col = "seagreen", lwd = 2, lty = 5,
     xlab = "", ylab = "", axes = FALSE, frame.plot = F)
axis(4)
mtext("Total Time (Sec)", side = 4, line = 3)

dev.off()



