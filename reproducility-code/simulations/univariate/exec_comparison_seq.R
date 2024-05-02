## TESTING SCRIPT ###################################################################################################################################
rm(list = ls())
gc()
setwd(".../src")

# Packages ----------------------------------------------------------------

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
library(ASMK)

# Data generation -------------------------------------------------------------

# loading fitting function for Bayesian transfer learning
Sys.setenv(PKG_CXXFLAGS = "-Ofast"); sourceCpp("code.cpp")

# dimensions
n <- 5000
m <- 1000
u <- 250
p <- 2

# parameters
B <- c(2.05, -0.55)
tau2 <- 0.5
sigma2 <- 2
delta <- tau2/sigma2
phi <- 10

set.seed(1997)
# generate sintethic data
crd <- matrix(runif((n+u) * 2), ncol = 2)
X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
D <- arma_dist(crd)
gc()
Rphi <- exp(-phi * D)
rm("D"); gc()
W_or <- matrix(0, n+u) + mniw::rmNorm(1, rep(0, n+u), sigma2*Rphi)
rm("Rphi"); gc()
Y_or <- X_or %*% B + W_or + mniw::rmNorm(1, rep(0, n+u), diag(delta*sigma2, n+u))
gc()

# sample data
crd_s <- crd[1:n, ]
X <- X_or[1:n, ]
W <- W_or[1:n, ]
Y <- matrix(Y_or[1:n, ])

# prediction data
crd_u <- crd[-(1:(n-m)), ]
X_u <- X_or[-(1:(n-m)), ]
W_u <- W_or[-(1:(n-m)), ]
Y_u <- matrix(Y_or[-(1:(n-m)), ])


# Data subsetting ---------------------------------------------------------

subset_size <- 500
K <- n/subset_size
data_part <- ASMK::subset_data(data = list(Y = Y, X = X, crd = crd_s), K = K)


# well specified - Sequantial Learning ------------------------------------

# prior inizialization
p <- ncol(X)
priors_list <- vector("list", K+1)
priors_list[[1]] <- list(mu_b = rep(0, p),
                         V_b = diag(10, p),
                         a = 2,
                         b = 2)
W_smp <- NULL
beta_smp <- NULL
sigma2_smp <- NULL

# setting the true hyperparameters
delta_hat <- delta
phi_hat <- phi
hpar <- list(delta = delta_hat, phi = phi_hat)

tic()
pb <- txtProgressBar(0, K, style = 3)
for (i in 1:K) {
  
  out <- fit_cpp(data = list(Y = data_part$Y_list[[i]], X = data_part$X_list[[i]]),
                 priors = priors_list[[i]],
                 coords = data_part$crd_list[[i]],
                 hyperpar = hpar)
  
  priors_list[[i+1]] <- list(mu_b = out$gamma_hat[1:p,],
                             V_b = out$M_star[1:p, 1:p],
                             a = out$a_star,
                             b = out$b_star)
  
  smp <- post_draws(out, 250, F, p)
  
  W_smp <- cbind(W_smp, smp$Betas[, -(1:p)])
  beta_smp <- rbind(beta_smp, smp$Betas[, (1:p)])
  sigma2_smp <- rbind(sigma2_smp, smp$Sigmas)
  
  setTxtProgressBar(pb, i)
  
}
fit_time_ws <- toc()

# compute MAP estimates
betahat <- beta_smp |> colMeans()
sigma2hat <- sigma2_smp |> mean()
What <- W_smp |> colMeans()
ordW <- order(do.call(rbind, data_part$sets))
What <- What[ordW]

# plot posterior inference results
true_par <- c(sigma2, B)
posterior_ws <- cbind(sigma2_smp, beta_smp)[(250*(K-1)+1):(250*K),]
colnames(posterior_ws) <- c("sigma^2", "beta[1]", "beta[2]")

# plotting asmk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("red")
postdens_ws <- mcmc_areas(posterior_ws,
                          prob = 0.95,
                          point_est = "mean") + 
  ggtitle("Posterior distributions - well specified - Sequantial Learning",
          "with posterior means, and 95% credible intervals") +
  scale_y_discrete(labels = c(expression(sigma^2), expression(beta[1]), expression(beta[2])))

color_scheme_set("viridis")
postint_ws <- mcmc_recover_intervals(posterior_ws, true_par,
                                     prob = 0.89,
                                     prob_outer = 0.95,
                                     point_est = "mean",
                                     size = 4,
                                     alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(sigma^2))) + 
  ggtitle("Credible intervals - well specified - Sequential Learning",
          "with posterior means, true values, and 89% - 95% credible intervals")


# predictions
YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
tic()
predictions <- r_pred_marg(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                           d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                           hyperpar = hpar, poster = out)
toc()

# statistics computations Y
pred_mat_Y <- predictions$Y_u
post_mean_Y <- colMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))

# Root Mean Square Prediction Error
# (rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_u, Y_u), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Yhat <- MBA::mba.surf(cbind(crd_u, post_mean_Y), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y$z, surf.Yhat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y$z))
zlimh <- range(surf.Yhat$z)

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("seq_pred_well.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 2))
iw <- as.image.SpatialGridDataFrame(surf.Y)
iy <- as.image.SpatialGridDataFrame(surf.Yhat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True spatial process") 
axis(2, las=1)
axis(1)
image.plot(iw, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated spatial process")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimh)
dev.off()

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

png("seq_post_well.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(postdens_ws, postint_ws)
dev.off()


# misspecified - Sequantial Learning ------------------------------------

# prior inizialization
p <- ncol(X)
priors_list <- vector("list", K+1)
priors_list[[1]] <- list(mu_b = rep(0, p),
                         V_b = diag(10, p),
                         a = 2,
                         b = 2)
W_smp <- NULL
beta_smp <- NULL
sigma2_smp <- NULL

# setting the true hyperparameters
set.seed(1997)
delta_hat <- delta + runif(1, -0.25, 0.25)
phi_hat <- phi + runif(1, -10, 10)
hpar <- list(delta = delta_hat, phi = phi_hat)

tic()
pb <- txtProgressBar(0, K, style = 3)
for (i in 1:K) {
  
  out <- fit_cpp(data = list(Y = data_part$Y_list[[i]], X = data_part$X_list[[i]]),
                 priors = priors_list[[i]],
                 coords = data_part$crd_list[[i]],
                 hyperpar = hpar)
  
  priors_list[[i+1]] <- list(mu_b = out$gamma_hat[1:p,],
                             V_b = out$M_star[1:p, 1:p],
                             a = out$a_star,
                             b = out$b_star)
  
  smp <- post_draws(out, 250, F, p)
  
  W_smp <- cbind(W_smp, smp$Betas[, -(1:p)])
  beta_smp <- rbind(beta_smp, smp$Betas[, (1:p)])
  sigma2_smp <- rbind(sigma2_smp, smp$Sigmas)
  
  setTxtProgressBar(pb, i)
  
}
fit_time_ms <- toc()

# compute MAP estimates
betahat <- beta_smp |> colMeans()
sigma2hat <- sigma2_smp |> mean()
What <- W_smp |> colMeans()
ordW <- order(do.call(rbind, data_part$sets))
What <- What[ordW]

# plot posterior inference results
true_par <- c(sigma2, B)
posterior_ms <- cbind(sigma2_smp, beta_smp)[(250*(K-1)+1):(250*K),]
colnames(posterior_ms) <- c("sigma^2", "beta[1]", "beta[2]")

# plotting asmk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("red")
postdens_ms <- mcmc_areas(posterior_ms,
                          prob = 0.95,
                          point_est = "mean") + 
  ggtitle("Posterior distributions - misspecified - Sequantial Learning",
          "with posterior means, and 95% credible intervals") +
  scale_y_discrete(labels = c(expression(sigma^2), expression(beta[1]), expression(beta[2])))

color_scheme_set("viridis")
postint_ms <- mcmc_recover_intervals(posterior_ms, true_par,
                                     prob = 0.89,
                                     prob_outer = 0.95,
                                     point_est = "mean",
                                     size = 4,
                                     alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(sigma^2))) + 
  ggtitle("Credible intervals - misspecified - Sequential Learning",
          "with posterior means, true values, and 89% - 95% credible intervals")


# predictions
YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
tic()
predictions <- r_pred_marg(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                           d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                           hyperpar = hpar, poster = out)
toc()

# statistics computations Y
pred_mat_Y <- predictions$Y_u
post_mean_Y <- colMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))

# Root Mean Square Prediction Error
# (rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_u[1:m,], Y_u[1:m]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Yhat <- MBA::mba.surf(cbind(crd_u[1:m,], post_mean_Y[1:m]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y$z, surf.Yhat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y$z))
zlimh <- range(surf.Yhat$z)

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("seq_pred_mis.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 2))
iw <- as.image.SpatialGridDataFrame(surf.Y)
iy <- as.image.SpatialGridDataFrame(surf.Yhat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True spatial process") 
axis(2, las=1)
axis(1)
image.plot(iw, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated spatial process")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimh)
dev.off()

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

png("seq_post_mis.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(postdens_ms, postint_ms)
dev.off()

# highly misspecified - Sequantial Learning ------------------------------------

# prior inizialization
p <- ncol(X)
priors_list <- vector("list", K+1)
priors_list[[1]] <- list(mu_b = rep(0, p),
                         V_b = diag(10, p),
                         a = 2,
                         b = 2)
W_smp <- NULL
beta_smp <- NULL
sigma2_smp <- NULL

# setting the true hyperparameters
delta_hat <- 0.95
phi_hat <- 50
hpar <- list(delta = delta_hat, phi = phi_hat)

tic()
pb <- txtProgressBar(0, K, style = 3)
for (i in 1:K) {
  
  out <- fit_cpp(data = list(Y = data_part$Y_list[[i]], X = data_part$X_list[[i]]),
                 priors = priors_list[[i]],
                 coords = data_part$crd_list[[i]],
                 hyperpar = hpar)
  
  priors_list[[i+1]] <- list(mu_b = out$gamma_hat[1:p,],
                             V_b = out$M_star[1:p, 1:p],
                             a = out$a_star,
                             b = out$b_star)
  
  smp <- post_draws(out, 250, F, p)
  
  W_smp <- cbind(W_smp, smp$Betas[, -(1:p)])
  beta_smp <- rbind(beta_smp, smp$Betas[, (1:p)])
  sigma2_smp <- rbind(sigma2_smp, smp$Sigmas)
  
  setTxtProgressBar(pb, i)
  
}
fit_time_hms <- toc()

# compute MAP estimates
betahat <- beta_smp |> colMeans()
sigma2hat <- sigma2_smp |> mean()
What <- W_smp |> colMeans()
ordW <- order(do.call(rbind, data_part$sets))
What <- What[ordW]

# plot posterior inference results
true_par <- c(sigma2, B)
posterior_hms <- cbind(sigma2_smp, beta_smp)[(250*(K-1)+1):(250*K),]
colnames(posterior_hms) <- c("sigma^2", "beta[1]", "beta[2]")

# plotting asmk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("red")
postdens_hms <- mcmc_areas(posterior_hms,
                          prob = 0.95,
                          point_est = "mean") + 
  ggtitle("Posterior distributions - highly misspecified - Sequantial Learning",
          "with posterior means, and 95% credible intervals") +
  scale_y_discrete(labels = c(expression(sigma^2), expression(beta[1]), expression(beta[2])))

color_scheme_set("viridis")
postint_hms <- mcmc_recover_intervals(posterior_hms, true_par,
                                     prob = 0.89,
                                     prob_outer = 0.95,
                                     point_est = "mean",
                                     size = 4,
                                     alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(sigma^2))) + 
  ggtitle("Credible intervals - highly misspecified - Sequential Learning",
          "with posterior means, true values, and 89% - 95% credible intervals")


# predictions
YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
tic()
predictions <- r_pred_marg(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                           d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                           hyperpar = hpar, poster = out)
toc()

# statistics computations Y
pred_mat_Y <- predictions$Y_u
post_mean_Y <- colMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))

# Root Mean Square Prediction Error
# (rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_u[1:m,], Y_u[1:m]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Yhat <- MBA::mba.surf(cbind(crd_u[1:m,], post_mean_Y[1:m]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y$z, surf.Yhat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y$z))
zlimh <- range(surf.Yhat$z)

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("seq_pred_hms.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 2))
iw <- as.image.SpatialGridDataFrame(surf.Y)
iy <- as.image.SpatialGridDataFrame(surf.Yhat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True spatial process") 
axis(2, las=1)
axis(1)
image.plot(iw, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated spatial process")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimh)
dev.off()

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

png("seq_post_hms.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(postdens_hms, postint_hms)
dev.off()


# ASMK - Parallel Learning ------------------------------------------------

delta_seq <- seq(0.1, 0.4, 0.15)
phi_seq <- seq(8, 12, 2)

# function for the fit loop
fit_loop <- function(i) {
  
  Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
  p <- ncol(Xi)
  bps <- ASMK::BPS_weights(data = list(Y = Yi, X = Xi),
                           priors = list(mu_b = matrix(rep(0, p)),
                                         V_b = diag(10, p),
                                         a = 2,
                                         b = 2), coords = crd_i,
                           hyperpar = list(delta = delta_seq, phi = phi_seq), K = 5)
  w_hat <- bps$W
  epd <- bps$epd
  
  result <- list(epd, w_hat)
  return(result)
  
}

# function for the pred loop
pred_loop <- function(r) {
  
  ind_s <- subset_ind[r]
  Ys <- matrix(data_part$Y_list[[ind_s]]); Xs <- data_part$X_list[[ind_s]]; crds <- data_part$crd_list[[ind_s]]; Ws <- W_list[[ind_s]]
  result <- ASMK::BPS_post(data = list(Y = Ys, X = Xs), coords = crds,
                           X_u = X_u, crd_u = crd_u,
                           priors = list(mu_b = matrix(rep(0, p)),
                                         V_b = diag(10, p),
                                         a = 2,
                                         b = 2),
                           hyperpar = list(delta = delta_seq, phi = phi_seq),
                           W = Ws, R = 1)
  
  return(result)
}

# number of clusters for parallel implementation
n.core <- parallel::detectCores(logical = F)-1

# list of function
funs_fit <- lsf.str()[which(lsf.str() != "fit_loop")]

# list of function
funs_pred <- lsf.str()[which(lsf.str() != "pred_loop")]

# starting cluster
cl <- makeCluster(n.core)  
registerDoParallel(cl)

# parallelized subset computation of GP in different cores
tic("fit_time_asmk")
obj_fit <- foreach(i = 1:K, .noexport = funs_fit) %dopar% { fit_loop(i) }

gc()
# Combination using double BPS
comb_bps <- BPS_combine(obj_fit, K, 1)
fit_time_asmk <- toc()

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

# expected values of hyperparamters
sapply(1:K, function(k) t(W_list[[k]]) %*% expand_grid_cpp(delta_seq, phi_seq)) %*% Wbps 

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]][[2]]})
post_mean_Y <- rowMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))
(CI_avlen_asmk <- mean(post_qnt_Y[2,]-post_qnt_Y[1,]))

# Root Mean Square Prediction Error
# (rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_u[1:m,], Y_u[1:m]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Yhat <- MBA::mba.surf(cbind(crd_u[1:m,], post_mean_Y[1:m]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y$z, surf.Yhat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y$z))
zlimh <- range(surf.Yhat$z)

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("seq_pred_amsk.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 2))
iw <- as.image.SpatialGridDataFrame(surf.Y)
iy <- as.image.SpatialGridDataFrame(surf.Yhat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True spatial process") 
axis(2, las=1)
axis(1)
image.plot(iw, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated spatial process")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimh)
dev.off()

# plot posterior inference results
true_par <- c(sigma2, B)
posterior_asmk <- t(sapply(1:R, function(r){c(predictions[[r]][[4]], predictions[[r]][[3]][,1:p])}))
colnames(posterior_asmk) <- c("sigma^2", "beta[1]", "beta[2]")

# plotting asmk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("red")
postdens_asmk <- mcmc_areas(posterior_asmk,
                            prob = 0.95,
                            point_est = "mean") + 
  ggtitle("Posterior distributions - ASMK - Parallel Learning",
          "with posterior means, and 95% credible intervals") +
  scale_y_discrete(labels = c(expression(sigma^2), expression(beta[1]), expression(beta[2])))

color_scheme_set("viridis")
postint_asmk <- mcmc_recover_intervals(posterior_asmk, true_par,
                                       prob = 0.89,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 4,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(sigma^2))) + 
  ggtitle("Credible intervals - ASMK - Parallel Learning",
          "with posterior means, true values, and 89% - 95% credible intervals")

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

png("seq_post_asmk.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(postdens_asmk, postint_asmk)
dev.off()
