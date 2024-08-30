## MULTIVARIATE TRANSFER LEARNING SIMULATION ########################################################################################################

rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-for-GeoAI")

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

# Data generation -------------------------------------------------------------

# loading fitting function for Bayesian transfer learning
Sys.setenv(PKG_CXXFLAGS = "-Ofast"); sourceCpp("code/src/code.cpp")

# dimensions
n <- 5000
m <- 500
u <- 250
p <- 2
q <- 3

# parameters
B <- matrix(c(-0.75, 2.20, 1.05, -1.1, -0.35, 0.45), p, q)
sigma2 <- matrix(c(2, 0.8, 0.2, 0.8, 2, -0.45, 0.2, -0.45, 2), q, q)
alfa <- 0.8
phi <- 4

set.seed(97)
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


# Data subsetting ---------------------------------------------------------

subset_size <- 500
K <- n/subset_size
data_part <- spBPS::subset_data(data = list(Y = Y, X = X, crd = crd_s), K = K)

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
pb <- txtProgressBar(0, K, style = 3)
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
  
  setTxtProgressBar(pb, i)
  
}
fit_time_ws <- toc()

# compute MAP estimates
betahat <- beta_smp[[K]] |> apply(c(1,2), mean)
sigma2hat <- sigma2_smp[[K]] |> apply(c(1,2), mean)
What <- do.call(rbind, lapply(1:K, function(k) W_smp[[k]] |> apply(c(1,2), mean)))
ordW <- order(do.call(rbind, data_part$sets))
What <- What[ordW,]

# plot posterior inference results
true_par <- c(sigma2[upper.tri(sigma2,T)], matrix(B))
sigma2_mat <- matrix(sigma2_smp[[K]], nrow = R, byrow = T)
beta_mat <- matrix(beta_smp[[K]], nrow = R, byrow = T)
posterior_ws <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
colnames(posterior_ws)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                           "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
postint_ws <- mcmc_recover_intervals(posterior_ws, true_par,
                                     prob = 0.95,
                                     prob_outer = 0.95,
                                     point_est = "mean",
                                     size = 4,
                                     alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[list(1,1)]), expression(beta[list(2,1)]), expression(beta[list(1,2)]),
                              expression(beta[list(2,2)]), expression(beta[list(1,3)]), expression(beta[list(2,3)]),
                              expression(Sigma[list(1,1)]), expression(Sigma[list(2,1)]), expression(Sigma[list(1,3)]),
                              expression(Sigma[list(2,2)]), expression(Sigma[list(2,3)]), expression(Sigma[list(3,3)]))) + 
  ggtitle("Credible intervals - well specified - Transfer Learning",
          "with posterior means, true values, and 95% credible intervals")


# predictions
YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
tic()
predictions_ws <- r_pred_marg_MvT(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                                              d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                                              hyperpar = hpar, poster = out)$Yu
toc()

# statistics computations Y
pred_mat_Y <- predictions_ws
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_low_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025))
post_upp_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.975))

# Empirical coverage for Y
coverage_Y <- colMeans(Y_u >= post_low_Y & Y_u <= post_upp_Y)
cat("Empirical average coverage for Response:", round(mean(coverage_Y), 3))

# Root Mean Square Prediction Error
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) ))

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y1 <- MBA::mba.surf(cbind(crd_u, Y_u[,1]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y1hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,1]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2 <- MBA::mba.surf(cbind(crd_u, Y_u[,2]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,2]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3 <- MBA::mba.surf(cbind(crd_u, Y_u[,3]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,3]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y1$z, surf.Y1hat$z, surf.Y2$z, surf.Y2hat$z, surf.Y3$z, surf.Y3hat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y1$z, surf.Y2$z, surf.Y3$z))
zlimp <- range(c(surf.Y1hat$z, surf.Y2hat$z, surf.Y3hat$z))

# Size for the mapping
width <- 360*3
height <- 360*2
pointsize <- 12

png("output/surface_M_TLwell.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(2, 3))

iy1 <- as.image.SpatialGridDataFrame(surf.Y1)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 1st response") 
axis(2, las=1)
axis(1)
image.plot(iy1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy2 <- as.image.SpatialGridDataFrame(surf.Y2)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 2nd response") 
axis(2, las=1)
axis(1)
image.plot(iy2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy3 <- as.image.SpatialGridDataFrame(surf.Y3)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 3rd response") 
axis(2, las=1)
axis(1)
image.plot(iy3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iyp1 <- as.image.SpatialGridDataFrame(surf.Y1hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 1st response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[1], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp2 <- as.image.SpatialGridDataFrame(surf.Y2hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 2nd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[2], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp3 <- as.image.SpatialGridDataFrame(surf.Y3hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 3rd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[3], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

dev.off()

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("output/CIpost_M_TLwell.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(postint_ws)
dev.off()


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
pb <- txtProgressBar(0, K, style = 3)
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
  
  setTxtProgressBar(pb, i)
  
}
fit_time_ms <- toc()

# compute MAP estimates
betahat <- beta_smp[[K]] |> apply(c(1,2), mean)
sigma2hat <- sigma2_smp[[K]] |> apply(c(1,2), mean)
What <- do.call(rbind, lapply(1:K, function(k) W_smp[[k]] |> apply(c(1,2), mean)))
ordW <- order(do.call(rbind, data_part$sets))
What <- What[ordW,]

# plot posterior inference results
true_par <- c(sigma2[upper.tri(sigma2,T)], matrix(B))
sigma2_mat <- matrix(sigma2_smp[[K]], nrow = R, byrow = T)
beta_mat <- matrix(beta_smp[[K]], nrow = R, byrow = T)
posterior_ms <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
colnames(posterior_ms)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                           "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
postint_ms <- mcmc_recover_intervals(posterior_ms, true_par,
                                     prob = 0.95,
                                     prob_outer = 0.95,
                                     point_est = "mean",
                                     size = 4,
                                     alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[list(1,1)]), expression(beta[list(2,1)]), expression(beta[list(1,2)]),
                              expression(beta[list(2,2)]), expression(beta[list(1,3)]), expression(beta[list(2,3)]),
                              expression(Sigma[list(1,1)]), expression(Sigma[list(2,1)]), expression(Sigma[list(1,3)]),
                              expression(Sigma[list(2,2)]), expression(Sigma[list(2,3)]), expression(Sigma[list(3,3)]))) + 
  ggtitle("Credible intervals - misspecified - Transfer Learning",
          "with posterior means, true values, and 95% credible intervals")


# predictions
YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
tic()
predictions_ms <- r_pred_marg_MvT(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                                  d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                                  hyperpar = hpar, poster = out)$Yu
toc()

# statistics computations Y
pred_mat_Y <- predictions_ms
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_low_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025))
post_upp_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.975))

# Empirical coverage for Y
coverage_Y <- colMeans(Y_u >= post_low_Y & Y_u <= post_upp_Y)
cat("Empirical average coverage for Response:", round(mean(coverage_Y), 3))

# Root Mean Square Prediction Error
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) ))

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y1 <- MBA::mba.surf(cbind(crd_u, Y_u[,1]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y1hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,1]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2 <- MBA::mba.surf(cbind(crd_u, Y_u[,2]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,2]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3 <- MBA::mba.surf(cbind(crd_u, Y_u[,3]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,3]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y1$z, surf.Y1hat$z, surf.Y2$z, surf.Y2hat$z, surf.Y3$z, surf.Y3hat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y1$z, surf.Y2$z, surf.Y3$z))
zlimp <- range(c(surf.Y1hat$z, surf.Y2hat$z, surf.Y3hat$z))

# Size for the mapping
width <- 360*3
height <- 360*2
pointsize <- 12

png("output/surface_M_TLmis.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(2, 3))

iy1 <- as.image.SpatialGridDataFrame(surf.Y1)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 1st response") 
axis(2, las=1)
axis(1)
image.plot(iy1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy2 <- as.image.SpatialGridDataFrame(surf.Y2)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 2nd response") 
axis(2, las=1)
axis(1)
image.plot(iy2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy3 <- as.image.SpatialGridDataFrame(surf.Y3)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 3rd response") 
axis(2, las=1)
axis(1)
image.plot(iy3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iyp1 <- as.image.SpatialGridDataFrame(surf.Y1hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 1st response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[1], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp2 <- as.image.SpatialGridDataFrame(surf.Y2hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 2nd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[2], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp3 <- as.image.SpatialGridDataFrame(surf.Y3hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 3rd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[3], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

dev.off()

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("output/CIpost_M_TLmis.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(postint_ms)
dev.off()



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
pb <- txtProgressBar(0, K, style = 3)
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
  
  setTxtProgressBar(pb, i)
  
}
fit_time_hms <- toc()

# compute MAP estimates
betahat <- beta_smp[[K]] |> apply(c(1,2), mean)
sigma2hat <- sigma2_smp[[K]] |> apply(c(1,2), mean)
What <- do.call(rbind, lapply(1:K, function(k) W_smp[[k]] |> apply(c(1,2), mean)))
ordW <- order(do.call(rbind, data_part$sets))
What <- What[ordW,]

# plot posterior inference results
true_par <- c(sigma2[upper.tri(sigma2,T)], matrix(B))
sigma2_mat <- matrix(sigma2_smp[[K]], nrow = R, byrow = T)
beta_mat <- matrix(beta_smp[[K]], nrow = R, byrow = T)
posterior_hms <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
colnames(posterior_hms)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                            "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
postint_hms <- mcmc_recover_intervals(posterior_hms, true_par,
                                     prob = 0.95,
                                     prob_outer = 0.95,
                                     point_est = "mean",
                                     size = 4,
                                     alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[list(1,1)]), expression(beta[list(2,1)]), expression(beta[list(1,2)]),
                              expression(beta[list(2,2)]), expression(beta[list(1,3)]), expression(beta[list(2,3)]),
                              expression(Sigma[list(1,1)]), expression(Sigma[list(2,1)]), expression(Sigma[list(1,3)]),
                              expression(Sigma[list(2,2)]), expression(Sigma[list(2,3)]), expression(Sigma[list(3,3)]))) + 
  ggtitle("Credible intervals - highly misspecified - Transfer Learning",
          "with posterior means, true values, and 95% credible intervals")


# predictions
YK <- data_part$Y_list[[K]]; XK <- data_part$X_list[[K]]; crd_K <- data_part$crd_list[[K]]
tic()
predictions_hms <- r_pred_marg_MvT(data = list(Y = YK, X = XK), X_u = X_u, R = 250,
                                  d_u = arma_dist(crd_u), d_us = arma_dist(rbind(crd_u, crd_K)),
                                  hyperpar = hpar, poster = out)$Yu
toc()

# statistics computations Y
pred_mat_Y <- predictions_hms
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_low_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025))
post_upp_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.975))

# Empirical coverage for Y
coverage_Y <- colMeans(Y_u >= post_low_Y & Y_u <= post_upp_Y)
cat("Empirical average coverage for Response:", round(mean(coverage_Y), 3))

# Root Mean Square Prediction Error
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) ))

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y1 <- MBA::mba.surf(cbind(crd_u, Y_u[,1]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y1hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,1]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2 <- MBA::mba.surf(cbind(crd_u, Y_u[,2]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,2]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3 <- MBA::mba.surf(cbind(crd_u, Y_u[,3]), no.X = 500, no.Y = 500, 
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,3]), no.X = 500, no.Y = 500, 
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y1$z, surf.Y1hat$z, surf.Y2$z, surf.Y2hat$z, surf.Y3$z, surf.Y3hat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y1$z, surf.Y2$z, surf.Y3$z))
zlimp <- range(c(surf.Y1hat$z, surf.Y2hat$z, surf.Y3hat$z))

# Size for the mapping
width <- 360*3
height <- 360*2
pointsize <- 12

png("output/surface_M_TLhms.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(2, 3))

iy1 <- as.image.SpatialGridDataFrame(surf.Y1)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 1st response") 
axis(2, las=1)
axis(1)
image.plot(iy1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy2 <- as.image.SpatialGridDataFrame(surf.Y2)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 2nd response") 
axis(2, las=1)
axis(1)
image.plot(iy2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy3 <- as.image.SpatialGridDataFrame(surf.Y3)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 3rd response") 
axis(2, las=1)
axis(1)
image.plot(iy3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iyp1 <- as.image.SpatialGridDataFrame(surf.Y1hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 1st response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[1], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp2 <- as.image.SpatialGridDataFrame(surf.Y2hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 2nd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[2], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp3 <- as.image.SpatialGridDataFrame(surf.Y3hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 3rd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[3], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

dev.off()

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("output/CIpost_M_TLhms.png", width = width, height = height, pointsize = pointsize, family = "sans")
# win.graph()
cowplot::plot_grid(postint_hms)
dev.off()


# BPS - Parallel Learning ------------------------------------------------

alfa_seq <- seq(0.75, 0.85, 0.05)
phi_seq <- seq(8, 12, 2)

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

# expected values of hyperparamters
sapply(1:K, function(k) t(W_list[[k]]) %*% expand_grid_cpp(alfa_seq, phi_seq)) %*% Wbps 

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]]$Pred[[1]]$Yu}, simplify = "array")
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(Y_u[,1] >= post_qnt_Y[1,,1] & Y_u[,1] <= post_qnt_Y[2,,1]),
                mean(Y_u[,2] >= post_qnt_Y[1,,2] & Y_u[,2] <= post_qnt_Y[2,,2]))
cat("Empirical average coverage for Spatial process:", round(mean(coverage_Y), 3))
(CI_avlen_mbps <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
# (rmspe_W <- sqrt( colMeans( (W_u[-(1:m),] - post_mean_W[-(1:m),])^2 ) )); mean(rmspe_W)
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) )); mean(rmspe_Y)

# plot the spatial process surfaces true vs. recovered
h <- 12
surf.Y1 <- MBA::mba.surf(cbind(crd_u, Y_u[,1]), no.X = 500, no.Y = 500, 
                         exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y1hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,1]), no.X = 500, no.Y = 500, 
                            exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2 <- MBA::mba.surf(cbind(crd_u, Y_u[,2]), no.X = 500, no.Y = 500, 
                         exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y2hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,2]), no.X = 500, no.Y = 500, 
                            exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3 <- MBA::mba.surf(cbind(crd_u, Y_u[,3]), no.X = 500, no.Y = 500, 
                         exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.Y3hat <- MBA::mba.surf(cbind(crd_u, post_mean_Y[,3]), no.X = 500, no.Y = 500, 
                            exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(c(surf.Y1$z, surf.Y1hat$z, surf.Y2$z, surf.Y2hat$z, surf.Y3$z, surf.Y3hat$z), 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)

zlim <- range(c(surf.Y1$z, surf.Y2$z, surf.Y3$z))
zlimp <- range(c(surf.Y1hat$z, surf.Y2hat$z, surf.Y3hat$z))

# Size for the mapping
width <- 360*3
height <- 360*2
pointsize <- 12

png("output/CIpost_M_TLbps.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(2, 3))

iy1 <- as.image.SpatialGridDataFrame(surf.Y1)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 1st response") 
axis(2, las=1)
axis(1)
image.plot(iy1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy2 <- as.image.SpatialGridDataFrame(surf.Y2)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 2nd response") 
axis(2, las=1)
axis(1)
image.plot(iy2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iy3 <- as.image.SpatialGridDataFrame(surf.Y3)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="True 3rd response") 
axis(2, las=1)
axis(1)
image.plot(iy3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

iyp1 <- as.image.SpatialGridDataFrame(surf.Y1hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 1st response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[1], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp2 <- as.image.SpatialGridDataFrame(surf.Y2hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 2nd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[2], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

iyp3 <- as.image.SpatialGridDataFrame(surf.Y3hat)
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting")
title(main="Estimated 3rd response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y[3], 3)))
axis(2, las=1)
axis(1)
image.plot(iyp3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

dev.off()

# plot posterior inference results
true_par <- c(sigma2[upper.tri(sigma2,T)], matrix(B))
sigma2_mat <- matrix(sapply(1:R, function(r){predictions[[r]]$Post[[1]]$sigma}, simplify = "array"), nrow = R, byrow = T)
beta_mat <- matrix(sapply(1:R, function(r){predictions[[r]]$Post[[1]]$beta[1:p,]}, simplify = "array"), nrow = R, byrow = T)
posterior_bps <- cbind(sigma2_mat[,c(1,2,5,3,6,9)], beta_mat)
colnames(posterior_bps)<- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[2,2]", "Sigma[1,3]", "Sigma[2,3]", "Sigma[3,3]",
                           "beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "beta[1,3]", "beta[2,3]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
postint_bps <- mcmc_recover_intervals(posterior_bps, true_par,
                                     prob = 0.95,
                                     prob_outer = 0.95,
                                     point_est = "mean",
                                     size = 4,
                                     alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[list(1,1)]), expression(beta[list(2,1)]), expression(beta[list(1,2)]),
                              expression(beta[list(2,2)]), expression(beta[list(1,3)]), expression(beta[list(2,3)]),
                              expression(Sigma[list(1,1)]), expression(Sigma[list(2,1)]), expression(Sigma[list(1,3)]),
                              expression(Sigma[list(2,2)]), expression(Sigma[list(2,3)]), expression(Sigma[list(3,3)]))) + 
  ggtitle("Credible intervals - BPS - Parallel Learning",
          "with posterior means, true values, and 95% credible intervals")

# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 12

png("output/CIpost_M_TLbps.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(postint_bps)
dev.off()
