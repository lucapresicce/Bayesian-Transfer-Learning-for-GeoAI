#####################################################################################################################################################
## BPS vs SMK COMPARISON - 5K - 500 ############################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets")

# Packages --------------------------------------------------------------------
library(spBPS)
library(Rcpp)
library(RcppArmadillo)
library(mniw)
library(MCMCpack)
library(ggplot2)
library(tictoc)
# devtools::install_github("david-dunson/Mposterior")
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
p <- 4

# parameters
B <- c(-0.75, 0.90, -1.1, 1.85)
tau2 <- 0.25
sigma2 <- 1
delta <- tau2/sigma2
phi <- 4

set.seed(97)
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
Y <- Y_or[1:n, ]

# prediction data
crd_u <- crd[-(1:(n-m)), ]
X_u <- X_or[-(1:(n-m)), ]
W_u <- W_or[-(1:(n-m)), ]
Y_u <- Y_or[-(1:(n-m)), ]

# Subset posterior models -----------------------------------------------------

# hyperparameters values
delta_seq <- c(0.2, 0.25, 0.3)
phi_seq <- c(3, 4, 5)

# function for the fit loop
fit_loop <- function(i) {

  Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
  p <- ncol(Xi)
  bps <- spBPS::BPS_weights(data = list(Y = Yi, X = Xi),
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
  result <- spBPS::BPS_post(data = list(Y = Ys, X = Xs), coords = crds,
                           X_u = X_u, crd_u = crd_u,
                           priors = list(mu_b = matrix(rep(0, p)),
                                         V_b = diag(10, p),
                                         a = 2,
                                         b = 2),
                           hyperpar = list(delta = delta_seq, phi = phi_seq),
                           W = Ws, R = 1)

  return(result)
}


# subsetting data
subset_size <- 500
K <- n/subset_size
data_part <- subset_data(data = list(Y = matrix(Y), X = X, crd = crd_s), K = K)

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

# statistics computations W
pred_mat_W <- sapply(1:R, function(r){predictions[[r]][[1]]})
post_mean_W <- rowMeans(pred_mat_W)
post_var_W <- apply(pred_mat_W, 1, sd)
post_qnt_W <- apply(pred_mat_W, 1, quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- mean(W_u >= post_qnt_W[1,] & W_u <= post_qnt_W[2,])
cat("Empirical coverage for Spatial process:", round(coverage_W, 3))

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]][[2]]})
post_mean_Y <- rowMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,]-post_qnt_Y[1,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                   "Combination" = as.numeric(comb_time$toc-comb_time$tic),
                   "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                   "Total time" = as.numeric(tot_time$toc-tot_time$tic))


# SMK from Guhaniyogi (2018) --------------------------------------------------
gc()
## Libraries
library(mcmc)
library(MASS)
library(KernSmooth)
library(fields)
library(pscl)
library(spBayes)
library(mvtnorm)
library(MCMCpack)
library(Mposterior)
library(parallel)
library(doParallel)
library(foreach)

# load common simulation data
dat.nomiss <- cbind(crd_s, Y, X)
dat.miss <- cbind(crd_u, Y_u, X_u)

N.core <- n.core
n.core <- K

## Useful Quantities
## Total number of training samples
n.sample <- nrow(dat.nomiss)
## Total number of test samples
n.test <- nrow(dat.miss)
## Number of cores or subsets
# n.core <- K
per.core <- floor(5000/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part)
## Number of MCMC iterations
mcmc.sample <- 2500
## Burn in
n.burn <- 0.5*mcmc.sample
## Divide predicted data and predict them
## independently in each subset
# p.sub <- c(0, seq(per.core, n.sample, by = per.core))
per.core2 <- floor(m/n.core)
p.sub <- c(0, seq(per.core2, per.core2*(n.core-1), by = per.core2), m)


## GP regression on data subsets

a <- 1:n.sample
sample.loc <- dat.nomiss[,1:2]
y <- dat.nomiss[,3]
x <- dat.nomiss[,-(1:3)]
## Training response
Y.train <- y
X.train <- x
## Training coordinates
coords.train <- sample.loc
## Test coordinates
coords.pred <- dat.miss[,1:2]
## Test sample size
n.test <- nrow(dat.miss)
index.part <- list()
X.part <- list()
Y.part <- list()
coords.part <- list()

# sample the index of subsets by removing each time the sampled set with setdiff()
# it removes the second argument to the first argument
for(i in 1:n.split){
  beg<-Sys.time()
  index.part[[i]] <- sample(a,n.part[i],replace=FALSE)
  ## Response in i th subset
  Y.part[[i]] <- Y.train[index.part[[i]]]
  ## Predictor in i th subset
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset
  coords.part[[i]] <- coords.train[index.part[[i]],]
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){
  ## Model fitting

  library(spBayes)
  starting <- list("phi"=3, "sigma.sq"=5, "tau.sq"=1)
  tuning <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01)
  priors.1 <- list("beta.Norm"=list(c(0,0,0,0), diag(1000, p)),
                   "phi.Unif"=c(3/10, 3/0.1), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 0.1))
  cov.model <- "exponential"
  ## Response in subset i
  ZZ <- Y.part[[i]]
  ## Predictor in subset i
  XX <- X.part[[i]]
  ## Coordinates in subset i
  CC <- coords.part[[i]]

  ## GP computation in each subset
  m.1 <- spLM(ZZ~XX-1, coords=CC, starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=mcmc.sample,verbose=FALSE)
  ## Recover all MCMC samples in each subset
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=FALSE)
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning
  gc()

  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()

  # tic()
  ## Prediction
  Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
                         pred.coords=coords.pred[-(1:m),],
                         start=0.5*mcmc.sample)$p.y.predictive.samples)
  Y.pred <- NULL
  for(l in 1:(length(p.sub)-1)){
    m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
                          pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
                          start=0.5*mcmc.sample)$p.y.predictive.samples
    Y.pred <- cbind(Y.pred,t(m.1.pred))
  }

  Y.pred <- cbind(Y.pred, Yu.pred)

  # toc()

  ## MCMC samples for both parameter estimates and prediction after burn-in
  hh <- list(subAtom,Y.pred,W.mcmc)
  names(hh) <- c("atoms","predictions","W.mcmc")
  return(hh)
}

####### Parallelization
## Number of clusters for parallel implementation
# cl<-makeCluster(n.core)
cl<-makeCluster(N.core)
registerDoParallel(cl)

## Start time
strt<-Sys.time()
## Parallelized subset computation of GP in different cores
# obj <- foreach(i=1:n.core) %dopar% partitioned_GP(i)
obj <- foreach(i=1:K) %dopar% partitioned_GP(i)
## Total time for parallelized inference
final.time <- Sys.time()-strt
stopCluster(cl)

####### Combine

subAtomList <- list()
for(i in 1:length(n.part)){
  ## MCMC samples to run Weiszfeld algorithm
  subAtomList[[i]] <- obj[[i]]$atoms[seq(1,mcmc.sample-n.burn,10),]
  attr(subAtomList[[i]], "class") <- "matrix"
}

beg.rec <- Sys.time()
## Combination using Weiszfeld algorithm
medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-5)
## Posterior wts. using Weiszfeld's algorithm
wts <- medPosterior$weiszfeldWts
## Time for combining subset posteriors
recomb.time <- Sys.time()-beg.rec

####### Prediction
Y.med <- numeric()
Y.lower_qnt <- numeric()
Y.upper_qnt <- numeric()
beg <- Sys.time()
## Number of posterior predictive samples
n.pred.sample <- nrow(obj[[1]]$predictions)

for(j in 1:ncol(obj[[1]]$predictions)){
  ## Posterior predictive samples for j th predicted point from first subsample
  a.sub.i<- obj[[1]]$predictions[,j]
  a.wts <- wts[1]*rep(1/n.pred.sample,n.pred.sample)
  for(i in 2:n.split){
    ## Posterior predictive samples for j th predicted point from i th subsample
    a.sub.i <- c(a.sub.i,obj[[i]]$predictions[,j])
    ## Corresponding empirical weights
    a.wts <- c(a.wts, wts[i]*rep(1/n.pred.sample,n.pred.sample))
  }
  new.atom <- sort(a.sub.i,decreasing=F,index.return=T)$ix
  ## MCMC samples for the meta posterior
  atoms <- a.sub.i[new.atom]
  ## Corresponding weights
  prob <- a.wts[new.atom]
  id11 <- min(which(cumsum(prob)>=0.025))
  id12 <- min(which(cumsum(prob)>=0.975))
  id13 <- min(which(cumsum(prob)>=0.5))
  ## Posterior 2.5% quantile for j th predicted point
  Y.lower_qnt[j] <- atoms[id11]
  ## Posterior 97.5% quantile for j th predicted point
  Y.upper_qnt[j] <- atoms[id12]
  ## Posterior median for j th predicted point
  Y.med[j] <- atoms[id13]
}

## Time for calculating quantiles
quantile.calculation.time <- Sys.time()-beg

## Total time for subset running, combining subset posteriors and calculating quantiles
Time <- final.time+recomb.time+quantile.calculation.time

# Posterior inference -----------------------------------------------------

# bps
smp <- sapply(1:R, function(r){c(predictions[[r]][[4]], predictions[[r]][[3]][,1:p])})
post_mean_smp <- rowMeans(smp)
post_var_smp <- apply(smp, 1, sd)
post_qnt_smp <- apply(smp, 1, quantile, c(0.05, 0.95))
post_mean_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps
post_var_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# smk
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:K, function(a)obj[[a]]$atoms[,j]) %*% wts
}
post_smp_smk |> colMeans()

# posterior inference - credible interval
true_par <- c(sigma2, B)

# posterior samples collection
posterior_bps <- t(smp)
posterior_smk <- cbind(post_smp_smk[,5], post_smp_smk[,1:4])
colnames(posterior_bps) <- colnames(posterior_smk) <- c("sigma^2", "beta[1]", "beta[2]", "beta[3]", "beta[4]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_bps <- mcmc_recover_intervals(posterior_bps, true_par,
                                        prob = 0.95,
                                        prob_outer = 0.95,
                                        point_est = "mean",
                                        size = 4,
                                        alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals BPS ",
          "with posterior means, true values, and 95% credible intervals")


# plotting smk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_smk <- mcmc_recover_intervals(posterior_smk, true_par,
                                       prob = 0.95,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 3,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals SMK",
          "with posterior means, true values, and 95% credible intervals")

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12
png("output/CIpost_5_500.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(post_int_bps, post_int_smk)
dev.off()


# Compare with SMK ------------------------------------------------------------

# load the SMK predictions
(rmspe_smk <- sqrt( mean( (Y_u[-(1:m)] - Y.med[-(1:m)])^2 ) ))
# (rmspe_smk <- sqrt( mean( (Y_u - Y.med)^2 ) ))

coverage_smk <- mean(Y_u >= Y.lower_qnt & Y_u <= Y.upper_qnt)
cat("Empirical coverage SMK :", round(coverage_smk, 3))
(CI_avlen_smk <- mean(Y.upper_qnt-Y.lower_qnt))

elapsed_times_smk <- c("Fitting" = final.time,
                       "Combination" = recomb.time,
                       "Prediction" = quantile.calculation.time,
                       "Total time" = Time)

cat("BPS  timing : \n"); round(elapsed_times/60, 2)
cat("SMK timing : \n"); round(elapsed_times_smk/60, 2)

# Plotting data --------------------------------------------------
gc()

# true surfaces interpolation
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_s, Y), no.X = 500, no.Y = 500,
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(surf.Y$z, 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)
zlim <- range(surf.Y$z)

# image for plot
iy <- as.image.SpatialGridDataFrame(surf.Y)

# BPS surfaces interpolation
h <- 12
surf.Yp <- MBA::mba.surf(cbind(crd_u, post_mean_Y), no.X = 500, no.Y = 500,
                         exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimp <- range(surf.Yp$z)

# image for plot
iyp <- as.image.SpatialGridDataFrame(surf.Yp)

# SMK surfaces interpolation
h <- 12
surf.Ysmk <- MBA::mba.surf(cbind(crd_u, Y.med), no.X = 500, no.Y = 500,
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimsmk <- range(surf.Ysmk$z)

# image for plotting
iysmk <- as.image.SpatialGridDataFrame(surf.Ysmk)

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

# Plotting
png("output/surface_5_500.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 3))
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="Response")
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="BPS  Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iyp, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="SMK Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_smk, 3)))
axis(2, las=1)
axis(1)
image.plot(iysmk, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimsmk)
dev.off()

# graphical UC comparison for Y (BPS  vs SMK - ordered) 
ord_y <- order(Y_u[-(1:m)])
df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Y_lower = post_qnt_Y[,-(1:m)][, ord_y][1,],
  CI_Y_upper = post_qnt_Y[,-(1:m)][, ord_y][2,],
  Ymap_ord = post_mean_Y[-(1:m)][ord_y])
# Create the ggplot
uc_Ybps <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Y_lower, ymax = CI_Y_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Y_lower, df$CI_Y_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "BPS  Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_Y, 3))) +
  geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) 

df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Ysmk_lower = Y.lower_qnt[-(1:m)][ord_y],
  CI_Ysmk_upper = Y.upper_qnt[-(1:m)][ord_y],
  Ysmk_ord = Y.med[-(1:m)][ord_y])
# Create the ggplot
uc_Ysmk <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Ysmk_lower, ymax = CI_Ysmk_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Ysmk_lower, df$CI_Ysmk_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "SMK Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_smk, 3))) +
  geom_point(aes(y = Ysmk_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

png("output/UC_5_500.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(uc_Ybps, uc_Ysmk)
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS " = list("time" = elapsed_times,
                              "RMSPE" = rmspe_Y),
                "SMK" = list("time" = elapsed_times_smk,
                             "RMSPE" = rmspe_smk))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "output/simulation_univariate_5_500.RData")

#####################################################################################################################################################
## BPS vs SMK COMPARISON - 5K - 1000 ############################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets")
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
p <- 4

# parameters
B <- c(-0.75, 0.90, -1.1, 1.85)
tau2 <- 0.25
sigma2 <- 1
delta <- tau2/sigma2
phi <- 4

set.seed(97)
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
Y <- Y_or[1:n, ]

# prediction data
crd_u <- crd[-(1:(n-m)), ]
X_u <- X_or[-(1:(n-m)), ]
W_u <- W_or[-(1:(n-m)), ]
Y_u <- Y_or[-(1:(n-m)), ]

# Subset posterior models -----------------------------------------------------

# hyperparameters values
delta_seq <- c(0.2, 0.25, 0.3)
phi_seq <- c(3, 4, 5)

# function for the fit loop
fit_loop <- function(i) {

  Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
  p <- ncol(Xi)
  bps <- spBPS::BPS_weights(data = list(Y = Yi, X = Xi),
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
  result <- spBPS::BPS_post(data = list(Y = Ys, X = Xs), coords = crds,
                           X_u = X_u, crd_u = crd_u,
                           priors = list(mu_b = matrix(rep(0, p)),
                                         V_b = diag(10, p),
                                         a = 2,
                                         b = 2),
                           hyperpar = list(delta = delta_seq, phi = phi_seq),
                           W = Ws, R = 1)

  return(result)
}

# subsetting data
subset_size <- 1000
K <- n/subset_size
data_part <- subset_data(data = list(Y = matrix(Y), X = X, crd = crd_s), K = K)

# BPS  parallel fit -------------------------------------------------------

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

# statistics computations W
pred_mat_W <- sapply(1:R, function(r){predictions[[r]][[1]]})
post_mean_W <- rowMeans(pred_mat_W)
post_var_W <- apply(pred_mat_W, 1, sd)
post_qnt_W <- apply(pred_mat_W, 1, quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- mean(W_u >= post_qnt_W[1,] & W_u <= post_qnt_W[2,])
cat("Empirical coverage for Spatial process:", round(coverage_W, 3))

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]][[2]]})
post_mean_Y <- rowMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,]-post_qnt_Y[1,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                   "Combination" = as.numeric(comb_time$toc-comb_time$tic),
                   "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                   # "Post sample" = as.numeric(post_time$toc-post_time$tic),
                   "Total time" = as.numeric(tot_time$toc-tot_time$tic))


# SMK from Guhaniyogi (2018) --------------------------------------------------
gc()
## Libraries
library(mcmc)
library(MASS)
library(KernSmooth)
library(fields)
library(pscl)
library(spBayes)
library(mvtnorm)
library(MCMCpack)
library(Mposterior)
library(parallel)
library(doParallel)
library(foreach)

# load common simulation data
dat.nomiss <- cbind(crd_s, Y, X)
dat.miss <- cbind(crd_u, Y_u, X_u)

N.core <- n.core
n.core <- K

## Useful Quantities
## Total number of training samples
n.sample <- nrow(dat.nomiss)
## Total number of test samples
n.test <- nrow(dat.miss)
## Number of cores or subsets
# n.core <- K
per.core <- floor(5000/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part)
## Number of MCMC iterations
mcmc.sample <- 2500
## Burn in
n.burn <- 0.5*mcmc.sample
## Divide predicted data and predict them
## independently in each subset
# p.sub <- c(0, seq(per.core, n.sample, by = per.core))
per.core2 <- floor(m/n.core)
p.sub <- c(0, seq(per.core2, per.core2*(n.core-1), by = per.core2), m)


## GP regression on data subsets

a <- 1:n.sample
sample.loc <- dat.nomiss[,1:2]
y <- dat.nomiss[,3]
x <- dat.nomiss[,-(1:3)]
## Training response
Y.train <- y
X.train <- x
## Training coordinates
coords.train <- sample.loc
## Test coordinates
coords.pred <- dat.miss[,1:2]
## Test sample size
n.test <- nrow(dat.miss)
index.part <- list()
X.part <- list()
Y.part <- list()
coords.part <- list()

# sample the index of subsets by removing each time the sampled set with setdiff()
# it removes the second argument to the first argument
for(i in 1:n.split){
  beg<-Sys.time()
  index.part[[i]] <- sample(a,n.part[i],replace=FALSE)
  ## Response in i th subset
  Y.part[[i]] <- Y.train[index.part[[i]]]
  ## Predictor in i th subset
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset
  coords.part[[i]] <- coords.train[index.part[[i]],]
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){
  ## Model fitting

  library(spBayes)
  starting <- list("phi"=3, "sigma.sq"=5, "tau.sq"=1)
  tuning <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01)
  priors.1 <- list("beta.Norm"=list(c(0,0,0,0), diag(1000, p)),
                   "phi.Unif"=c(3/10, 3/0.1), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 0.1))
  cov.model <- "exponential"
  ## Response in subset i
  ZZ <- Y.part[[i]]
  ## Predictor in subset i
  XX <- X.part[[i]]
  ## Coordinates in subset i
  CC <- coords.part[[i]]

  ## GP computation in each subset
  m.1 <- spLM(ZZ~XX-1, coords=CC, starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=mcmc.sample,verbose=FALSE)
  ## Recover all MCMC samples in each subset
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=FALSE)
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning
  gc()

  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()

  # tic()
  ## Prediction
  Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
                         pred.coords=coords.pred[-(1:m),],
                         start=0.5*mcmc.sample)$p.y.predictive.samples)
  Y.pred <- NULL
  for(l in 1:(length(p.sub)-1)){
    m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
                          pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
                          start=0.5*mcmc.sample)$p.y.predictive.samples
    Y.pred <- cbind(Y.pred,t(m.1.pred))
  }

  Y.pred <- cbind(Y.pred, Yu.pred)

  # toc()

  ## MCMC samples for both parameter estimates and prediction after burn-in
  hh <- list(subAtom,Y.pred,W.mcmc)
  names(hh) <- c("atoms","predictions","W.mcmc")
  return(hh)
}

####### Parallelization
## Number of clusters for parallel implementation
# cl<-makeCluster(n.core)
cl<-makeCluster(N.core)
registerDoParallel(cl)

## Start time
strt<-Sys.time()
## Parallelized subset computation of GP in different cores
# obj <- foreach(i=1:n.core) %dopar% partitioned_GP(i)
obj <- foreach(i=1:K) %dopar% partitioned_GP(i)
## Total time for parallelized inference
final.time <- Sys.time()-strt
stopCluster(cl)

####### Combine

subAtomList <- list()
for(i in 1:length(n.part)){
  ## MCMC samples to run Weiszfeld algorithm
  subAtomList[[i]] <- obj[[i]]$atoms[seq(1,mcmc.sample-n.burn,10),]
  attr(subAtomList[[i]], "class") <- "matrix"
}

beg.rec <- Sys.time()
## Combination using Weiszfeld algorithm
medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-5)
## Posterior wts. using Weiszfeld's algorithm
wts <- medPosterior$weiszfeldWts
## Time for combining subset posteriors
recomb.time <- Sys.time()-beg.rec

####### Prediction
Y.med <- numeric()
Y.lower_qnt <- numeric()
Y.upper_qnt <- numeric()
beg <- Sys.time()
## Number of posterior predictive samples
n.pred.sample <- nrow(obj[[1]]$predictions)

for(j in 1:ncol(obj[[1]]$predictions)){
  ## Posterior predictive samples for j th predicted point from first subsample
  a.sub.i<- obj[[1]]$predictions[,j]
  a.wts <- wts[1]*rep(1/n.pred.sample,n.pred.sample)
  for(i in 2:n.split){
    ## Posterior predictive samples for j th predicted point from i th subsample
    a.sub.i <- c(a.sub.i,obj[[i]]$predictions[,j])
    ## Corresponding empirical weights
    a.wts <- c(a.wts, wts[i]*rep(1/n.pred.sample,n.pred.sample))
  }
  new.atom <- sort(a.sub.i,decreasing=F,index.return=T)$ix
  ## MCMC samples for the meta posterior
  atoms <- a.sub.i[new.atom]
  ## Corresponding weights
  prob <- a.wts[new.atom]
  id11 <- min(which(cumsum(prob)>=0.025))
  id12 <- min(which(cumsum(prob)>=0.975))
  id13 <- min(which(cumsum(prob)>=0.5))
  ## Posterior 2.5% quantile for j th predicted point
  Y.lower_qnt[j] <- atoms[id11]
  ## Posterior 97.5% quantile for j th predicted point
  Y.upper_qnt[j] <- atoms[id12]
  ## Posterior median for j th predicted point
  Y.med[j] <- atoms[id13]
}

## Time for calculating quantiles
quantile.calculation.time <- Sys.time()-beg

## Total time for subset running, combining subset posteriors and calculating quantiles
Time <- final.time+recomb.time+quantile.calculation.time

# Posterior inference -----------------------------------------------------

# bps
smp <- sapply(1:R, function(r){c(predictions[[r]][[4]], predictions[[r]][[3]][,1:p])})
post_mean_smp <- rowMeans(smp)
post_var_smp <- apply(smp, 1, sd)
post_qnt_smp <- apply(smp, 1, quantile, c(0.05, 0.95))
post_mean_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps
post_var_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# smk
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:K, function(a)obj[[a]]$atoms[,j]) %*% wts
}
post_smp_smk |> colMeans()

# posterior inference - credible interval
true_par <- c(sigma2, B)

# posterior samples collection
posterior_bps <- t(smp)
posterior_smk <- cbind(post_smp_smk[,5], post_smp_smk[,1:4])
colnames(posterior_bps) <- colnames(posterior_smk) <- c("sigma^2", "beta[1]", "beta[2]", "beta[3]", "beta[4]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_bps <- mcmc_recover_intervals(posterior_bps, true_par,
                                       prob = 0.95,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 4,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals BPS ",
          "with posterior means, true values, and 95% credible intervals")


# plotting smk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_smk <- mcmc_recover_intervals(posterior_smk, true_par,
                                       prob = 0.95,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 3,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals SMK",
          "with posterior means, true values, and 95% credible intervals")

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12
png("output/CIpost_5_1000.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(post_int_bps, post_int_smk)
dev.off()


# Compare with SMK ------------------------------------------------------------

# load the SMK predictions
(rmspe_smk <- sqrt( mean( (Y_u[-(1:m)] - Y.med[-(1:m)])^2 ) ))
# (rmspe_smk <- sqrt( mean( (Y_u - Y.med)^2 ) ))

coverage_smk <- mean(Y_u >= Y.lower_qnt & Y_u <= Y.upper_qnt)
cat("Empirical coverage SMK :", round(coverage_smk, 3))
(CI_avlen_smk <- mean(Y.upper_qnt-Y.lower_qnt))

elapsed_times_smk <- c("Fitting" = final.time,
                       "Combination" = recomb.time,
                       "Prediction" = quantile.calculation.time,
                       "Total time" = Time)

cat("BPS  timing : \n"); round(elapsed_times/60, 2)
cat("SMK timing : \n"); round(elapsed_times_smk/60, 2)

# Plotting data --------------------------------------------------
gc()

# true surfaces interpolation
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_s, Y), no.X = 500, no.Y = 500,
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(surf.Y$z, 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)
zlim <- range(surf.Y$z)

# image for plot
iy <- as.image.SpatialGridDataFrame(surf.Y)

# BPS surfaces interpolation
h <- 12
surf.Yp <- MBA::mba.surf(cbind(crd_u, post_mean_Y), no.X = 500, no.Y = 500,
                         exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimp <- range(surf.Yp$z)

# image for plot
iyp <- as.image.SpatialGridDataFrame(surf.Yp)

# SMK surfaces interpolation
h <- 12
surf.Ysmk <- MBA::mba.surf(cbind(crd_u, Y.med), no.X = 500, no.Y = 500,
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimsmk <- range(surf.Ysmk$z)

# image for plotting
iysmk <- as.image.SpatialGridDataFrame(surf.Ysmk)

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

# Plotting
png("output/surface_5_1000.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 3))
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="Response")
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="BPS  Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iyp, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="SMK Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_smk, 3)))
axis(2, las=1)
axis(1)
image.plot(iysmk, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimsmk)
dev.off()

# graphical UC comparison for Y (BPS  vs SMK - ordered) 
ord_y <- order(Y_u[-(1:m)])
df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Y_lower = post_qnt_Y[,-(1:m)][, ord_y][1,],
  CI_Y_upper = post_qnt_Y[,-(1:m)][, ord_y][2,],
  Ymap_ord = post_mean_Y[-(1:m)][ord_y])
# Create the ggplot
uc_Ybps <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Y_lower, ymax = CI_Y_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Y_lower, df$CI_Y_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "BPS  Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_Y, 3))) +
  geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) 

df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Ysmk_lower = Y.lower_qnt[-(1:m)][ord_y],
  CI_Ysmk_upper = Y.upper_qnt[-(1:m)][ord_y],
  Ysmk_ord = Y.med[-(1:m)][ord_y])
# Create the ggplot
uc_Ysmk <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Ysmk_lower, ymax = CI_Ysmk_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Ysmk_lower, df$CI_Ysmk_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "SMK Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_smk, 3))) +
  geom_point(aes(y = Ysmk_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

png("output/UC_5_1000.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(uc_Ybps, uc_Ysmk)
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS " = list("time" = elapsed_times,
                              "RMSPE" = rmspe_Y),
                "SMK" = list("time" = elapsed_times_smk,
                             "RMSPE" = rmspe_smk))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "output/simulation_univariate_5_1000.RData")

#####################################################################################################################################################
## BPS vs SMK COMPARISON - 10K - 500 ############################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets")
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
n <- 10000
m <- 500
u <- 250
p <- 4

# parameters
B <- c(-0.75, 0.90, -1.1, 1.85)
tau2 <- 0.25
sigma2 <- 1
delta <- tau2/sigma2
phi <- 4

set.seed(97)
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
Y <- Y_or[1:n, ]

# prediction data
crd_u <- crd[-(1:(n-m)), ]
X_u <- X_or[-(1:(n-m)), ]
W_u <- W_or[-(1:(n-m)), ]
Y_u <- Y_or[-(1:(n-m)), ]

# Subset posterior models -----------------------------------------------------

# hyperparameters values
delta_seq <- c(0.2, 0.25, 0.3)
phi_seq <- c(3, 4, 5)

# function for the fit loop
fit_loop <- function(i) {
  
  Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
  p <- ncol(Xi)
  bps <- spBPS::BPS_weights(data = list(Y = Yi, X = Xi),
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
  result <- spBPS::BPS_post(data = list(Y = Ys, X = Xs), coords = crds,
                           X_u = X_u, crd_u = crd_u,
                           priors = list(mu_b = matrix(rep(0, p)),
                                         V_b = diag(10, p),
                                         a = 2,
                                         b = 2),
                           hyperpar = list(delta = delta_seq, phi = phi_seq),
                           W = Ws, R = 1)
  
  return(result)
}

# subsetting data
subset_size <- 500
K <- n/subset_size
data_part <- subset_data(data = list(Y = matrix(Y), X = X, crd = crd_s), K = K)

# BPS  parallel fit -------------------------------------------------------

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

# statistics computations W
pred_mat_W <- sapply(1:R, function(r){predictions[[r]][[1]]})
post_mean_W <- rowMeans(pred_mat_W)
post_var_W <- apply(pred_mat_W, 1, sd)
post_qnt_W <- apply(pred_mat_W, 1, quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- mean(W_u >= post_qnt_W[1,] & W_u <= post_qnt_W[2,])
cat("Empirical coverage for Spatial process:", round(coverage_W, 3))

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]][[2]]})
post_mean_Y <- rowMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,]-post_qnt_Y[1,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                   "Combination" = as.numeric(comb_time$toc-comb_time$tic),
                   "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                   # "Post sample" = as.numeric(post_time$toc-post_time$tic),
                   "Total time" = as.numeric(tot_time$toc-tot_time$tic))


# SMK from Guhaniyogi (2018) --------------------------------------------------
gc()
## Libraries
library(mcmc)
library(MASS)
library(KernSmooth)
library(fields)
library(pscl)
library(spBayes)
library(mvtnorm)
library(MCMCpack)
library(Mposterior)
library(parallel)
library(doParallel)
library(foreach)

# load common simulation data
dat.nomiss <- cbind(crd_s, Y, X)
dat.miss <- cbind(crd_u, Y_u, X_u)

N.core <- n.core
n.core <- K

## Useful Quantities
## Total number of training samples
n.sample <- nrow(dat.nomiss)
## Total number of test samples 
n.test <- nrow(dat.miss)  
## Number of cores or subsets
# n.core <- K
per.core <- floor(5000/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part) 
## Number of MCMC iterations
mcmc.sample <- 2500
## Burn in 
n.burn <- 0.5*mcmc.sample  
## Divide predicted data and predict them
## independently in each subset
# p.sub <- c(0, seq(per.core, n.sample, by = per.core))
per.core2 <- floor(m/n.core)
p.sub <- c(0, seq(per.core2, per.core2*(n.core-1), by = per.core2), m)


## GP regression on data subsets

a <- 1:n.sample 
sample.loc <- dat.nomiss[,1:2] 
y <- dat.nomiss[,3]
x <- dat.nomiss[,-(1:3)]
## Training response
Y.train <- y  
X.train <- x
## Training coordinates
coords.train <- sample.loc
## Test coordinates  
coords.pred <- dat.miss[,1:2] 
## Test sample size 
n.test <- nrow(dat.miss) 
index.part <- list() 
X.part <- list()
Y.part <- list()
coords.part <- list()

# sample the index of subsets by removing each time the sampled set with setdiff()
# it removes the second argument to the first argument
for(i in 1:n.split){
  beg<-Sys.time()
  index.part[[i]] <- sample(a,n.part[i],replace=FALSE) 
  ## Response in i th subset 
  Y.part[[i]] <- Y.train[index.part[[i]]]  
  ## Predictor in i th subset 
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset       
  coords.part[[i]] <- coords.train[index.part[[i]],]  
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){
  ## Model fitting
  
  library(spBayes)
  starting <- list("phi"=3, "sigma.sq"=5, "tau.sq"=1) 
  tuning <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01)
  priors.1 <- list("beta.Norm"=list(c(0,0,0,0), diag(1000, p)),
                   "phi.Unif"=c(3/10, 3/0.1), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 0.1))
  cov.model <- "exponential"
  ## Response in subset i
  ZZ <- Y.part[[i]]   
  ## Predictor in subset i    
  XX <- X.part[[i]] 
  ## Coordinates in subset i      
  CC <- coords.part[[i]]  
  
  ## GP computation in each subset
  m.1 <- spLM(ZZ~XX-1, coords=CC, starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=mcmc.sample,verbose=FALSE)           
  ## Recover all MCMC samples in each subset
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=FALSE)  
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning  
  gc()          
  
  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()
  
  # tic()
  ## Prediction
  Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
                         pred.coords=coords.pred[-(1:m),],
                         start=0.5*mcmc.sample)$p.y.predictive.samples)   
  Y.pred <- NULL
  for(l in 1:(length(p.sub)-1)){
    m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
                          pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
                          start=0.5*mcmc.sample)$p.y.predictive.samples
    Y.pred <- cbind(Y.pred,t(m.1.pred))
  }
  
  Y.pred <- cbind(Y.pred, Yu.pred)
  
  # toc()
  
  ## MCMC samples for both parameter estimates and prediction after burn-in
  hh <- list(subAtom,Y.pred,W.mcmc)
  names(hh) <- c("atoms","predictions","W.mcmc") 
  return(hh)
}

####### Parallelization 
## Number of clusters for parallel implementation
# cl<-makeCluster(n.core)  
cl<-makeCluster(N.core)   
registerDoParallel(cl)

## Start time
strt<-Sys.time()
## Parallelized subset computation of GP in different cores
# obj <- foreach(i=1:n.core) %dopar% partitioned_GP(i)  
obj <- foreach(i=1:K) %dopar% partitioned_GP(i) 
## Total time for parallelized inference
final.time <- Sys.time()-strt  
stopCluster(cl)

####### Combine 

subAtomList <- list()
for(i in 1:length(n.part)){
  ## MCMC samples to run Weiszfeld algorithm
  subAtomList[[i]] <- obj[[i]]$atoms[seq(1,mcmc.sample-n.burn,10),]
  attr(subAtomList[[i]], "class") <- "matrix"
}

beg.rec <- Sys.time()
## Combination using Weiszfeld algorithm
medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-5) 
## Posterior wts. using Weiszfeld's algorithm
wts <- medPosterior$weiszfeldWts  
## Time for combining subset posteriors
recomb.time <- Sys.time()-beg.rec 

####### Prediction 
Y.med <- numeric()
Y.lower_qnt <- numeric()
Y.upper_qnt <- numeric()
beg <- Sys.time()
## Number of posterior predictive samples
n.pred.sample <- nrow(obj[[1]]$predictions) 

for(j in 1:ncol(obj[[1]]$predictions)){
  ## Posterior predictive samples for j th predicted point from first subsample
  a.sub.i<- obj[[1]]$predictions[,j]    
  a.wts <- wts[1]*rep(1/n.pred.sample,n.pred.sample)
  for(i in 2:n.split){
    ## Posterior predictive samples for j th predicted point from i th subsample
    a.sub.i <- c(a.sub.i,obj[[i]]$predictions[,j])  
    ## Corresponding empirical weights
    a.wts <- c(a.wts, wts[i]*rep(1/n.pred.sample,n.pred.sample)) 
  }
  new.atom <- sort(a.sub.i,decreasing=F,index.return=T)$ix
  ## MCMC samples for the meta posterior
  atoms <- a.sub.i[new.atom]
  ## Corresponding weights  
  prob <- a.wts[new.atom]     
  id11 <- min(which(cumsum(prob)>=0.025))
  id12 <- min(which(cumsum(prob)>=0.975))
  id13 <- min(which(cumsum(prob)>=0.5))
  ## Posterior 2.5% quantile for j th predicted point
  Y.lower_qnt[j] <- atoms[id11]     
  ## Posterior 97.5% quantile for j th predicted point      
  Y.upper_qnt[j] <- atoms[id12] 
  ## Posterior median for j th predicted point           
  Y.med[j] <- atoms[id13]                
}

## Time for calculating quantiles
quantile.calculation.time <- Sys.time()-beg  

## Total time for subset running, combining subset posteriors and calculating quantiles
Time <- final.time+recomb.time+quantile.calculation.time  

# Posterior inference -----------------------------------------------------

# bps
smp <- sapply(1:R, function(r){c(predictions[[r]][[4]], predictions[[r]][[3]][,1:p])})
post_mean_smp <- rowMeans(smp)
post_var_smp <- apply(smp, 1, sd)
post_qnt_smp <- apply(smp, 1, quantile, c(0.05, 0.95))
post_mean_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps
post_var_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# smk
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:K, function(a)obj[[a]]$atoms[,j]) %*% wts
}
post_smp_smk |> colMeans()

# posterior inference - credible interval
true_par <- c(sigma2, B)

# posterior samples collection
posterior_bps <- t(smp)
posterior_smk <- cbind(post_smp_smk[,5], post_smp_smk[,1:4])
colnames(posterior_bps) <- colnames(posterior_smk) <- c("sigma^2", "beta[1]", "beta[2]", "beta[3]", "beta[4]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_bps <- mcmc_recover_intervals(posterior_bps, true_par,
                                       prob = 0.95,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 4,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals BPS ",
          "with posterior means, true values, and 95% credible intervals")


# plotting smk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_smk <- mcmc_recover_intervals(posterior_smk, true_par,
                                       prob = 0.95,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 3,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals SMK",
          "with posterior means, true values, and 95% credible intervals")

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12
png("output/CIpost_10_500.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(post_int_bps, post_int_smk)
dev.off()


# Compare with SMK ------------------------------------------------------------

# load the SMK predictions
(rmspe_smk <- sqrt( mean( (Y_u[-(1:m)] - Y.med[-(1:m)])^2 ) ))
# (rmspe_smk <- sqrt( mean( (Y_u - Y.med)^2 ) ))

coverage_smk <- mean(Y_u >= Y.lower_qnt & Y_u <= Y.upper_qnt)
cat("Empirical coverage SMK :", round(coverage_smk, 3))
(CI_avlen_smk <- mean(Y.upper_qnt-Y.lower_qnt))

elapsed_times_smk <- c("Fitting" = final.time,
                       "Combination" = recomb.time,
                       "Prediction" = quantile.calculation.time,
                       "Total time" = Time)

cat("BPS  timing : \n"); round(elapsed_times/60, 2)
cat("SMK timing : \n"); round(elapsed_times_smk/60, 2)

# Plotting data --------------------------------------------------
gc()

# true surfaces interpolation
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_s, Y), no.X = 500, no.Y = 500,
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(surf.Y$z, 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)
zlim <- range(surf.Y$z)

# image for plot
iy <- as.image.SpatialGridDataFrame(surf.Y)

# BPS surfaces interpolation
h <- 12
surf.Yp <- MBA::mba.surf(cbind(crd_u, post_mean_Y), no.X = 500, no.Y = 500,
                         exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimp <- range(surf.Yp$z)

# image for plot
iyp <- as.image.SpatialGridDataFrame(surf.Yp)

# SMK surfaces interpolation
h <- 12
surf.Ysmk <- MBA::mba.surf(cbind(crd_u, Y.med), no.X = 500, no.Y = 500,
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimsmk <- range(surf.Ysmk$z)

# image for plotting
iysmk <- as.image.SpatialGridDataFrame(surf.Ysmk)

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

# Plotting
png("output/surface_10_500.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 3))
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="Response")
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="BPS  Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iyp, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="SMK Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_smk, 3)))
axis(2, las=1)
axis(1)
image.plot(iysmk, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimsmk)
dev.off()

# graphical UC comparison for Y (BPS  vs SMK - ordered) 
ord_y <- order(Y_u[-(1:m)])
df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Y_lower = post_qnt_Y[,-(1:m)][, ord_y][1,],
  CI_Y_upper = post_qnt_Y[,-(1:m)][, ord_y][2,],
  Ymap_ord = post_mean_Y[-(1:m)][ord_y])
# Create the ggplot
uc_Ybps <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Y_lower, ymax = CI_Y_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Y_lower, df$CI_Y_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "BPS  Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_Y, 3))) +
  geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) 

df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Ysmk_lower = Y.lower_qnt[-(1:m)][ord_y],
  CI_Ysmk_upper = Y.upper_qnt[-(1:m)][ord_y],
  Ysmk_ord = Y.med[-(1:m)][ord_y])
# Create the ggplot
uc_Ysmk <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Ysmk_lower, ymax = CI_Ysmk_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Ysmk_lower, df$CI_Ysmk_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "SMK Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_smk, 3))) +
  geom_point(aes(y = Ysmk_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

png("output/UC_10_500.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(uc_Ybps, uc_Ysmk)
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS " = list("time" = elapsed_times,
                              "RMSPE" = rmspe_Y),
                "SMK" = list("time" = elapsed_times_smk,
                             "RMSPE" = rmspe_smk))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "output/simulation_univariate_10_500.RData")

#####################################################################################################################################################
## BPS vs SMK COMPARISON - 10K - 1000 ############################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets")
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
n <- 10000
m <- 500
u <- 250
p <- 4

# parameters
B <- c(-0.75, 0.90, -1.1, 1.85)
tau2 <- 0.25
sigma2 <- 1
delta <- tau2/sigma2
phi <- 4

set.seed(97)
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
Y <- Y_or[1:n, ]

# prediction data
crd_u <- crd[-(1:(n-m)), ]
X_u <- X_or[-(1:(n-m)), ]
W_u <- W_or[-(1:(n-m)), ]
Y_u <- Y_or[-(1:(n-m)), ]

# Subset posterior models -----------------------------------------------------

# hyperparameters values
delta_seq <- c(0.2, 0.25, 0.3)
phi_seq <- c(3, 4, 5)

# function for the fit loop
fit_loop <- function(i) {
  
  Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
  p <- ncol(Xi)
  bps <- spBPS::BPS_weights(data = list(Y = Yi, X = Xi),
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
  result <- spBPS::BPS_post(data = list(Y = Ys, X = Xs), coords = crds,
                           X_u = X_u, crd_u = crd_u,
                           priors = list(mu_b = matrix(rep(0, p)),
                                         V_b = diag(10, p),
                                         a = 2,
                                         b = 2),
                           hyperpar = list(delta = delta_seq, phi = phi_seq),
                           W = Ws, R = 1)
  
  return(result)
}

# subsetting data
subset_size <- 1000
K <- n/subset_size
data_part <- subset_data(data = list(Y = matrix(Y), X = X, crd = crd_s), K = K)

# BPS  parallel fit -------------------------------------------------------

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

# statistics computations W
pred_mat_W <- sapply(1:R, function(r){predictions[[r]][[1]]})
post_mean_W <- rowMeans(pred_mat_W)
post_var_W <- apply(pred_mat_W, 1, sd)
post_qnt_W <- apply(pred_mat_W, 1, quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- mean(W_u >= post_qnt_W[1,] & W_u <= post_qnt_W[2,])
cat("Empirical coverage for Spatial process:", round(coverage_W, 3))

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]][[2]]})
post_mean_Y <- rowMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(Y_u >= post_qnt_Y[1,] & Y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,]-post_qnt_Y[1,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( mean( (W_u[-(1:m)] - post_mean_W[-(1:m)])^2 ) ))
(rmspe_Y <- sqrt( mean( (Y_u[-(1:m)] - post_mean_Y[-(1:m)])^2 ) ))

# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                   "Combination" = as.numeric(comb_time$toc-comb_time$tic),
                   "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                   # "Post sample" = as.numeric(post_time$toc-post_time$tic),
                   "Total time" = as.numeric(tot_time$toc-tot_time$tic))


# SMK from Guhaniyogi (2018) --------------------------------------------------
gc()
## Libraries
library(mcmc)
library(MASS)
library(KernSmooth)
library(fields)
library(pscl)
library(spBayes)
library(mvtnorm)
library(MCMCpack)
library(Mposterior)
library(parallel)
library(doParallel)
library(foreach)

# load common simulation data
dat.nomiss <- cbind(crd_s, Y, X)
dat.miss <- cbind(crd_u, Y_u, X_u)

N.core <- n.core
n.core <- K

## Useful Quantities
## Total number of training samples
n.sample <- nrow(dat.nomiss)
## Total number of test samples 
n.test <- nrow(dat.miss)  
## Number of cores or subsets
# n.core <- K
per.core <- floor(5000/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part) 
## Number of MCMC iterations
mcmc.sample <- 2500
## Burn in 
n.burn <- 0.5*mcmc.sample  
## Divide predicted data and predict them
## independently in each subset
# p.sub <- c(0, seq(per.core, n.sample, by = per.core))
per.core2 <- floor(m/n.core)
p.sub <- c(0, seq(per.core2, per.core2*(n.core-1), by = per.core2), m)


## GP regression on data subsets

a <- 1:n.sample 
sample.loc <- dat.nomiss[,1:2] 
y <- dat.nomiss[,3]
x <- dat.nomiss[,-(1:3)]
## Training response
Y.train <- y  
X.train <- x
## Training coordinates
coords.train <- sample.loc
## Test coordinates  
coords.pred <- dat.miss[,1:2] 
## Test sample size 
n.test <- nrow(dat.miss) 
index.part <- list() 
X.part <- list()
Y.part <- list()
coords.part <- list()

# sample the index of subsets by removing each time the sampled set with setdiff()
# it removes the second argument to the first argument
for(i in 1:n.split){
  beg<-Sys.time()
  index.part[[i]] <- sample(a,n.part[i],replace=FALSE) 
  ## Response in i th subset 
  Y.part[[i]] <- Y.train[index.part[[i]]]  
  ## Predictor in i th subset 
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset       
  coords.part[[i]] <- coords.train[index.part[[i]],]  
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){
  ## Model fitting
  
  library(spBayes)
  starting <- list("phi"=3, "sigma.sq"=5, "tau.sq"=1) 
  tuning <- list("phi"=0.01, "sigma.sq"=0.01, "tau.sq"=0.01)
  priors.1 <- list("beta.Norm"=list(c(0,0,0,0), diag(1000, p)),
                   "phi.Unif"=c(3/10, 3/0.1), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 0.1))
  cov.model <- "exponential"
  ## Response in subset i
  ZZ <- Y.part[[i]]   
  ## Predictor in subset i    
  XX <- X.part[[i]] 
  ## Coordinates in subset i      
  CC <- coords.part[[i]]  
  
  ## GP computation in each subset
  m.1 <- spLM(ZZ~XX-1, coords=CC, starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=mcmc.sample,verbose=FALSE)           
  ## Recover all MCMC samples in each subset
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=FALSE)  
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning  
  gc()          
  
  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()
  
  # tic()
  ## Prediction
  Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
                         pred.coords=coords.pred[-(1:m),],
                         start=0.5*mcmc.sample)$p.y.predictive.samples)   
  Y.pred <- NULL
  for(l in 1:(length(p.sub)-1)){
    m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
                          pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
                          start=0.5*mcmc.sample)$p.y.predictive.samples
    Y.pred <- cbind(Y.pred,t(m.1.pred))
  }
  
  Y.pred <- cbind(Y.pred, Yu.pred)
  
  # toc()
  
  ## MCMC samples for both parameter estimates and prediction after burn-in
  hh <- list(subAtom,Y.pred,W.mcmc)
  names(hh) <- c("atoms","predictions","W.mcmc") 
  return(hh)
}

####### Parallelization 
## Number of clusters for parallel implementation
# cl<-makeCluster(n.core)  
cl<-makeCluster(N.core)   
registerDoParallel(cl)

## Start time
strt<-Sys.time()
## Parallelized subset computation of GP in different cores
# obj <- foreach(i=1:n.core) %dopar% partitioned_GP(i)  
obj <- foreach(i=1:K) %dopar% partitioned_GP(i) 
## Total time for parallelized inference
final.time <- Sys.time()-strt  
stopCluster(cl)

####### Combine 

subAtomList <- list()
for(i in 1:length(n.part)){
  ## MCMC samples to run Weiszfeld algorithm
  subAtomList[[i]] <- obj[[i]]$atoms[seq(1,mcmc.sample-n.burn,10),]
  attr(subAtomList[[i]], "class") <- "matrix"
}

beg.rec <- Sys.time()
## Combination using Weiszfeld algorithm
medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-5) 
## Posterior wts. using Weiszfeld's algorithm
wts <- medPosterior$weiszfeldWts  
## Time for combining subset posteriors
recomb.time <- Sys.time()-beg.rec 

####### Prediction 
Y.med <- numeric()
Y.lower_qnt <- numeric()
Y.upper_qnt <- numeric()
beg <- Sys.time()
## Number of posterior predictive samples
n.pred.sample <- nrow(obj[[1]]$predictions) 

for(j in 1:ncol(obj[[1]]$predictions)){
  ## Posterior predictive samples for j th predicted point from first subsample
  a.sub.i<- obj[[1]]$predictions[,j]    
  a.wts <- wts[1]*rep(1/n.pred.sample,n.pred.sample)
  for(i in 2:n.split){
    ## Posterior predictive samples for j th predicted point from i th subsample
    a.sub.i <- c(a.sub.i,obj[[i]]$predictions[,j])  
    ## Corresponding empirical weights
    a.wts <- c(a.wts, wts[i]*rep(1/n.pred.sample,n.pred.sample)) 
  }
  new.atom <- sort(a.sub.i,decreasing=F,index.return=T)$ix
  ## MCMC samples for the meta posterior
  atoms <- a.sub.i[new.atom]
  ## Corresponding weights  
  prob <- a.wts[new.atom]     
  id11 <- min(which(cumsum(prob)>=0.025))
  id12 <- min(which(cumsum(prob)>=0.975))
  id13 <- min(which(cumsum(prob)>=0.5))
  ## Posterior 2.5% quantile for j th predicted point
  Y.lower_qnt[j] <- atoms[id11]     
  ## Posterior 97.5% quantile for j th predicted point      
  Y.upper_qnt[j] <- atoms[id12] 
  ## Posterior median for j th predicted point           
  Y.med[j] <- atoms[id13]                
}

## Time for calculating quantiles
quantile.calculation.time <- Sys.time()-beg  

## Total time for subset running, combining subset posteriors and calculating quantiles
Time <- final.time+recomb.time+quantile.calculation.time  

# Posterior inference -----------------------------------------------------

# bps
smp <- sapply(1:R, function(r){c(predictions[[r]][[4]], predictions[[r]][[3]][,1:p])})
post_mean_smp <- rowMeans(smp)
post_var_smp <- apply(smp, 1, sd)
post_qnt_smp <- apply(smp, 1, quantile, c(0.05, 0.95))
post_mean_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps
post_var_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# smk
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:K, function(a)obj[[a]]$atoms[,j]) %*% wts
}
post_smp_smk |> colMeans()

# posterior inference - credible interval
true_par <- c(sigma2, B)

# posterior samples collection
posterior_bps <- t(smp)
posterior_smk <- cbind(post_smp_smk[,5], post_smp_smk[,1:4])
colnames(posterior_bps) <- colnames(posterior_smk) <- c("sigma^2", "beta[1]", "beta[2]", "beta[3]", "beta[4]")

# plotting bps
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_bps <- mcmc_recover_intervals(posterior_bps, true_par,
                                       prob = 0.95,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 4,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals BPS ",
          "with posterior means, true values, and 95% credible intervals")


# plotting smk
bayesplot_theme_set(theme_default(base_size = 14, base_family = "sans"))
color_scheme_set("viridis")
post_int_smk <- mcmc_recover_intervals(posterior_smk, true_par,
                                       prob = 0.95,
                                       prob_outer = 0.95,
                                       point_est = "mean",
                                       size = 3,
                                       alpha = 0.75) +
  scale_x_discrete(labels = c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(beta[4]), expression(sigma^2))) + 
  ggtitle("", "") + 
  ggtitle("Credible intervals SMK",
          "with posterior means, true values, and 95% credible intervals")

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12
png("output/CIpost_10_1000.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(post_int_bps, post_int_smk)
dev.off()


# Compare with SMK ------------------------------------------------------------

# load the SMK predictions
(rmspe_smk <- sqrt( mean( (Y_u[-(1:m)] - Y.med[-(1:m)])^2 ) ))
# (rmspe_smk <- sqrt( mean( (Y_u - Y.med)^2 ) ))

coverage_smk <- mean(Y_u >= Y.lower_qnt & Y_u <= Y.upper_qnt)
cat("Empirical coverage SMK :", round(coverage_smk, 3))
(CI_avlen_smk <- mean(Y.upper_qnt-Y.lower_qnt))

elapsed_times_smk <- c("Fitting" = final.time,
                       "Combination" = recomb.time,
                       "Prediction" = quantile.calculation.time,
                       "Total time" = Time)

cat("BPS  timing : \n"); round(elapsed_times/60, 2)
cat("SMK timing : \n"); round(elapsed_times_smk/60, 2)

# Plotting data --------------------------------------------------
gc()

# true surfaces interpolation
h <- 12
surf.Y <- MBA::mba.surf(cbind(crd_s, Y), no.X = 500, no.Y = 500,
                        exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.brks <- classIntervals(surf.Y$z, 100, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1)
zlim <- range(surf.Y$z)

# image for plot
iy <- as.image.SpatialGridDataFrame(surf.Y)

# BPS surfaces interpolation
h <- 12
surf.Yp <- MBA::mba.surf(cbind(crd_u, post_mean_Y), no.X = 500, no.Y = 500,
                         exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimp <- range(surf.Yp$z)

# image for plot
iyp <- as.image.SpatialGridDataFrame(surf.Yp)

# SMK surfaces interpolation
h <- 12
surf.Ysmk <- MBA::mba.surf(cbind(crd_u, Y.med), no.X = 500, no.Y = 500,
                           exten = TRUE, sp = TRUE, h = h)$xyz.est
zlimsmk <- range(surf.Ysmk$z)

# image for plotting
iysmk <- as.image.SpatialGridDataFrame(surf.Ysmk)

# Size for the mapping
width <- 360*3
height <- 360
pointsize <- 12

# Plotting
png("output/surface_10_1000.png", width = width, height = height, pointsize = pointsize, family = "sans")
par(mfrow = c(1, 3))
plot(crd, type="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="Northing", xlab="Easting",
     main="Response")
axis(2, las=1)
axis(1)
image.plot(iy, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="BPS  Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_Y, 3)))
axis(2, las=1)
axis(1)
image.plot(iyp, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimp)

plot(crd, type="n", cex=0.5, xlim=xlim, axes=F, ylab="Northing", xlab="Easting")
title(main="SMK Prediction (Posterior Mean) for Response")
mtext(side = 3, paste("RMSPE :", round(rmspe_smk, 3)))
axis(2, las=1)
axis(1)
image.plot(iysmk, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlimsmk)
dev.off()

# graphical UC comparison for Y (BPS  vs SMK - ordered) 
ord_y <- order(Y_u[-(1:m)])
df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Y_lower = post_qnt_Y[,-(1:m)][, ord_y][1,],
  CI_Y_upper = post_qnt_Y[,-(1:m)][, ord_y][2,],
  Ymap_ord = post_mean_Y[-(1:m)][ord_y])
# Create the ggplot
uc_Ybps <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Y_lower, ymax = CI_Y_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Y_lower, df$CI_Y_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "BPS  Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_Y, 3))) +
  geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) 

df <- data.frame(
  x_ax = 1:u,
  Yu_ord = Y_u[-(1:m)][ord_y],
  CI_Ysmk_lower = Y.lower_qnt[-(1:m)][ord_y],
  CI_Ysmk_upper = Y.upper_qnt[-(1:m)][ord_y],
  Ysmk_ord = Y.med[-(1:m)][ord_y])
# Create the ggplot
uc_Ysmk <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
  geom_point(pch = 18, size = 2.5, col = "#1A85FF") +
  geom_errorbar(aes(ymin = CI_Ysmk_lower, ymax = CI_Ysmk_upper), 
                width = 1, 
                linetype = "dashed",
                linewidth = 0.25,
                color = "#D41159") +
  ylim(range(c(df$CI_Ysmk_lower, df$CI_Ysmk_upper))) +
  labs(x = "Ordered locations", y = "Response values", 
       title = "SMK Response Credible Intervals for new locations",
       subtitle = paste("Empirical coverage :", round(coverage_smk, 3))) +
  geom_point(aes(y = Ysmk_ord), pch = 15, size = 1.5, col = "#D41159") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

png("output/UC_10_1000.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(uc_Ybps, uc_Ysmk)
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS " = list("time" = elapsed_times,
                              "RMSPE" = rmspe_Y),
                "SMK" = list("time" = elapsed_times_smk,
                             "RMSPE" = rmspe_smk))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "output/simulation_univariate_10_1000.RData")

#####################################################################################################################################################