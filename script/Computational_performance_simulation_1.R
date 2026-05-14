#####################################################################################################################################################
## BPS vs MSMK COMPARISON - 5k - 500 #################################################################################################################################
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
p <- 2
q <- 2

# parameters
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- diag(q)
alfa <- 0.8
phi <- 4

set.seed(97)
# generate sintethic data
crd <- matrix(runif((n+u) * 2), ncol = 2)
X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
D <- spBPS:::arma_dist(crd)
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
hyperpar <- list(alpha = alfa_seq, phi = phi_seq)

# subsetting data
subset_size <- 500
K <- n/subset_size

# number of clusters for parallel implementation
n.core <- parallel::detectCores(logical=F)-1

# Multivariate BPS parallel fit -------------------------------------------------------

# fitting DBPS
out <- spBPS(data = list(Y = Y, X = X),
                priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                              V_r = diag(10, p),
                              Psi = diag(1, q),
                              nu = 3),
                coords = crd_s,
                hyperpar = hyperpar,
                subset_size = subset_size,
                cv_folds = 5,
                combine_method = "bps",
                draws = 200,
                n_cores = n.core,
                newdata = list(X = X_u, coords = crd_u),
                pred_batch_size = 50)

# BPS Results collection ----------------------------------------------------------

# statistics computations W
pred_mat_W <- out$predictive$Wu
post_mean_W <- apply(pred_mat_W, c(1,2), mean)
post_var_W <- apply(pred_mat_W, c(1,2), sd)
post_qnt_W <- apply(pred_mat_W, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- c(mean(W_u[,1] >= post_qnt_W[1,,1] & W_u[,1] <= post_qnt_W[2,,1]),
                mean(W_u[,2] >= post_qnt_W[1,,2] & W_u[,2] <= post_qnt_W[2,,2]))
cat("Empirical average coverage for Spatial process:", round(mean(coverage_W), 3))

# statistics computations Y
pred_mat_Y <- out$predictive$Yu
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(Y_u[,1] >= post_qnt_Y[1,,1] & Y_u[,1] <= post_qnt_Y[2,,1]),
                mean(Y_u[,2] >= post_qnt_Y[1,,2] & Y_u[,2] <= post_qnt_Y[2,,2]))
cat("Empirical average BPS coverage for Responses:", round(mean(coverage_Y), 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( colMeans( (W_u[-(1:m),] - post_mean_W[-(1:m),])^2 ) )); mean(rmspe_W)
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) )); mean(rmspe_Y)


# BPS Posterior inference -----------------------------------------------------

beta_smp <- out$posterior$beta
(post_mean_beta <-  apply(beta_smp, c(1,2), mean))
post_var_beta <- apply(beta_smp, c(1,2), sd)
post_low_beta <- apply(beta_smp, c(1,2), quantile, c(0.025))
post_upp_beta <- apply(beta_smp, c(1,2), quantile, c(0.975))

sigma_smp <- out$posterior$sigma
(post_mean_sigma <- apply(sigma_smp, c(1,2), mean))
post_var_sigma <- apply(sigma_smp, c(1,2), sd)
post_low_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025))
post_upp_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.975))

Wbps <- matrix(out$weights_global)
W_list <- out$weights_local
(post_mean_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps)
post_var_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# collecting
R <- 200
posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                        t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3],
                        cbind(rep(post_mean_hyp[2], R), rep(post_mean_hyp[2], R)))
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Multivariate SMK parallel fit --------------------------------------------------------

# ## Libraries
# library(mcmc)
# library(MASS)
# library(KernSmooth)
# library(fields)
# library(pscl)
# library(spBayes)
# library(mvtnorm)
# library(MCMCpack)
# library(Mposterior)
# library(parallel)
# library(doParallel)
# library(foreach)

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
per.core <- floor(n/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part)
## Number of MCMC iterations
mcmc.sample <- 500
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
y <- dat.nomiss[,3:4]
x <- dat.nomiss[,-(1:4)]
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
  Y.part[[i]] <- Y.train[index.part[[i]],]
  ## Predictor in i th subset
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset
  coords.part[[i]] <- coords.train[index.part[[i]],]
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){

  ## library
  library(spBayes)

  ## Response in subset i
  ZZ <- Y.part[[i]]
  ## Predictor in subset i
  XX <- X.part[[i]]
  ## Coordinates in subset i
  CC <- coords.part[[i]]

  ## priors
  q <- ncol(ZZ)
  A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

  starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
  tuning <- list("phi"=rep(1,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.01,q))
  priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
                 "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))
  cov.model <- "exponential"

  ## Model fitting
  # tictoc::tic()
  m.1 <- spMvLM(list(ZZ[,1]~XX-1, ZZ[,2]~XX-1),
                coords=CC, starting=starting, tuning=tuning, priors=priors,
                n.samples=mcmc.sample, cov.model=cov.model, verbose = F)
  # tictoc::toc()


  ## Recover all MCMC samples in each subset
  # tictoc::tic()
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=F)
  # tictoc::toc()

  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning
  gc()

  ## prediction covariates
  XXp <- dat.miss[,-c(1:4)]
  pred.covars <- mkMvX(list(XXp, XXp))

  # ## Prediction
  # tictoc::tic()
  Y.pred <- t(spPredict(m.1, coords.pred, pred.covars, start=0.5*mcmc.sample,verbose=T)$p.y.predictive.samples)
  # tictoc::toc()

  # Y.pred <- NULL
  # Y.pred[, 1] <- Y.predQ[,seq(1,ncol(Y.predQ),q)]
  # Y.pred[, 2] <- Y.predQ[,seq(2,ncol(Y.predQ),q)]

  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()

  # tic()

  ## Prediction
  # Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
  #                        pred.coords=coords.pred[-(1:m),],
  #                        start=0.5*mcmc.sample)$p.y.predictive.samples)
  # Y.pred <- NULL
  # for(l in 1:(length(p.sub)-1)){
  #   m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
  #                         pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
  #                         start=0.5*mcmc.sample)$p.y.predictive.samples
  #   Y.pred <- cbind(Y.pred,t(m.1.pred))
  # }
  #
  # Y.pred <- cbind(Y.pred, Yu.pred)

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

# MSMK Results collection ------------------------------------------------------------

# collecting predictions
Y.medQ <- Y.med
Y.med <- matrix(0, nrow(dat.miss), q)
Y.med[, 1] <- Y.medQ[seq(1,length(Y.medQ),q)]
Y.med[, 2] <- Y.medQ[seq(2,length(Y.medQ),q)]

Y.lower_qntQ <- Y.lower_qnt
Y.lower_qnt <- matrix(0, nrow(dat.miss), q)
Y.lower_qnt[, 1] <- Y.lower_qntQ[seq(1,length(Y.medQ),q)]
Y.lower_qnt[, 2] <- Y.lower_qntQ[seq(2,length(Y.medQ),q)]

Y.upper_qntQ <- Y.upper_qnt
Y.upper_qnt <- matrix(0, nrow(dat.miss), q)
Y.upper_qnt[, 1] <- Y.upper_qntQ[seq(1,length(Y.medQ),q)]
Y.upper_qnt[, 2] <- Y.upper_qntQ[seq(2,length(Y.medQ),q)]

# average coverage
coverage_smk <- c(mean(Y_u[,1] >= Y.lower_qnt[,1] & Y_u[,1] <= Y.upper_qnt[,1]),
                  mean(Y_u[,2] >= Y.lower_qnt[,2] & Y_u[,2] <= Y.upper_qnt[,2]))
cat("Empirical average MSMK coverage for Responses:", round(mean(coverage_smk), 3))
(CI_avlen_smk <- mean(Y.upper_qnt[,]-Y.lower_qnt[,]))

# Root Mean Square Prediction Error
(rmspe_smk <- sqrt( colMeans( (Y_u[-(1:m),] - Y.med[-(1:m),])^2 ) )); mean(rmspe_smk)


# MSMK Posterior inference -----------------------------------------------------

# collecting
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:n.core, function(a)obj[[a]]$atoms[,j]) %*% wts
}
posterior_smk <- post_smp_smk[,c(1:7,10:11)]
colnames(posterior_smk) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = out$timings[1],
                   "Combination" = out$timings[2],
                   "Prediction" = out$timings[3],
                   "Total time" = sum(out$timings))

elapsed_smk <- c("Fitting" = final.time,
                 "Combination" = recomb.time,
                 "Prediction" = quantile.calculation.time,
                 "Total time" = Time)

cat("minutes elapsed for BPS fully model-based uncertainty quantification : \n"); round(elapsed_times, 2)
cat("minutes elapsed for MSMK fully model-based uncertainty quantification : \n"); round(elapsed_smk, 2)


# Plotting data --------------------------------------------------

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(bayesplot)

mba_to_df <- function(coords, values, no.X = 500, no.Y = 500, h = 12) {
  surf <- MBA::mba.surf(cbind(coords, values), no.X = no.X, no.Y = no.Y,
                        exten = TRUE, sp = FALSE, h = h)$xyz.est
  data.frame(
    x = rep(surf$x, times = length(surf$y)),
    y = rep(surf$y, each  = length(surf$x)),
    z = as.vector(surf$z)
  )
}

df_Y1_true <- mba_to_df(crd_s, Y[, 1])
df_Y2_true <- mba_to_df(crd_s, Y[, 2])
df_Y1_bps  <- mba_to_df(crd_u, post_mean_Y[, 1])
df_Y2_bps  <- mba_to_df(crd_u, post_mean_Y[, 2])
df_Y1_smk  <- mba_to_df(crd_u, Y.med[, 1])
df_Y2_smk  <- mba_to_df(crd_u, Y.med[, 2])

lims_Y1 <- range(c(df_Y1_true$z, df_Y1_bps$z, df_Y1_smk$z), na.rm = TRUE)
lims_Y2 <- range(c(df_Y2_true$z, df_Y2_bps$z, df_Y2_smk$z), na.rm = TRUE)

rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu")[11:1])(100))

surf_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text         = element_text(size = 16, face = "bold"),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5, vjust = -0.35),
    plot.subtitle     = element_text(size = 14, hjust = 0.5, color = "grey30"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 14),
    legend.key.height = unit(1.2, "cm")
  )

plot_surf <- function(df, title, subtitle = NULL, limits) {
  ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colors = rdbu_pal, limits = limits,
                         name = "Value", na.value = NA) +
    coord_equal() +
    labs(x = "Easting", y = "Northing", title = title, subtitle = subtitle) +
    surf_theme
}

p_y1_true <- plot_surf(df_Y1_true, title = expression(Y[1] ~ "– True"),                                          limits = lims_Y1)
p_y1_bps  <- plot_surf(df_Y1_bps,  title = bquote(Y[1] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[1],  3)) ~ ")"),   limits = lims_Y1)
p_y1_smk  <- plot_surf(df_Y1_smk,  title = bquote(Y[1] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[1], 3)) ~ ")"),  limits = lims_Y1)
p_y2_true <- plot_surf(df_Y2_true, title = expression(Y[2] ~ "– True"),                                          limits = lims_Y2)
p_y2_bps  <- plot_surf(df_Y2_bps,  title = bquote(Y[2] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[2],  3)) ~ ")"),   limits = lims_Y2)
p_y2_smk  <- plot_surf(df_Y2_smk,  title = bquote(Y[2] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[2], 3)) ~ ")"),  limits = lims_Y2)

row1 <- (p_y1_true | p_y1_bps | p_y1_smk) + plot_layout(guides = "collect")
row2 <- (p_y2_true | p_y2_bps | p_y2_smk) + plot_layout(guides = "collect")
final_surf <- row1 / row2

png("surface_M_5_500.png", width = 16, height = 10, units = "in", res = 300, family = "sans")
print(final_surf)
dev.off()

uc_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text     = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 18, face = "bold", color = "grey30")
  )

make_uc <- function(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord, ylab, subtitle) {
  df <- data.frame(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord)
  ggplot(df, aes(x = x_ax, y = Yu_ord)) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
              width = 0, linetype = "dashed", linewidth = 0.3, color = "#D41159") +
    geom_point(pch = 18, size = 3.5, color = "#1A85FF") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, color = "#D41159") +
    ylim(range(c(CI_lower, CI_upper))) +
    labs(x = "Ordered locations", y = ylab, subtitle = subtitle) +
    uc_theme
}

ord1 <- order(Y_u[-(1:m), 1])
ord2 <- order(Y_u[-(1:m), 2])

uc_Y1bps <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  post_qnt_Y[1,-(1:m),1][ord1], post_qnt_Y[2,-(1:m),1][ord1],
  post_mean_Y[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[1], 3),")")
)
uc_Y2bps <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  post_qnt_Y[1,-(1:m),2][ord2], post_qnt_Y[2,-(1:m),2][ord2],
  post_mean_Y[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[2], 3),")")
)
uc_Y1smk <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  Y.lower_qnt[-(1:m),1][ord1], Y.upper_qnt[-(1:m),1][ord1],
  Y.med[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[1], 3),")")
)
uc_Y2smk <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  Y.lower_qnt[-(1:m),2][ord2], Y.upper_qnt[-(1:m),2][ord2],
  Y.med[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[2], 3),")")
)

png("UC_M_5_500.png", width = 16, height = 10,
    units = "in", res = 300, family = "sans")
print((uc_Y1bps | uc_Y2bps) / (uc_Y1smk | uc_Y2smk))
dev.off()


bayesplot_theme_set(theme_classic(base_size = 16))
color_scheme_set("viridis")

ci_labels <- c(
  expression(beta[list(1,1)]), expression(beta[list(2,1)]),
  expression(beta[list(1,2)]), expression(beta[list(2,2)]),
  expression(phi[1]),          expression(phi[2]),
  expression(Sigma[list(1,1)]),expression(Sigma[list(2,1)]),expression(Sigma[list(2,2)])
)

ci_theme_extra <- theme(
  axis.text.x     = element_text(size = 18, face = "bold", angle = 45, hjust = 1),
  axis.text.y     = element_text(size = 16, face = "bold"),
  axis.title      = element_text(size = 16, face = "bold"),
  legend.text     = element_text(size = 14, face = "bold"),
  legend.title    = element_blank(),
  legend.position = "right",
  plot.title      = element_blank(),
  plot.subtitle   = element_text(size = 18, face = "bold", color = "grey30")
)

post_int_bps <- mcmc_recover_intervals(
  posterior_bps, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "DBPS (posterior means and 95% CI)") +
  ci_theme_extra

post_int_smk <- mcmc_recover_intervals(
  posterior_smk, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "MSMK (posterior means and 95% CI)") +
  ci_theme_extra

png("CIpost_M_5_500.png", width = 14, height = 5, units = "in", res = 300, family = "sans")
print((post_int_bps | post_int_smk) + plot_layout(guides = "collect"))
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS" = list("time"    = elapsed_times,
                             "metrics" = c("RMSPE" = rmspe_Y, "CI_len" = CI_avlen_bps, "ECoverage" = coverage_Y)),
                "SMK" = list("time"    = elapsed_smk,
                             "metrics" = c("RMSPE" = rmspe_smk, "CI_len" = CI_avlen_smk, "ECoverage" = coverage_smk)))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "simulation_multivariate_5_500.RData")

#####################################################################################################################################################
## BPS vs MSMK COMPARISON - 5k - 1000 #################################################################################################################################
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

# Sys.setenv(PKG_CXXFLAGS = "-Ofast"); sourceCpp("new_code.cpp")
# ADD TO THE PACKAGE NEW FUNCTION IN new_code.cpp

# Data generation -------------------------------------------------------------

# dimensions
n <- 5000
m <- 500
u <- 250
p <- 2
q <- 2

# parameters
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- diag(q)
alfa <- 0.8
phi <- 4

set.seed(97)
# generate sintethic data
crd <- matrix(runif((n+u) * 2), ncol = 2)
X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
D <- spBPS:::arma_dist(crd)
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
hyperpar <- list(alpha = alfa_seq, phi = phi_seq)

# subsetting data
subset_size <- 1000
K <- n/subset_size

# number of clusters for parallel implementation
n.core <- parallel::detectCores(logical=F)-1


# Multivariate BPS parallel fit -------------------------------------------------------

# fitting DBPS
out <- spBPS(data = list(Y = Y, X = X),
             priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                           V_r = diag(10, p),
                           Psi = diag(1, q),
                           nu = 3),
             coords = crd_s,
             hyperpar = hyperpar,
             subset_size = subset_size,
             cv_folds = 5,
             combine_method = "bps",
             draws = 200,
             n_cores = n.core,
             newdata = list(X = X_u, coords = crd_u),
             pred_batch_size = 50)

# BPS Results collection ----------------------------------------------------------

# statistics computations W
pred_mat_W <- out$predictive$Wu
post_mean_W <- apply(pred_mat_W, c(1,2), mean)
post_var_W <- apply(pred_mat_W, c(1,2), sd)
post_qnt_W <- apply(pred_mat_W, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- c(mean(W_u[,1] >= post_qnt_W[1,,1] & W_u[,1] <= post_qnt_W[2,,1]),
                mean(W_u[,2] >= post_qnt_W[1,,2] & W_u[,2] <= post_qnt_W[2,,2]))
cat("Empirical average coverage for Spatial process:", round(mean(coverage_W), 3))

# statistics computations Y
pred_mat_Y <- out$predictive$Yu
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(Y_u[,1] >= post_qnt_Y[1,,1] & Y_u[,1] <= post_qnt_Y[2,,1]),
                mean(Y_u[,2] >= post_qnt_Y[1,,2] & Y_u[,2] <= post_qnt_Y[2,,2]))
cat("Empirical average BPS coverage for Responses:", round(mean(coverage_Y), 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( colMeans( (W_u[-(1:m),] - post_mean_W[-(1:m),])^2 ) )); mean(rmspe_W)
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) )); mean(rmspe_Y)


# BPS Posterior inference -----------------------------------------------------

beta_smp <- out$posterior$beta
(post_mean_beta <-  apply(beta_smp, c(1,2), mean))
post_var_beta <- apply(beta_smp, c(1,2), sd)
post_low_beta <- apply(beta_smp, c(1,2), quantile, c(0.025))
post_upp_beta <- apply(beta_smp, c(1,2), quantile, c(0.975))

sigma_smp <- out$posterior$sigma
(post_mean_sigma <- apply(sigma_smp, c(1,2), mean))
post_var_sigma <- apply(sigma_smp, c(1,2), sd)
post_low_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025))
post_upp_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.975))

Wbps <- matrix(out$weights_global)
W_list <- out$weights_local
(post_mean_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps)
post_var_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# collecting
R <- 200
posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                       t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3],
                       cbind(rep(post_mean_hyp[2], R), rep(post_mean_hyp[2], R)))
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Multivariate SMK parallel fit --------------------------------------------------------

# ## Libraries
# library(mcmc)
# library(MASS)
# library(KernSmooth)
# library(fields)
# library(pscl)
# library(spBayes)
# library(mvtnorm)
# library(MCMCpack)
# library(Mposterior)
# library(parallel)
# library(doParallel)
# library(foreach)

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
per.core <- floor(n/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part)
## Number of MCMC iterations
mcmc.sample <- 500
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
y <- dat.nomiss[,3:4]
x <- dat.nomiss[,-(1:4)]
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
  Y.part[[i]] <- Y.train[index.part[[i]],]
  ## Predictor in i th subset
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset
  coords.part[[i]] <- coords.train[index.part[[i]],]
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){
  
  ## library
  library(spBayes)
  
  ## Response in subset i
  ZZ <- Y.part[[i]]
  ## Predictor in subset i
  XX <- X.part[[i]]
  ## Coordinates in subset i
  CC <- coords.part[[i]]
  
  ## priors
  q <- ncol(ZZ)
  A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
  
  starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
  tuning <- list("phi"=rep(1,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.01,q))
  priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
                 "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))
  cov.model <- "exponential"
  
  ## Model fitting
  # tictoc::tic()
  m.1 <- spMvLM(list(ZZ[,1]~XX-1, ZZ[,2]~XX-1),
                coords=CC, starting=starting, tuning=tuning, priors=priors,
                n.samples=mcmc.sample, cov.model=cov.model, verbose = F)
  # tictoc::toc()
  
  
  ## Recover all MCMC samples in each subset
  # tictoc::tic()
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=F)
  # tictoc::toc()
  
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning
  gc()
  
  ## prediction covariates
  XXp <- dat.miss[,-c(1:4)]
  pred.covars <- mkMvX(list(XXp, XXp))
  
  # ## Prediction
  # tictoc::tic()
  Y.pred <- t(spPredict(m.1, coords.pred, pred.covars, start=0.5*mcmc.sample,verbose=T)$p.y.predictive.samples)
  # tictoc::toc()
  
  # Y.pred <- NULL
  # Y.pred[, 1] <- Y.predQ[,seq(1,ncol(Y.predQ),q)]
  # Y.pred[, 2] <- Y.predQ[,seq(2,ncol(Y.predQ),q)]
  
  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()
  
  # tic()
  
  ## Prediction
  # Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
  #                        pred.coords=coords.pred[-(1:m),],
  #                        start=0.5*mcmc.sample)$p.y.predictive.samples)
  # Y.pred <- NULL
  # for(l in 1:(length(p.sub)-1)){
  #   m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
  #                         pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
  #                         start=0.5*mcmc.sample)$p.y.predictive.samples
  #   Y.pred <- cbind(Y.pred,t(m.1.pred))
  # }
  #
  # Y.pred <- cbind(Y.pred, Yu.pred)
  
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

# MSMK Results collection ------------------------------------------------------------

# collecting predictions
Y.medQ <- Y.med
Y.med <- matrix(0, nrow(dat.miss), q)
Y.med[, 1] <- Y.medQ[seq(1,length(Y.medQ),q)]
Y.med[, 2] <- Y.medQ[seq(2,length(Y.medQ),q)]

Y.lower_qntQ <- Y.lower_qnt
Y.lower_qnt <- matrix(0, nrow(dat.miss), q)
Y.lower_qnt[, 1] <- Y.lower_qntQ[seq(1,length(Y.medQ),q)]
Y.lower_qnt[, 2] <- Y.lower_qntQ[seq(2,length(Y.medQ),q)]

Y.upper_qntQ <- Y.upper_qnt
Y.upper_qnt <- matrix(0, nrow(dat.miss), q)
Y.upper_qnt[, 1] <- Y.upper_qntQ[seq(1,length(Y.medQ),q)]
Y.upper_qnt[, 2] <- Y.upper_qntQ[seq(2,length(Y.medQ),q)]

# average coverage
coverage_smk <- c(mean(Y_u[,1] >= Y.lower_qnt[,1] & Y_u[,1] <= Y.upper_qnt[,1]),
                  mean(Y_u[,2] >= Y.lower_qnt[,2] & Y_u[,2] <= Y.upper_qnt[,2]))
cat("Empirical average MSMK coverage for Responses:", round(mean(coverage_smk), 3))
(CI_avlen_smk <- mean(Y.upper_qnt[,]-Y.lower_qnt[,]))

# Root Mean Square Prediction Error
(rmspe_smk <- sqrt( colMeans( (Y_u[-(1:m),] - Y.med[-(1:m),])^2 ) )); mean(rmspe_smk)


# MSMK Posterior inference -----------------------------------------------------

# collecting
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:n.core, function(a)obj[[a]]$atoms[,j]) %*% wts
}
posterior_smk <- post_smp_smk[,c(1:7,10:11)]
colnames(posterior_smk) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = out$timings[1],
                   "Combination" = out$timings[2],
                   "Prediction" = out$timings[3],
                   "Total time" = sum(out$timings))

elapsed_smk <- c("Fitting" = final.time,
                 "Combination" = recomb.time,
                 "Prediction" = quantile.calculation.time,
                 "Total time" = Time)

cat("minutes elapsed for BPS fully model-based uncertainty quantification : \n"); round(elapsed_times, 2)
cat("minutes elapsed for MSMK fully model-based uncertainty quantification : \n"); round(elapsed_smk, 2)


# Plotting data --------------------------------------------------

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(bayesplot)

mba_to_df <- function(coords, values, no.X = 500, no.Y = 500, h = 12) {
  surf <- MBA::mba.surf(cbind(coords, values), no.X = no.X, no.Y = no.Y,
                        exten = TRUE, sp = FALSE, h = h)$xyz.est
  data.frame(
    x = rep(surf$x, times = length(surf$y)),
    y = rep(surf$y, each  = length(surf$x)),
    z = as.vector(surf$z)
  )
}

df_Y1_true <- mba_to_df(crd_s, Y[, 1])
df_Y2_true <- mba_to_df(crd_s, Y[, 2])
df_Y1_bps  <- mba_to_df(crd_u, post_mean_Y[, 1])
df_Y2_bps  <- mba_to_df(crd_u, post_mean_Y[, 2])
df_Y1_smk  <- mba_to_df(crd_u, Y.med[, 1])
df_Y2_smk  <- mba_to_df(crd_u, Y.med[, 2])

lims_Y1 <- range(c(df_Y1_true$z, df_Y1_bps$z, df_Y1_smk$z), na.rm = TRUE)
lims_Y2 <- range(c(df_Y2_true$z, df_Y2_bps$z, df_Y2_smk$z), na.rm = TRUE)

rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu")[11:1])(100))

surf_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text         = element_text(size = 16, face = "bold"),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5, vjust = -0.35),
    plot.subtitle     = element_text(size = 14, hjust = 0.5, color = "grey30"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 14),
    legend.key.height = unit(1.2, "cm")
  )

plot_surf <- function(df, title, subtitle = NULL, limits) {
  ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colors = rdbu_pal, limits = limits,
                         name = "Value", na.value = NA) +
    coord_equal() +
    labs(x = "Easting", y = "Northing", title = title, subtitle = subtitle) +
    surf_theme
}

p_y1_true <- plot_surf(df_Y1_true, title = expression(Y[1] ~ "– True"),                                          limits = lims_Y1)
p_y1_bps  <- plot_surf(df_Y1_bps,  title = bquote(Y[1] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[1],  3)) ~ ")"),   limits = lims_Y1)
p_y1_smk  <- plot_surf(df_Y1_smk,  title = bquote(Y[1] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[1], 3)) ~ ")"),  limits = lims_Y1)
p_y2_true <- plot_surf(df_Y2_true, title = expression(Y[2] ~ "– True"),                                          limits = lims_Y2)
p_y2_bps  <- plot_surf(df_Y2_bps,  title = bquote(Y[2] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[2],  3)) ~ ")"),   limits = lims_Y2)
p_y2_smk  <- plot_surf(df_Y2_smk,  title = bquote(Y[2] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[2], 3)) ~ ")"),  limits = lims_Y2)

row1 <- (p_y1_true | p_y1_bps | p_y1_smk) + plot_layout(guides = "collect")
row2 <- (p_y2_true | p_y2_bps | p_y2_smk) + plot_layout(guides = "collect")
final_surf <- row1 / row2

png("surface_M_5_1000.png", width = 16, height = 10, units = "in", res = 300, family = "sans")
print(final_surf)
dev.off()

uc_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text     = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 18, face = "bold", color = "grey30")
  )

make_uc <- function(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord, ylab, subtitle) {
  df <- data.frame(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord)
  ggplot(df, aes(x = x_ax, y = Yu_ord)) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  width = 0, linetype = "dashed", linewidth = 0.3, color = "#D41159") +
    geom_point(pch = 18, size = 3.5, color = "#1A85FF") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, color = "#D41159") +
    ylim(range(c(CI_lower, CI_upper))) +
    labs(x = "Ordered locations", y = ylab, subtitle = subtitle) +
    uc_theme
}

ord1 <- order(Y_u[-(1:m), 1])
ord2 <- order(Y_u[-(1:m), 2])

uc_Y1bps <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  post_qnt_Y[1,-(1:m),1][ord1], post_qnt_Y[2,-(1:m),1][ord1],
  post_mean_Y[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[1], 3),")")
)
uc_Y2bps <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  post_qnt_Y[1,-(1:m),2][ord2], post_qnt_Y[2,-(1:m),2][ord2],
  post_mean_Y[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[2], 3),")")
)
uc_Y1smk <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  Y.lower_qnt[-(1:m),1][ord1], Y.upper_qnt[-(1:m),1][ord1],
  Y.med[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[1], 3),")")
)
uc_Y2smk <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  Y.lower_qnt[-(1:m),2][ord2], Y.upper_qnt[-(1:m),2][ord2],
  Y.med[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[2], 3),")")
)

png("UC_M_5_1000.png", width = 16, height = 10,
    units = "in", res = 300, family = "sans")
print((uc_Y1bps | uc_Y2bps) / (uc_Y1smk | uc_Y2smk))
dev.off()

bayesplot_theme_set(theme_classic(base_size = 16))
color_scheme_set("viridis")

ci_labels <- c(
  expression(beta[list(1,1)]), expression(beta[list(2,1)]),
  expression(beta[list(1,2)]), expression(beta[list(2,2)]),
  expression(phi[1]),          expression(phi[2]),
  expression(Sigma[list(1,1)]),expression(Sigma[list(2,1)]),expression(Sigma[list(2,2)])
)

ci_theme_extra <- theme(
  axis.text.x     = element_text(size = 18, face = "bold", angle = 45, hjust = 1),
  axis.text.y     = element_text(size = 16, face = "bold"),
  axis.title      = element_text(size = 16, face = "bold"),
  legend.text     = element_text(size = 14, face = "bold"),
  legend.title    = element_blank(),
  legend.position = "right",
  plot.title      = element_blank(),
  plot.subtitle   = element_text(size = 18, face = "bold", color = "grey30")
)

post_int_bps <- mcmc_recover_intervals(
  posterior_bps, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "DBPS (posterior means and 95% CI)") +
  ci_theme_extra

post_int_smk <- mcmc_recover_intervals(
  posterior_smk, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "MSMK (posterior means and 95% CI)") +
  ci_theme_extra

png("CIpost_M_5_1000.png", width = 14, height = 5, units = "in", res = 300, family = "sans")
print((post_int_bps | post_int_smk) + plot_layout(guides = "collect"))
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS" = list("time"    = elapsed_times,
                             "metrics" = c("RMSPE" = rmspe_Y, "CI_len" = CI_avlen_bps, "ECoverage" = coverage_Y)),
                "SMK" = list("time"    = elapsed_smk,
                             "metrics" = c("RMSPE" = rmspe_smk, "CI_len" = CI_avlen_smk, "ECoverage" = coverage_smk)))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "simulation_multivariate_5_1000.RData")

#####################################################################################################################################################
## BPS vs MSMK COMPARISON - 10k - 500 #################################################################################################################################
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

# Sys.setenv(PKG_CXXFLAGS = "-Ofast"); sourceCpp("new_code.cpp")
# ADD TO THE PACKAGE NEW FUNCTION IN new_code.cpp

# Data generation -------------------------------------------------------------

# dimensions
n <- 10000
m <- 500
u <- 250
p <- 2
q <- 2

# parameters
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- diag(q)
alfa <- 0.8
phi <- 4

set.seed(97)
# generate sintethic data
crd <- matrix(runif((n+u) * 2), ncol = 2)
X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
D <- spBPS:::arma_dist(crd)
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
hyperpar <- list(alpha = alfa_seq, phi = phi_seq)

# subsetting data
subset_size <- 500
K <- n/subset_size

# number of clusters for parallel implementation
n.core <- parallel::detectCores(logical=F)-1

# Multivariate BPS parallel fit -------------------------------------------------------

# fitting DBPS
out <- spBPS(data = list(Y = Y, X = X),
             priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                           V_r = diag(10, p),
                           Psi = diag(1, q),
                           nu = 3),
             coords = crd_s,
             hyperpar = hyperpar,
             subset_size = subset_size,
             cv_folds = 5,
             combine_method = "bps",
             draws = 200,
             n_cores = n.core,
             newdata = list(X = X_u, coords = crd_u),
             pred_batch_size = 50)

# BPS Results collection ----------------------------------------------------------

# statistics computations W
pred_mat_W <- out$predictive$Wu
post_mean_W <- apply(pred_mat_W, c(1,2), mean)
post_var_W <- apply(pred_mat_W, c(1,2), sd)
post_qnt_W <- apply(pred_mat_W, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- c(mean(W_u[,1] >= post_qnt_W[1,,1] & W_u[,1] <= post_qnt_W[2,,1]),
                mean(W_u[,2] >= post_qnt_W[1,,2] & W_u[,2] <= post_qnt_W[2,,2]))
cat("Empirical average coverage for Spatial process:", round(mean(coverage_W), 3))

# statistics computations Y
pred_mat_Y <- out$predictive$Yu
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(Y_u[,1] >= post_qnt_Y[1,,1] & Y_u[,1] <= post_qnt_Y[2,,1]),
                mean(Y_u[,2] >= post_qnt_Y[1,,2] & Y_u[,2] <= post_qnt_Y[2,,2]))
cat("Empirical average BPS coverage for Responses:", round(mean(coverage_Y), 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( colMeans( (W_u[-(1:m),] - post_mean_W[-(1:m),])^2 ) )); mean(rmspe_W)
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) )); mean(rmspe_Y)


# BPS Posterior inference -----------------------------------------------------

beta_smp <- out$posterior$beta
(post_mean_beta <-  apply(beta_smp, c(1,2), mean))
post_var_beta <- apply(beta_smp, c(1,2), sd)
post_low_beta <- apply(beta_smp, c(1,2), quantile, c(0.025))
post_upp_beta <- apply(beta_smp, c(1,2), quantile, c(0.975))

sigma_smp <- out$posterior$sigma
(post_mean_sigma <- apply(sigma_smp, c(1,2), mean))
post_var_sigma <- apply(sigma_smp, c(1,2), sd)
post_low_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025))
post_upp_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.975))

Wbps <- matrix(out$weights_global)
W_list <- out$weights_local
(post_mean_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps)
post_var_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# collecting
R <- 200
posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                       t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3],
                       cbind(rep(post_mean_hyp[2], R), rep(post_mean_hyp[2], R)))
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Multivariate SMK parallel fit --------------------------------------------------------

# ## Libraries
# library(mcmc)
# library(MASS)
# library(KernSmooth)
# library(fields)
# library(pscl)
# library(spBayes)
# library(mvtnorm)
# library(MCMCpack)
# library(Mposterior)
# library(parallel)
# library(doParallel)
# library(foreach)

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
per.core <- floor(n/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part)
## Number of MCMC iterations
mcmc.sample <- 500
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
y <- dat.nomiss[,3:4]
x <- dat.nomiss[,-(1:4)]
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
  Y.part[[i]] <- Y.train[index.part[[i]],]
  ## Predictor in i th subset
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset
  coords.part[[i]] <- coords.train[index.part[[i]],]
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){
  
  ## library
  library(spBayes)
  
  ## Response in subset i
  ZZ <- Y.part[[i]]
  ## Predictor in subset i
  XX <- X.part[[i]]
  ## Coordinates in subset i
  CC <- coords.part[[i]]
  
  ## priors
  q <- ncol(ZZ)
  A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
  
  starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
  tuning <- list("phi"=rep(1,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.01,q))
  priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
                 "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))
  cov.model <- "exponential"
  
  ## Model fitting
  # tictoc::tic()
  m.1 <- spMvLM(list(ZZ[,1]~XX-1, ZZ[,2]~XX-1),
                coords=CC, starting=starting, tuning=tuning, priors=priors,
                n.samples=mcmc.sample, cov.model=cov.model, verbose = F)
  # tictoc::toc()
  
  
  ## Recover all MCMC samples in each subset
  # tictoc::tic()
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=F)
  # tictoc::toc()
  
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning
  gc()
  
  ## prediction covariates
  XXp <- dat.miss[,-c(1:4)]
  pred.covars <- mkMvX(list(XXp, XXp))
  
  # ## Prediction
  # tictoc::tic()
  Y.pred <- t(spPredict(m.1, coords.pred, pred.covars, start=0.5*mcmc.sample,verbose=T)$p.y.predictive.samples)
  # tictoc::toc()
  
  # Y.pred <- NULL
  # Y.pred[, 1] <- Y.predQ[,seq(1,ncol(Y.predQ),q)]
  # Y.pred[, 2] <- Y.predQ[,seq(2,ncol(Y.predQ),q)]
  
  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()
  
  # tic()
  
  ## Prediction
  # Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
  #                        pred.coords=coords.pred[-(1:m),],
  #                        start=0.5*mcmc.sample)$p.y.predictive.samples)
  # Y.pred <- NULL
  # for(l in 1:(length(p.sub)-1)){
  #   m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
  #                         pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
  #                         start=0.5*mcmc.sample)$p.y.predictive.samples
  #   Y.pred <- cbind(Y.pred,t(m.1.pred))
  # }
  #
  # Y.pred <- cbind(Y.pred, Yu.pred)
  
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

# MSMK Results collection ------------------------------------------------------------

# collecting predictions
Y.medQ <- Y.med
Y.med <- matrix(0, nrow(dat.miss), q)
Y.med[, 1] <- Y.medQ[seq(1,length(Y.medQ),q)]
Y.med[, 2] <- Y.medQ[seq(2,length(Y.medQ),q)]

Y.lower_qntQ <- Y.lower_qnt
Y.lower_qnt <- matrix(0, nrow(dat.miss), q)
Y.lower_qnt[, 1] <- Y.lower_qntQ[seq(1,length(Y.medQ),q)]
Y.lower_qnt[, 2] <- Y.lower_qntQ[seq(2,length(Y.medQ),q)]

Y.upper_qntQ <- Y.upper_qnt
Y.upper_qnt <- matrix(0, nrow(dat.miss), q)
Y.upper_qnt[, 1] <- Y.upper_qntQ[seq(1,length(Y.medQ),q)]
Y.upper_qnt[, 2] <- Y.upper_qntQ[seq(2,length(Y.medQ),q)]

# average coverage
coverage_smk <- c(mean(Y_u[,1] >= Y.lower_qnt[,1] & Y_u[,1] <= Y.upper_qnt[,1]),
                  mean(Y_u[,2] >= Y.lower_qnt[,2] & Y_u[,2] <= Y.upper_qnt[,2]))
cat("Empirical average MSMK coverage for Responses:", round(mean(coverage_smk), 3))
(CI_avlen_smk <- mean(Y.upper_qnt[,]-Y.lower_qnt[,]))

# Root Mean Square Prediction Error
(rmspe_smk <- sqrt( colMeans( (Y_u[-(1:m),] - Y.med[-(1:m),])^2 ) )); mean(rmspe_smk)


# MSMK Posterior inference -----------------------------------------------------

# collecting
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:n.core, function(a)obj[[a]]$atoms[,j]) %*% wts
}
posterior_smk <- post_smp_smk[,c(1:7,10:11)]
colnames(posterior_smk) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = out$timings[1],
                   "Combination" = out$timings[2],
                   "Prediction" = out$timings[3],
                   "Total time" = sum(out$timings))

elapsed_smk <- c("Fitting" = final.time,
                 "Combination" = recomb.time,
                 "Prediction" = quantile.calculation.time,
                 "Total time" = Time)

cat("minutes elapsed for BPS fully model-based uncertainty quantification : \n"); round(elapsed_times, 2)
cat("minutes elapsed for MSMK fully model-based uncertainty quantification : \n"); round(elapsed_smk, 2)


# Plotting data --------------------------------------------------

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(bayesplot)

mba_to_df <- function(coords, values, no.X = 500, no.Y = 500, h = 12) {
  surf <- MBA::mba.surf(cbind(coords, values), no.X = no.X, no.Y = no.Y,
                        exten = TRUE, sp = FALSE, h = h)$xyz.est
  data.frame(
    x = rep(surf$x, times = length(surf$y)),
    y = rep(surf$y, each  = length(surf$x)),
    z = as.vector(surf$z)
  )
}

df_Y1_true <- mba_to_df(crd_s, Y[, 1])
df_Y2_true <- mba_to_df(crd_s, Y[, 2])
df_Y1_bps  <- mba_to_df(crd_u, post_mean_Y[, 1])
df_Y2_bps  <- mba_to_df(crd_u, post_mean_Y[, 2])
df_Y1_smk  <- mba_to_df(crd_u, Y.med[, 1])
df_Y2_smk  <- mba_to_df(crd_u, Y.med[, 2])

lims_Y1 <- range(c(df_Y1_true$z, df_Y1_bps$z, df_Y1_smk$z), na.rm = TRUE)
lims_Y2 <- range(c(df_Y2_true$z, df_Y2_bps$z, df_Y2_smk$z), na.rm = TRUE)

rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu")[11:1])(100))

surf_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text         = element_text(size = 16, face = "bold"),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5, vjust = -0.35),
    plot.subtitle     = element_text(size = 14, hjust = 0.5, color = "grey30"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 14),
    legend.key.height = unit(1.2, "cm")
  )

plot_surf <- function(df, title, subtitle = NULL, limits) {
  ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colors = rdbu_pal, limits = limits,
                         name = "Value", na.value = NA) +
    coord_equal() +
    labs(x = "Easting", y = "Northing", title = title, subtitle = subtitle) +
    surf_theme
}

p_y1_true <- plot_surf(df_Y1_true, title = expression(Y[1] ~ "– True"),                                          limits = lims_Y1)
p_y1_bps  <- plot_surf(df_Y1_bps,  title = bquote(Y[1] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[1],  3)) ~ ")"),   limits = lims_Y1)
p_y1_smk  <- plot_surf(df_Y1_smk,  title = bquote(Y[1] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[1], 3)) ~ ")"),  limits = lims_Y1)
p_y2_true <- plot_surf(df_Y2_true, title = expression(Y[2] ~ "– True"),                                          limits = lims_Y2)
p_y2_bps  <- plot_surf(df_Y2_bps,  title = bquote(Y[2] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[2],  3)) ~ ")"),   limits = lims_Y2)
p_y2_smk  <- plot_surf(df_Y2_smk,  title = bquote(Y[2] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[2], 3)) ~ ")"),  limits = lims_Y2)

row1 <- (p_y1_true | p_y1_bps | p_y1_smk) + plot_layout(guides = "collect")
row2 <- (p_y2_true | p_y2_bps | p_y2_smk) + plot_layout(guides = "collect")
final_surf <- row1 / row2

png("surface_M_10_500.png", width = 16, height = 10, units = "in", res = 300, family = "sans")
print(final_surf)
dev.off()

uc_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text     = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 18, face = "bold", color = "grey30")
  )

make_uc <- function(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord, ylab, subtitle) {
  df <- data.frame(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord)
  ggplot(df, aes(x = x_ax, y = Yu_ord)) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  width = 0, linetype = "dashed", linewidth = 0.3, color = "#D41159") +
    geom_point(pch = 18, size = 3.5, color = "#1A85FF") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, color = "#D41159") +
    ylim(range(c(CI_lower, CI_upper))) +
    labs(x = "Ordered locations", y = ylab, subtitle = subtitle) +
    uc_theme
}

ord1 <- order(Y_u[-(1:m), 1])
ord2 <- order(Y_u[-(1:m), 2])

uc_Y1bps <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  post_qnt_Y[1,-(1:m),1][ord1], post_qnt_Y[2,-(1:m),1][ord1],
  post_mean_Y[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[1], 3),")")
)
uc_Y2bps <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  post_qnt_Y[1,-(1:m),2][ord2], post_qnt_Y[2,-(1:m),2][ord2],
  post_mean_Y[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[2], 3),")")
)
uc_Y1smk <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  Y.lower_qnt[-(1:m),1][ord1], Y.upper_qnt[-(1:m),1][ord1],
  Y.med[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[1], 3),")")
)
uc_Y2smk <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  Y.lower_qnt[-(1:m),2][ord2], Y.upper_qnt[-(1:m),2][ord2],
  Y.med[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[2], 3),")")
)

png("UC_M_10_500.png", width = 16, height = 10,
    units = "in", res = 300, family = "sans")
print((uc_Y1bps | uc_Y2bps) / (uc_Y1smk | uc_Y2smk))
dev.off()

bayesplot_theme_set(theme_classic(base_size = 16))
color_scheme_set("viridis")

ci_labels <- c(
  expression(beta[list(1,1)]), expression(beta[list(2,1)]),
  expression(beta[list(1,2)]), expression(beta[list(2,2)]),
  expression(phi[1]),          expression(phi[2]),
  expression(Sigma[list(1,1)]),expression(Sigma[list(2,1)]),expression(Sigma[list(2,2)])
)

ci_theme_extra <- theme(
  axis.text.x     = element_text(size = 18, face = "bold", angle = 45, hjust = 1),
  axis.text.y     = element_text(size = 16, face = "bold"),
  axis.title      = element_text(size = 16, face = "bold"),
  legend.text     = element_text(size = 14, face = "bold"),
  legend.title    = element_blank(),
  legend.position = "right",
  plot.title      = element_blank(),
  plot.subtitle   = element_text(size = 18, face = "bold", color = "grey30")
)

post_int_bps <- mcmc_recover_intervals(
  posterior_bps, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "DBPS (posterior means and 95% CI)") +
  ci_theme_extra

post_int_smk <- mcmc_recover_intervals(
  posterior_smk, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "MSMK (posterior means and 95% CI)") +
  ci_theme_extra

png("CIpost_M_10_500.png", width = 14, height = 5, units = "in", res = 300, family = "sans")
print((post_int_bps | post_int_smk) + plot_layout(guides = "collect"))
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS" = list("time"    = elapsed_times,
                             "metrics" = c("RMSPE" = rmspe_Y, "CI_len" = CI_avlen_bps, "ECoverage" = coverage_Y)),
                "SMK" = list("time"    = elapsed_smk,
                             "metrics" = c("RMSPE" = rmspe_smk, "CI_len" = CI_avlen_smk, "ECoverage" = coverage_smk)))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "simulation_multivariate_10_500.RData")

#####################################################################################################################################################
## BPS vs MSMK COMPARISON - 10k - 1000 #################################################################################################################################
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

# Sys.setenv(PKG_CXXFLAGS = "-Ofast"); sourceCpp("new_code.cpp")
# ADD TO THE PACKAGE NEW FUNCTION IN new_code.cpp

# Data generation -------------------------------------------------------------

# dimensions
n <- 10000
m <- 500
u <- 250
p <- 2
q <- 2

# parameters
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- diag(q)
alfa <- 0.8
phi <- 4

set.seed(97)
# generate sintethic data
crd <- matrix(runif((n+u) * 2), ncol = 2)
X_or <- cbind(rep(1, n+u), matrix(runif((p-1)*(n+u)), ncol = (p-1)))
D <- spBPS:::arma_dist(crd)
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
hyperpar <- list(alpha = alfa_seq, phi = phi_seq)

# subsetting data
subset_size <- 1000
K <- n/subset_size

# number of clusters for parallel implementation
n.core <- parallel::detectCores(logical=F)-1


# Multivariate BPS parallel fit -------------------------------------------------------

# fitting DBPS
out <- spBPS(data = list(Y = Y, X = X),
             priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                           V_r = diag(10, p),
                           Psi = diag(1, q),
                           nu = 3),
             coords = crd_s,
             hyperpar = hyperpar,
             subset_size = subset_size,
             cv_folds = 5,
             combine_method = "bps",
             draws = 200,
             n_cores = n.core,
             newdata = list(X = X_u, coords = crd_u),
             pred_batch_size = 50)

# BPS Results collection ----------------------------------------------------------

# statistics computations W
pred_mat_W <- out$predictive$Wu
post_mean_W <- apply(pred_mat_W, c(1,2), mean)
post_var_W <- apply(pred_mat_W, c(1,2), sd)
post_qnt_W <- apply(pred_mat_W, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for W
coverage_W <- c(mean(W_u[,1] >= post_qnt_W[1,,1] & W_u[,1] <= post_qnt_W[2,,1]),
                mean(W_u[,2] >= post_qnt_W[1,,2] & W_u[,2] <= post_qnt_W[2,,2]))
cat("Empirical average coverage for Spatial process:", round(mean(coverage_W), 3))

# statistics computations Y
pred_mat_Y <- out$predictive$Yu
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(Y_u[,1] >= post_qnt_Y[1,,1] & Y_u[,1] <= post_qnt_Y[2,,1]),
                mean(Y_u[,2] >= post_qnt_Y[1,,2] & Y_u[,2] <= post_qnt_Y[2,,2]))
cat("Empirical average BPS coverage for Responses:", round(mean(coverage_Y), 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
(rmspe_W <- sqrt( colMeans( (W_u[-(1:m),] - post_mean_W[-(1:m),])^2 ) )); mean(rmspe_W)
(rmspe_Y <- sqrt( colMeans( (Y_u[-(1:m),] - post_mean_Y[-(1:m),])^2 ) )); mean(rmspe_Y)


# BPS Posterior inference -----------------------------------------------------

beta_smp <- out$posterior$beta
(post_mean_beta <-  apply(beta_smp, c(1,2), mean))
post_var_beta <- apply(beta_smp, c(1,2), sd)
post_low_beta <- apply(beta_smp, c(1,2), quantile, c(0.025))
post_upp_beta <- apply(beta_smp, c(1,2), quantile, c(0.975))

sigma_smp <- out$posterior$sigma
(post_mean_sigma <- apply(sigma_smp, c(1,2), mean))
post_var_sigma <- apply(sigma_smp, c(1,2), sd)
post_low_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025))
post_upp_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.975))

Wbps <- matrix(out$weights_global)
W_list <- out$weights_local
(post_mean_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps)
post_var_hyp <- sapply(1:K, function(k)t(spBPS:::expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# collecting
R <- 200
posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                       t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3],
                       cbind(rep(post_mean_hyp[2], R), rep(post_mean_hyp[2], R)))
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Multivariate SMK parallel fit --------------------------------------------------------

# ## Libraries
# library(mcmc)
# library(MASS)
# library(KernSmooth)
# library(fields)
# library(pscl)
# library(spBayes)
# library(mvtnorm)
# library(MCMCpack)
# library(Mposterior)
# library(parallel)
# library(doParallel)
# library(foreach)

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
per.core <- floor(n/n.core)
## Number of observations in different subsets
n.part <- c(rep(per.core,n.core-1),n.sample-per.core*(n.core-1))
## This is same as n.core
n.split <- length(n.part)
## Number of MCMC iterations
mcmc.sample <- 500
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
y <- dat.nomiss[,3:4]
x <- dat.nomiss[,-(1:4)]
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
  Y.part[[i]] <- Y.train[index.part[[i]],]
  ## Predictor in i th subset
  X.part[[i]] <- X.train[index.part[[i]],]
  ## Coordinates in i th subset
  coords.part[[i]] <- coords.train[index.part[[i]],]
  a <- setdiff(a,index.part[[i]])
}


##Partitioned GP function
##Works with subset i, for i=1,...,n.core


partitioned_GP <- function(i){
  
  ## library
  library(spBayes)
  
  ## Response in subset i
  ZZ <- Y.part[[i]]
  ## Predictor in subset i
  XX <- X.part[[i]]
  ## Coordinates in subset i
  CC <- coords.part[[i]]
  
  ## priors
  q <- ncol(ZZ)
  A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
  
  starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
  tuning <- list("phi"=rep(1,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.01,q))
  priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
                 "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))
  cov.model <- "exponential"
  
  ## Model fitting
  # tictoc::tic()
  m.1 <- spMvLM(list(ZZ[,1]~XX-1, ZZ[,2]~XX-1),
                coords=CC, starting=starting, tuning=tuning, priors=priors,
                n.samples=mcmc.sample, cov.model=cov.model, verbose = F)
  # tictoc::toc()
  
  
  ## Recover all MCMC samples in each subset
  # tictoc::tic()
  m.1.samp <- spRecover(m.1, start=n.burn+1, verbose=F)
  # tictoc::toc()
  
  ## List of MCMC iterates from subset i
  subAtom <- cbind(m.1.samp$p.beta.recover.samples,m.1.samp$p.theta.recover.samples)
  W.mcmc <- m.1.samp$p.w.recover.samples
  ## Delete this quantity from the memory
  rm(m.1.samp)
  ## Garbage cleaning
  gc()
  
  ## prediction covariates
  XXp <- dat.miss[,-c(1:4)]
  pred.covars <- mkMvX(list(XXp, XXp))
  
  # ## Prediction
  # tictoc::tic()
  Y.pred <- t(spPredict(m.1, coords.pred, pred.covars, start=0.5*mcmc.sample,verbose=T)$p.y.predictive.samples)
  # tictoc::toc()
  
  # Y.pred <- NULL
  # Y.pred[, 1] <- Y.predQ[,seq(1,ncol(Y.predQ),q)]
  # Y.pred[, 2] <- Y.predQ[,seq(2,ncol(Y.predQ),q)]
  
  # tic()
  # ## Prediction
  # Y.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[,-(1:3)]),
  #                       pred.coords=coords.pred[,],
  #                       start=0.5*mcmc.sample)$p.y.predictive.samples)
  # toc()
  
  # tic()
  
  ## Prediction
  # Yu.pred <- t(spPredict(m.1, pred.covars=as.matrix(dat.miss[-(1:m),-(1:3)]),
  #                        pred.coords=coords.pred[-(1:m),],
  #                        start=0.5*mcmc.sample)$p.y.predictive.samples)
  # Y.pred <- NULL
  # for(l in 1:(length(p.sub)-1)){
  #   m.1.pred <- spPredict(m.1, pred.covars=as.matrix(dat.miss[(p.sub[l]+1):p.sub[l+1], -(1:3)]),
  #                         pred.coords=coords.pred[(p.sub[l]+1):p.sub[l+1],],
  #                         start=0.5*mcmc.sample)$p.y.predictive.samples
  #   Y.pred <- cbind(Y.pred,t(m.1.pred))
  # }
  #
  # Y.pred <- cbind(Y.pred, Yu.pred)
  
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

# MSMK Results collection ------------------------------------------------------------

# collecting predictions
Y.medQ <- Y.med
Y.med <- matrix(0, nrow(dat.miss), q)
Y.med[, 1] <- Y.medQ[seq(1,length(Y.medQ),q)]
Y.med[, 2] <- Y.medQ[seq(2,length(Y.medQ),q)]

Y.lower_qntQ <- Y.lower_qnt
Y.lower_qnt <- matrix(0, nrow(dat.miss), q)
Y.lower_qnt[, 1] <- Y.lower_qntQ[seq(1,length(Y.medQ),q)]
Y.lower_qnt[, 2] <- Y.lower_qntQ[seq(2,length(Y.medQ),q)]

Y.upper_qntQ <- Y.upper_qnt
Y.upper_qnt <- matrix(0, nrow(dat.miss), q)
Y.upper_qnt[, 1] <- Y.upper_qntQ[seq(1,length(Y.medQ),q)]
Y.upper_qnt[, 2] <- Y.upper_qntQ[seq(2,length(Y.medQ),q)]

# average coverage
coverage_smk <- c(mean(Y_u[,1] >= Y.lower_qnt[,1] & Y_u[,1] <= Y.upper_qnt[,1]),
                  mean(Y_u[,2] >= Y.lower_qnt[,2] & Y_u[,2] <= Y.upper_qnt[,2]))
cat("Empirical average MSMK coverage for Responses:", round(mean(coverage_smk), 3))
(CI_avlen_smk <- mean(Y.upper_qnt[,]-Y.lower_qnt[,]))

# Root Mean Square Prediction Error
(rmspe_smk <- sqrt( colMeans( (Y_u[-(1:m),] - Y.med[-(1:m),])^2 ) )); mean(rmspe_smk)


# MSMK Posterior inference -----------------------------------------------------

# collecting
J <- ncol(obj[[1]]$atoms)
nn <- nrow(obj[[1]]$atoms)
post_smp_smk <- matrix(0, nn, J)
for (j in 1:J) {
  post_smp_smk[,j] <- sapply(1:n.core, function(a)obj[[a]]$atoms[,j]) %*% wts
}
posterior_smk <- post_smp_smk[,c(1:7,10:11)]
colnames(posterior_smk) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]", "phi[1]", "phi[2]")

# fixing true parameter for plotting
true_par <- c(matrix(B), matrix(sigma2)[-3], rep(phi, 2))


# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = out$timings[1],
                   "Combination" = out$timings[2],
                   "Prediction" = out$timings[3],
                   "Total time" = sum(out$timings))

elapsed_smk <- c("Fitting" = final.time,
                 "Combination" = recomb.time,
                 "Prediction" = quantile.calculation.time,
                 "Total time" = Time)

cat("minutes elapsed for BPS fully model-based uncertainty quantification : \n"); round(elapsed_times, 2)
cat("minutes elapsed for MSMK fully model-based uncertainty quantification : \n"); round(elapsed_smk, 2)


# Plotting data --------------------------------------------------

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(bayesplot)

mba_to_df <- function(coords, values, no.X = 500, no.Y = 500, h = 12) {
  surf <- MBA::mba.surf(cbind(coords, values), no.X = no.X, no.Y = no.Y,
                        exten = TRUE, sp = FALSE, h = h)$xyz.est
  data.frame(
    x = rep(surf$x, times = length(surf$y)),
    y = rep(surf$y, each  = length(surf$x)),
    z = as.vector(surf$z)
  )
}

df_Y1_true <- mba_to_df(crd_s, Y[, 1])
df_Y2_true <- mba_to_df(crd_s, Y[, 2])
df_Y1_bps  <- mba_to_df(crd_u, post_mean_Y[, 1])
df_Y2_bps  <- mba_to_df(crd_u, post_mean_Y[, 2])
df_Y1_smk  <- mba_to_df(crd_u, Y.med[, 1])
df_Y2_smk  <- mba_to_df(crd_u, Y.med[, 2])

lims_Y1 <- range(c(df_Y1_true$z, df_Y1_bps$z, df_Y1_smk$z), na.rm = TRUE)
lims_Y2 <- range(c(df_Y2_true$z, df_Y2_bps$z, df_Y2_smk$z), na.rm = TRUE)

rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu")[11:1])(100))

surf_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text         = element_text(size = 16, face = "bold"),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5, vjust = -0.35),
    plot.subtitle     = element_text(size = 14, hjust = 0.5, color = "grey30"),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 14),
    legend.key.height = unit(1.2, "cm")
  )

plot_surf <- function(df, title, subtitle = NULL, limits) {
  ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(colors = rdbu_pal, limits = limits,
                         name = "Value", na.value = NA) +
    coord_equal() +
    labs(x = "Easting", y = "Northing", title = title, subtitle = subtitle) +
    surf_theme
}

p_y1_true <- plot_surf(df_Y1_true, title = expression(Y[1] ~ "– True"),                                          limits = lims_Y1)
p_y1_bps  <- plot_surf(df_Y1_bps,  title = bquote(Y[1] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[1],  3)) ~ ")"),   limits = lims_Y1)
p_y1_smk  <- plot_surf(df_Y1_smk,  title = bquote(Y[1] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[1], 3)) ~ ")"),  limits = lims_Y1)
p_y2_true <- plot_surf(df_Y2_true, title = expression(Y[2] ~ "– True"),                                          limits = lims_Y2)
p_y2_bps  <- plot_surf(df_Y2_bps,  title = bquote(Y[2] ~ "– DBPS ( RMSPE:" ~ .(round(rmspe_Y[2],  3)) ~ ")"),   limits = lims_Y2)
p_y2_smk  <- plot_surf(df_Y2_smk,  title = bquote(Y[2] ~ "– MSMK ( RMSPE:" ~ .(round(rmspe_smk[2], 3)) ~ ")"),  limits = lims_Y2)

row1 <- (p_y1_true | p_y1_bps | p_y1_smk) + plot_layout(guides = "collect")
row2 <- (p_y2_true | p_y2_bps | p_y2_smk) + plot_layout(guides = "collect")
final_surf <- row1 / row2

png("surface_M_10_1000.png", width = 16, height = 10, units = "in", res = 300, family = "sans")
print(final_surf)
dev.off()

uc_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text     = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 18, face = "bold", color = "grey30")
  )

make_uc <- function(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord, ylab, subtitle) {
  df <- data.frame(x_ax, Yu_ord, CI_lower, CI_upper, Ymap_ord)
  ggplot(df, aes(x = x_ax, y = Yu_ord)) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  width = 0, linetype = "dashed", linewidth = 0.3, color = "#D41159") +
    geom_point(pch = 18, size = 3.5, color = "#1A85FF") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, color = "#D41159") +
    ylim(range(c(CI_lower, CI_upper))) +
    labs(x = "Ordered locations", y = ylab, subtitle = subtitle) +
    uc_theme
}

ord1 <- order(Y_u[-(1:m), 1])
ord2 <- order(Y_u[-(1:m), 2])

uc_Y1bps <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  post_qnt_Y[1,-(1:m),1][ord1], post_qnt_Y[2,-(1:m),1][ord1],
  post_mean_Y[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[1], 3),")")
)
uc_Y2bps <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  post_qnt_Y[1,-(1:m),2][ord2], post_qnt_Y[2,-(1:m),2][ord2],
  post_mean_Y[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("DBPS ( Empirical coverage: ", round(coverage_Y[2], 3),")")
)
uc_Y1smk <- make_uc(
  1:u, Y_u[-(1:m),1][ord1],
  Y.lower_qnt[-(1:m),1][ord1], Y.upper_qnt[-(1:m),1][ord1],
  Y.med[-(1:m),1][ord1],
  ylab = expression(Y[1]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[1], 3),")")
)
uc_Y2smk <- make_uc(
  1:u, Y_u[-(1:m),2][ord2],
  Y.lower_qnt[-(1:m),2][ord2], Y.upper_qnt[-(1:m),2][ord2],
  Y.med[-(1:m),2][ord2],
  ylab = expression(Y[2]),
  subtitle = paste("MSMK ( Empirical coverage: ", round(coverage_smk[2], 3),")")
)

png("UC_M_10_1000.png", width = 16, height = 10,
    units = "in", res = 300, family = "sans")
print((uc_Y1bps | uc_Y2bps) / (uc_Y1smk | uc_Y2smk))
dev.off()

bayesplot_theme_set(theme_classic(base_size = 16))
color_scheme_set("viridis")

ci_labels <- c(
  expression(beta[list(1,1)]), expression(beta[list(2,1)]),
  expression(beta[list(1,2)]), expression(beta[list(2,2)]),
  expression(phi[1]),          expression(phi[2]),
  expression(Sigma[list(1,1)]),expression(Sigma[list(2,1)]),expression(Sigma[list(2,2)])
)

ci_theme_extra <- theme(
  axis.text.x     = element_text(size = 18, face = "bold", angle = 45, hjust = 1),
  axis.text.y     = element_text(size = 16, face = "bold"),
  axis.title      = element_text(size = 16, face = "bold"),
  legend.text     = element_text(size = 14, face = "bold"),
  legend.title    = element_blank(),
  legend.position = "right",
  plot.title      = element_blank(),
  plot.subtitle   = element_text(size = 18, face = "bold", color = "grey30")
)

post_int_bps <- mcmc_recover_intervals(
  posterior_bps, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "DBPS (posterior means and 95% CI)") +
  ci_theme_extra

post_int_smk <- mcmc_recover_intervals(
  posterior_smk, true_par,
  prob = 0.95, prob_outer = 0.95, point_est = "mean", size = 4, alpha = 0.75
) +
  scale_x_discrete(labels = ci_labels) +
  labs(subtitle = "MSMK (posterior means and 95% CI)") +
  ci_theme_extra

png("CIpost_M_10_1000.png", width = 14, height = 5, units = "in", res = 300, family = "sans")
print((post_int_bps | post_int_smk) + plot_layout(guides = "collect"))
dev.off()


# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("BPS" = list("time"    = elapsed_times,
                             "metrics" = c("RMSPE" = rmspe_Y, "CI_len" = CI_avlen_bps, "ECoverage" = coverage_Y)),
                "SMK" = list("time"    = elapsed_smk,
                             "metrics" = c("RMSPE" = rmspe_smk, "CI_len" = CI_avlen_smk, "ECoverage" = coverage_smk)))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "simulation_multivariate_10_1000.RData")

#####################################################################################################################################################
