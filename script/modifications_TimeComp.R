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
library(bayesplot)
####
library(MASS)
library(BART)
library(spNNGP)
library(spBayes)
library(subart)


# Data generation - 5k ----------------------------------------------------

# dimensions
n <- 5000
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
crd_u <- crd[-(1:n), ]
X_u <- X_or[-(1:n), ]
W_u <- W_or[-(1:n), ]
Y_u <- Y_or[-(1:n), ]

# mvGP --------------------------------------------------------------------

# Arrange data
X_train <- mkMvX(list(X, X))
Y_train_1 <- Y[seq(1, length(Y), q)]
Y_train_2 <- Y[seq(2, length(Y), q)]
coords <- crd_s

# Priors
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
starting <- list("phi" = rep(4, q),
                 "A" = A.starting,
                 "Psi" = rep(1, q)
                 )
tuning <- list("phi" = rep(1, q),
               "A" = rep(0.01, length(A.starting)),
               "Psi" = rep(0.01, q)
               )
priors <- list("beta.Flat",
               "phi.Unif" = list( rep(4/10, q), rep(4/0.1, q)),
               "K.IW" = list(q+1, diag(0.1, q)),
               "Psi.ig"= list(c(2, 2), c(0.1, 0.1))
               )

tic("total")
tic("fit")
# Fit
n.samples <- 2000
burn.in <- 0.5 * n.samples
mvgp_fit <- spMvLM(list(Y_train_1 ~ X - 1, Y_train_2 ~ X - 1),
              coords = coords, starting = starting, tuning = tuning, priors = priors,
              n.samples = n.samples, cov.model = "exponential", n.report = 100)
fit_time <- toc()

# New locations and covariates
coords.u <- crd_u
X_test <- mkMvX(X = list(X_u, X_u))
Y_u_1 <- Y_u[seq(1, length(Y_u), q)]
Y_u_2 <- Y_u[seq(2, length(Y_u), q)]

tic("prediction")
# Prediction
mvgp_preds <- spPredict(mvgp_fit, pred.covars = X_test, pred.coords = coords.u,
                      start = burn.in)
mvgp_preds_mean <- apply(mvgp_preds$p.y.predictive.samples, 1, mean)
mvgp_preds_mean_1 <- mvgp_preds_mean[seq(1, length(mvgp_preds_mean), q)]
mvgp_preds_mean_2 <- mvgp_preds_mean[seq(2, length(mvgp_preds_mean), q)]
prd_time <- toc()

# timing
tot_time <- toc()

# Evaluation
mvgp_time <- tot_time
(mvgp_rmspe_1 <- sqrt( mean( (mvgp_preds_mean_1 - Y_u[, 1])^2 ) ))
(mvgp_rmspe_2 <- sqrt( mean( (mvgp_preds_mean_2 - Y_u[, 2])^2 ) ))
(mvgp_cor_1 <- cor(Y_u[, 1], mvgp_preds_mean_1))
(mvgp_cor_2 <- cor(Y_u[, 2], mvgp_preds_mean_2))

mvgp_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                "Total time" = as.numeric(tot_time$toc-tot_time$tic))

mvgp_rmspe <- c(mvgp_rmspe_1, mvgp_rmspe_2)
mvgp_corr <- c(mvgp_cor_1, mvgp_cor_2)

# suBART ------------------------------------------------------------------

# devtools::install_github("MateusMaiaDS/subart")

library(subart)

tic()
# Fit
subart_fit <- subart(
  x_train = as.data.frame(X[,-1]),
  y_mat = Y,
  x_test = as.data.frame(X_u[,-1]),
  n_tree = 100,
  node_min_size = 5,
  n_mcmc = 2000,
  n_burn = 1000,
  varimportance = F
)
fit_time <- toc()

# Posterior predictive samples on new data
subart_preds <- subart_fit$y_hat_test_mean

# Evaluation
subart_time <- fit_time
(subart_rmspe <- sqrt( colMeans( (Y_u - subart_preds)^2 ) )); mean(subart_rmspe)
(subart_cor_1 <- cor(Y_u[, 1], subart_preds[, 1]))
(subart_cor_2 <- cor(Y_u[, 2], subart_preds[, 2]))
subart_corr <- c(subart_cor_1, subart_cor_2)

# save results ------------------------------------------------------------

gc()
# Save the entire environment
competitor_results_5k_M <- list("mvgp"   = list("time"  = mvgp_times,
                                                "rmspe" = mvgp_rmspe,
                                                "corr"  = mvgp_corr),
                                "subart" = list("time"  = subart_time,
                                                "rmspe" = subart_rmspe,
                                                "corr"  = subart_corr))

rm(list = ls()[which(!(ls() %in% c("competitor_results_5k_M")))])
save.image(file = "review_TimeComp_5k_M.RData")

# to read results
# load("review_TimeComp_5k_M.RData")
# competitor_results_5k_M$mvgp$time
# competitor_results_5k_M$mvgp$rmspe
# competitor_results_5k_M$mvgp$corr
# competitor_results_5k_M$subart$time
# competitor_results_5k_M$subart$rmspe
# competitor_results_5k_M$subart$corr


#################################################################
# Data generation - 10k ---------------------------------------------------

# dimensions
n <- 10000
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
crd_u <- crd[-(1:n), ]
X_u <- X_or[-(1:n), ]
W_u <- W_or[-(1:n), ]
Y_u <- Y_or[-(1:n), ]

# mvGP --------------------------------------------------------------------


# DO NOT RUN - INFEASIBLE
# # Arrange data
# X_train <- mkMvX(list(X, X))
# Y_train_1 <- Y[seq(1, length(Y), q)]
# Y_train_2 <- Y[seq(2, length(Y), q)]
# coords <- crd_s
# 
# # Priors
# A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
# starting <- list("phi" = rep(4, q),
#                  "A" = A.starting,
#                  "Psi" = rep(1, q)
#                  )
# tuning <- list("phi" = rep(1, q),
#                "A" = rep(0.01, length(A.starting)),
#                "Psi" = rep(0.01, q)
#                )
# priors <- list("beta.Flat",
#                "phi.Unif" = list( rep(4/10, q), rep(4/0.1, q)),
#                "K.IW" = list(q+1, diag(0.1, q)),
#                "Psi.ig"= list(c(2, 2), c(0.1, 0.1))
#                )
# 
# tic("total")
# tic("fit")
# # Fit
# n.samples <- 2000
# burn.in <- 0.5 * n.samples
# mvgp_fit <- spMvLM(list(Y_train_1 ~ X - 1, Y_train_2 ~ X - 1),
#               coords = coords, starting = starting, tuning = tuning, priors = priors,
#               n.samples = n.samples, cov.model = "exponential", n.report = 100)
# fit_time <- toc()
# 
# # New locations and covariates
# coords.u <- crd_u
# X_test <- mkMvX(X = list(X_u, X_u))
# Y_u_1 <- Y_u[seq(1, length(Y_u), q)]
# Y_u_2 <- Y_u[seq(2, length(Y_u), q)]
# 
# tic("prediction")
# # Prediction
# mvgp_preds <- spPredict(mvgp_fit, pred.covars = X_test, pred.coords = coords.u,
#                       start = burn.in)
# mvgp_preds_mean <- apply(mvgp_preds$p.y.predictive.samples, 1, mean)
# mvgp_preds_mean_1 <- mvgp_preds_mean[seq(1, length(mvgp_preds_mean), q)]
# mvgp_preds_mean_2 <- mvgp_preds_mean[seq(2, length(mvgp_preds_mean), q)]
# prd_time <- toc()
# 
# # timing
# tot_time <- toc()
# 
# # Evaluation
# mvgp_time <- tot_time
# (mvgp_rmspe_1 <- sqrt( mean( (mvgp_preds_mean_1 - Y_u[, 1])^2 ) ))
# (mvgp_rmspe_2 <- sqrt( mean( (mvgp_preds_mean_2 - Y_u[, 2])^2 ) ))
# (mvgp_cor_1 <- cor(Y_u[, 1], mvgp_preds_mean_1))
# (mvgp_cor_2 <- cor(Y_u[, 2], mvgp_preds_mean_2))
# 
# # par(mfrow = c(1, 2))
# # plot(Y_u[, 1], mvgp_preds_mean_1)
# # abline(0, 1)
# # plot(Y_u[, 2], mvgp_preds_mean_2)
# # abline(0, 1)
# # par(mfrow = c(1, 1))
# 
# mvgp_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
#                 "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
#                 "Total time" = as.numeric(tot_time$toc-tot_time$tic))
# 
# mvgp_rmspe <- c(mvgp_rmspe_1, mvgp_rmspe_2)
# mvgp_corr <- c(mvgp_cor_1, mvgp_cor_2)

mvgp_times <- NULL
mvgp_rmspe <- NULL
mvgp_corr  <- NULL

# suBART ------------------------------------------------------------------

# devtools::install_github("MateusMaiaDS/subart")

library(subart)

tic()
# Fit
subart_fit <- subart(
  x_train = as.data.frame(X[,-1]),
  y_mat = Y,
  x_test = as.data.frame(X_u[,-1]),
  n_tree = 100,
  node_min_size = 5,
  n_mcmc = 2000,
  n_burn = 1000,
  varimportance = F
)
fit_time <- toc()

# Posterior predictive samples on new data
subart_preds <- subart_fit$y_hat_test_mean

# Evaluation
subart_time <- fit_time
(subart_rmspe <- sqrt( colMeans( (Y_u - subart_preds)^2 ) )); mean(subart_rmspe)
(subart_cor_1 <- cor(Y_u[, 1], subart_preds[, 1]))
(subart_cor_2 <- cor(Y_u[, 2], subart_preds[, 2]))
subart_corr <- c(subart_cor_1, subart_cor_2)


# save results ------------------------------------------------------------

gc()
# Save the entire environment
competitor_results_10k_M <- list("mvgp"   = list("time"  = mvgp_times,
                                                 "rmspe" = mvgp_rmspe,
                                                 "corr"  = mvgp_corr),
                                 "subart" = list("time"  = subart_time,
                                                 "rmspe" = subart_rmspe,
                                                 "corr"  = subart_corr))

rm(list = ls()[which(!(ls() %in% c("competitor_results_10k_M")))])
save.image(file = "review_TimeComp_10k_M.RData")

# to read results
# load("review_TimeComp_10k_M.RData")
# competitor_results_10k_M$mvgp$time
# competitor_results_10k_M$mvgp$rmspe
# competitor_results_10k_M$mvgp$corr
# competitor_results_10k_M$subart$time
# competitor_results_10k_M$subart$rmspe
# competitor_results_10k_M$subart$corr
