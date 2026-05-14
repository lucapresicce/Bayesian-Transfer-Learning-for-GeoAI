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
u <- 500
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
crd_u <- crd[-(1:n), ]
X_u <- X_or[-(1:n), ]
W_u <- W_or[-(1:n), ]
Y_u <- Y_or[-(1:n), ]

# Subset posterior models -----------------------------------------------------

# hyperparameters values
alfa_seq <- c(0.7, 0.8, 0.9)
phi_seq <- c(3, 4, 5)
hyperpar <- list(alpha = alfa_seq, phi = phi_seq)

# number of available computational cores
n.core <- parallel::detectCores(logical = F)-1


# Subset sensitivity ------------------------------------------------------

res0 <- c(25, 50, 100, 250, 500, 1000, 1250) # subset sizes
res1 <- NULL                                 # result container for rmspe
res2 <- NULL                                 # result container for fit time

pb <- txtProgressBar(0, length(res0), style = 3)
for (i in 1:length(res0)) {
  
  subset_size <- res0[i]
  
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
  
  # Results collection ----------------------------------------------------------
  
  # statistics computations Y
  pred_mat_Y  <- out$predictive$Yu
  post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
  
  # Root Mean Square Prediction Error
  rmspe_Y <- sqrt( colMeans( (Y_u - post_mean_Y)^2 ) )
  
  # Save result -----------------------------------------------------------------
  
  res1 <- c(res1, mean(rmspe_Y))
  res2 <- c(res2, sum(out$timings))
  
  setTxtProgressBar(pb = pb, value = i)
  
}

# Min-Max normalization (to better comparison)
res1_norm <- (res1-min(res1))/(max(res1)-min(res1))
res2_norm <- (res2-min(res2))/(max(res2)-min(res2))

df <- data.frame(
  x     = res0,
  rmspe = res1_norm,
  time  = res2_norm
)

col_rmspe <- "cornflowerblue"
col_time  <- "seagreen"

p_sens <- ggplot(df, aes(x = x)) +
  geom_line(aes(y = rmspe, color = "RMSPE",       linetype = "RMSPE"),
            linewidth = 1, key_glyph = "path") +
  geom_line(aes(y = time,  color = "Time", linetype = "Time"),
            linewidth = 1, key_glyph = "path") +
  scale_color_manual(name = NULL,
                     values = c("RMSPE"       = col_rmspe,
                                "Time" = col_time)) +
  scale_linetype_manual(name = NULL,
                        values = c("RMSPE"       = "solid",
                                   "Time" = "dashed")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Subset size", y = "Normalized value") +
  theme_classic(base_size = 14) +
  theme(
    axis.text        = element_text(size = 16, face = "bold"),
    axis.title       = element_text(size = 16, face = "bold"),
    legend.position  = "right",
    legend.text      = element_text(size = 13),
    legend.key.width = unit(1.25, "cm")
  )

png("outout/subset_sens.png", width = 10, height = 5, units = "in", res = 300, family = "sans")
print(p_sens)
dev.off()

