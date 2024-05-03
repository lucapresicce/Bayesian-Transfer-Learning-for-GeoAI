#####################################################################################################################################################
## MASMK DATA ANALYSIS - multivariate model #########################################
rm(list = ls())
gc()
setwd(".../reproducibility-code/analyses/multivariate")

# Packages --------------------------------------------------------------------
library(ASMK)
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
library(ggplot2)
library(sf)
library(geoR)
library(spBayes)
library(bayesplot)
library(raster)
library(corrplot)


# Data loading ----------------------------------------------------------------

# load preprocessed RData (unzip before execute)
load("cleaned_data2_expanded.RData")
full_data <- data_cleaned2
rm(list = c("data_cleaned2"))
names(full_data)

# take a look to data structure and responses variables
head(full_data)
dim(full_data)

# train data and test data
set.seed(1997)
test_ind <- sample.int(nrow(full_data), floor(0.025 * nrow(full_data)) )
train_data <- full_data[-test_ind, ]
test_data<- full_data[test_ind, ]

# select responses and predictors sets
response_set <- c("NDVI", "red reflectance")
q <- length(response_set)
predictor_set <- c("Non_Vegetated_or_Builtup_Lands")
p <- length(predictor_set)+1

# define train dimensions
N <- nrow(train_data)
crd_S <- matrix(as.matrix(train_data[, c("scaled_x","scaled_y")]), ncol = 2)
Y_S   <- matrix(as.matrix(train_data[, response_set]), ncol = q)
X_S   <- cbind(1, matrix(as.matrix(train_data[, predictor_set]), ncol = (p-1)))

# define test dimensions
U <- nrow(test_data)
crd_U <- matrix(as.matrix(test_data[, c("scaled_x","scaled_y")]), ncol = 2)
Y_U   <- matrix(as.matrix(test_data[, response_set]), ncol = q)
X_U   <- cbind(1, matrix(as.matrix(test_data[, predictor_set]), ncol = (p-1)))

# remove full dataset and free memory
rm(list = c("full_data", "train_data", "test_data"))
gc()


# EDA -------------------------------------------------------------------------

# return to the original scale
Y_S[,1] <- (exp(Y_S[,1])-1)

# linear model to collect residual
lin_res <- Y_S - X_S %*% solve(crossprod(X_S))%*%(crossprod(X_S, Y_S))
summary(lin_res); cov(lin_res)

# subsample for feasible EDA
set.seed(1997)
eda_ind <- sample.int(N, round(N*0.001))

# computing the maximum distance
d.max <- sqrt((max(crd_S[,1]) - min(crd_S[,1]))^2 +
                (max(crd_S[,2]) - min(crd_S[,2]))^2)
d.max # around 1.572

# check the variogram for the first response
v.res_1 <- variog(coords = crd_S[eda_ind, ], data = lin_res[eda_ind, 1],
                  uvec = (seq(0, 0.675, length = 30))) # 30

par(mfrow=c(1,1))
vario.fit_1 <- variofit(v.res_1, cov.model="exponential")
summary(vario.fit_1)

variofitphi.resid1 <- 1 / vario.fit_1$cov.pars[2]; variofitphi.resid1
variofitalpha.resid1 <- vario.fit_1$cov.pars[1] / (vario.fit_1$nugget+vario.fit_1$cov.pars[1]); variofitalpha.resid1

# check the variogram for the second response
v.res_2 <- variog(coords = crd_S[eda_ind, ], data = lin_res[eda_ind, 2],
                  uvec = (seq(0, 0.5, length = 30))) # 30

par(mfrow=c(1,1))
vario.fit_2 <- variofit(v.res_2, cov.model="exponential")
summary(vario.fit_2)

variofitphi.resid2 <- 1 / vario.fit_2$cov.pars[2]; variofitphi.resid2
variofitalpha.resid2 <- vario.fit_2$cov.pars[1] / (vario.fit_2$nugget+vario.fit_2$cov.pars[1]); variofitalpha.resid2
# # variofitalpha.resid2 <- vario.fit_2$cov.pars[1] / (v.res_2$v[1]-1e-4+vario.fit_2$cov.pars[1]); variofitalpha.resid2

# check the variogram for the first response
v.res_y1 <- variog(coords = crd_S[eda_ind, ], data = Y_S[eda_ind, 1],
                   # uvec = (seq(0, 3*d.max, length = 50)[-c(7:11)])) # 30
                   uvec = (seq(0, 3*d.max, length = 50)[c(1:6,12:13)])) # 30

par(mfrow=c(1,1))
vario.fit_y1 <- variofit(v.res_y1, cov.model="exponential")
summary(vario.fit_y1)

# check the variogram for the second response
v.res_y2 <- variog(coords = crd_S[eda_ind, ], data = Y_S[eda_ind, 2],
                   # uvec = (seq(0, 3*d.max, length = 50)[-c(6:11)])) # 30
                   uvec = (seq(0, 3*d.max, length = 50)[c(1:5, 12:13)])) # 30

par(mfrow=c(1,1))
vario.fit_y2 <- variofit(v.res_y2, cov.model="exponential")
summary(vario.fit_y2)

# free memory
rm(list = c("lin_res"))
gc()

# SubSubsample for model testing ----------------------------------------------

set.seed(1997)
n <- 1000000
subsample <- sample(1:N, n)
crd_s <- crd_S[subsample, ]
y     <- Y_S[subsample, ]
x     <- X_S[subsample, ]

u <- 2500
subsampleu <- sample(1:U, u)
crd_u <- crd_U[subsampleu, ]
y_u   <- Y_U[subsampleu, ]
x_u   <- X_U[subsampleu, ]


# Fit linear model----------------------------------------------------

# dimension
n <- nrow(x)
p <- ncol(x)

naive <- x_u[,1] %*% solve(crossprod(x[,1]))%*%(crossprod(x[,1], y))

(rmspe_naive <- sqrt(colMeans((y_u - naive)^2))); mean(rmspe_naive)
(mape_naive <- colMeans( abs(y_u - naive) ) ); mean(mape_naive)

#####################################################################################################################################################

# Subset posterior models -----------------------------------------------------

# hyperparameters values by looking at variograms
(alfa_seq <- sort(c(variofitalpha.resid1, variofitalpha.resid2)))
(phi_seq <- sort(2*c(variofitphi.resid1, variofitphi.resid2)))

# function for the fit loop
fit_loop <- function(i) {
  
  Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
  p <- ncol(Xi); q <- ncol(Yi)
  bps <- ASMK::BPS_weights_MvT(data = list(Y = Yi, X = Xi),
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
  result <- ASMK::BPS_post_MvT(data = list(Y = Ys, X = Xs), coords = crds,
                               X_u = x_u, crd_u = crd_u,
                               priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                             V_r = diag(10, p),
                                             Psi = diag(1, q),
                                             nu = 3),
                               hyperpar = list(alpha = alfa_seq, phi = phi_seq),
                               W = Ws, R = 1)
  
  return(result)
}

# subsetting data
subset_size <- 500
K <- n/subset_size
data_part <- subset_data(data = list(Y = y, X = x, crd = crd_s), K = K)

# ASMK parallel fit -------------------------------------------------------

# number of clusters for parallel implementation
n.core <- parallel::detectCores(logical = F)-1

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
comb_bps <- BPS_PseudoBMA(obj_fit)
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
pred_mat_W <- sapply(1:R, function(r){predictions[[r]]$Pred[[1]]$Wu}, simplify = "array")
post_mean_W <- apply(pred_mat_W, c(1,2), mean)
post_var_W <- apply(pred_mat_W, c(1,2), sd)
post_qnt_W <- apply(pred_mat_W, c(1,2), quantile, c(0.025, 0.975))

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]]$Pred[[1]]$Yu}, simplify = "array")
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(y_u[,1] >= post_qnt_Y[1,,1] & y_u[,1] <= post_qnt_Y[2,,1]),
                mean(y_u[,2] >= post_qnt_Y[1,,2] & y_u[,2] <= post_qnt_Y[2,,2]))
# mean(Y_u >= post_qnt_Y[1,,] & Y_u <= post_qnt_Y[2,,])
cat("Empirical average coverage for Spatial process:", round(mean(coverage_Y), 3))
(CI_avlen_masmk <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
(rmspe_Y <- sqrt( colMeans( (y_u - post_mean_Y)^2 ) )); mean(rmspe_Y)
(mape_Y <- colMeans( abs(y_u - post_mean_Y) ) ); mean(mape_Y)


# Posterior inference -----------------------------------------------------

beta_smp <- sapply(1:R, function(r){predictions[[r]]$Post[[1]]$beta[1:p,]}, simplify = "array")
post_mean_beta <-  apply(beta_smp, c(1,2), mean)
post_var_beta <- apply(beta_smp, c(1,2), sd)
post_qnt_beta <- apply(beta_smp, c(1,2), quantile, c(0.05, 0.95))

sigma_smp <- sapply(1:R, function(r){predictions[[r]]$Post[[1]]$sigma}, simplify = "array")
post_mean_sigma <- apply(sigma_smp, c(1,2), mean)
post_var_sigma <- apply(sigma_smp, c(1,2), sd)
post_qnt_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.05, 0.95))

(post_mean_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps)
post_var_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                   "Combination" = as.numeric(comb_time$toc-comb_time$tic),
                   "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                   "Total time" = as.numeric(tot_time$toc-tot_time$tic))

cat("minutes elapsed for fully model-based uncertainty quantification : \n"); round(elapsed_times/60, 2)


# Plotting data --------------------------------------------------

# subsetting for plotting interpolation
set.seed(1997)
plot_ind <- sample(1:N, 5000)

# interpolation
h <- 12
surf.raw.NDVI <- mba.surf(cbind(crd_S[plot_ind,], Y_S[plot_ind, 1]), no.X = 300,
                          no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

# Color palettes
col.pal1 <- colorRampPalette(RColorBrewer::brewer.pal(11, 'RdBu')[1:11])
colors1 <- col.pal1(5)

# plot limits
xlim <- range(crd_S[plot_ind, 1])
zlimN <- range(surf.raw.NDVI[["z"]][which(!is.na(surf.raw.NDVI[["z"]]))])

# ggplot version plot
surf_df <- as.data.frame(surf.raw.NDVI)
(ggNDVI <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = colors1, limits = zlimN) +
    labs(x = "Easting", y = "Northing") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

# interpolation
h <- 12
surf.raw.RR <- mba.surf(cbind(crd_S[plot_ind,], Y_S[plot_ind, 2]), no.X = 300,
                        no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est
# Color palettes
col.pal2 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGn')[1:9])
colors2 <- rev(col.pal2(5))

# plot limits
xlim <- range(crd_S[plot_ind, 1])
zlimR <- range(surf.raw.RR[["z"]][which(!is.na(surf.raw.RR[["z"]]))])

# ggplot version plot
surf_df <- as.data.frame(surf.raw.RR)
(ggRR <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = colors2, limits = zlimR) +
    labs(x = "Easting", y = "Northing") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

# Plotting test data --------------------------------------------------

# interpolation
h <- 12
surf.raw.NDVIu <- mba.surf(cbind(crd_u, y_u[, 1]), no.X = 300,
                           no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

# ggplot version plot
surf_df <- as.data.frame(surf.raw.NDVIu)
(ggNDVIu <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = colors1, limits = zlimN*.8) +
    labs(x = "Easting", y = "Northing") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

# interpolation
h <- 12
surf.raw.RRu <- mba.surf(cbind(crd_u, y_u[, 2]), no.X = 300,
                         no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

# ggplot version plot
surf_df <- as.data.frame(surf.raw.RRu)
(ggRRu <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = colors2, limits = zlimR) +
    labs(x = "Easting", y = "Northing") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

# Plotting results ------------------------------------------------------------

# interpolation
h <- 12
surf.raw.NDVIhat <- mba.surf(cbind(crd_u, post_mean_Y[, 1]), no.X = 300,
                             no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

# ggplot version plot
surf_df <- as.data.frame(surf.raw.NDVIhat)
(ggNDVIhat <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = colors1, limits = zlimN) +
    labs(x = "Easting", y = "Northing") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

# interpolation
h <- 12
surf.raw.RRhat <- mba.surf(cbind(crd_u, post_mean_Y[, 2]), no.X = 300,
                           no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

# ggplot version plot
surf_df <- as.data.frame(surf.raw.RRhat)
(ggRRhat <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = colors2, limits = zlimR) +
    labs(x = "Easting", y = "Northing") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

# graphical UC for Y1 (ordered)
ord_y <- order(y_u[,1])
plot_ind <- sample(1:u, 250)
df <- data.frame(
  x_ax = (1:u)[plot_ind],
  Yu_ord = y_u[ord_y, 1][plot_ind],
  CI_Y_lower = post_qnt_Y[1, ord_y, 1][plot_ind],
  CI_Y_upper = post_qnt_Y[2, ord_y, 1][plot_ind],
  Ymap_ord = post_mean_Y[ord_y, 1][plot_ind])

# Create the ggplot
(uc_Y1 <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
    geom_point(pch = 18, size = 3.5, col = "#1A85FF") +
    geom_errorbar(aes(ymin = CI_Y_lower, ymax = CI_Y_upper),
                  width = 0.05,
                  linetype = "dashed",
                  linewidth = 0.05,
                  color = "#D41159") +
    ylim(range(c(df$CI_Y_lower, df$CI_Y_upper))) +
    labs(x = "Ordered locations", y = "Response values") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black")))


# graphical UC for Y2 (ordered)
ord_y <- order(y_u[,2])
plot_ind <- sample(1:u, 250)
df <- data.frame(
  x_ax = (1:u)[plot_ind],
  Yu_ord = y_u[ord_y, 2][plot_ind],
  CI_Y_lower = post_qnt_Y[1, ord_y, 2][plot_ind],
  CI_Y_upper = post_qnt_Y[2, ord_y, 2][plot_ind],
  Ymap_ord = post_mean_Y[ord_y, 2][plot_ind])

# Create the ggplot
(uc_Y2 <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
    geom_point(pch = 18, size = 3.5, col = "#1A85FF") +
    geom_errorbar(aes(ymin = CI_Y_lower, ymax = CI_Y_upper),
                  width = 0.05,
                  linetype = "dashed",
                  linewidth = 0.05,
                  color = "#D41159") +
    ylim(range(c(df$CI_Y_lower, df$CI_Y_upper))) +
    labs(x = "Ordered locations", y = "Response values") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black")))


# one line
width <- 360*5
height <- 360
pointsize <- 12
png("dataanalysis_multivariate_pred.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(ggNDVIu, ggRRu, ggNDVIhat, ggRRhat, nrow = 1)
dev.off()

# all together
width <- 360*5
height <- 360*2
pointsize <- 12
png("dataanalysis_multivariate_all.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(ggNDVI, ggRR, ggNDVIu, ggRRu, uc_Y1, uc_Y2, ggNDVIhat, ggRRhat, nrow = 2, ncol = 4)
dev.off()

# all together

# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("time"    = elapsed_times,
                "fit"     = obj_fit,
                "comb"    = comb_bps,
                "pred"    = predictions,
                "post"    = list(beta_smp, sigma_smp),
                "metrics" = c("RMSPE" = rmspe_Y, "naiveR" = rmspe_naive, "CI_len" = CI_avlen_masmk, "ECoverage" = coverage_Y, "MAPE" = mape_Y, "naiveM" = mape_naive))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "datanalysis_multivariate.RData")
# load("datanalysis_multivariate.RData")
# results$time; cat("minutes elapsed for fully model-based uncertainty quantification : \n"); round(results$time/60, 2)
# results$metrics
# results$comb

