#####################################################################################################################################################
## Univariate DATA ANALYSIS - univariate model #########################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets/")

# Packages --------------------------------------------------------------------
library(spBPS)
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
library(mapproj)
library(ggplot2)
library(rworldmap)
library(sf)
library(geoR)
library(spBayes)
library(bayesplot)

# Data loading ----------------------------------------------------------------

# load preprocessed RData
load("data/SST_data_2022_06_21.RData")

# offset on coordinates for plotting
SSTdata$lon <- SSTdata$lon-180

# sinusoidally projected coordinates (scaled to 1000km units) as explanatory variables
knots.sinusoidal <- mapproject(SSTdata$lon, SSTdata$lat, 
                               projection = "sinusoidal")

radius.of.earth = 6.371            ## 6.371 * 1000 kilometers 
knots.sinusoidal = radius.of.earth * (cbind(knots.sinusoidal$x, 
                                            knots.sinusoidal$y))
SSTdata$projX <- knots.sinusoidal[, 1]
SSTdata$projY <- knots.sinusoidal[, 2]

# take a look to data structure and response variable
head(SSTdata)
dim(SSTdata)

# train data
set.seed(1997)
ds_ind <- sample.int(dim(SSTdata)[1], 5e6)
SST_ds <- SSTdata[ds_ind,]

# test data
set.seed(1997)
test_ind <- sample.int(dim(SST_ds)[1], floor(0.025 * length(SST_ds$sst)))
SST_train <- SST_ds[-test_ind, ]
SST_test <- SST_ds[test_ind, ]

# define train dimensions
N <- nrow(SST_train)
crd_S <- SST_train[,c("lon","lat")]
Y <- SST_train$sst
X <- SST_train[,c("projX","projY")]
p <- ncol(X)

# define test dimensions
U <- nrow(SST_test)
crd_U <- SST_test[,c("lon","lat")]
Y_U <- SST_test$sst
X_U <- SST_test[,c("projX","projY")]

# remove full dataset and free memory
rm(list = c("SSTdata", "knots.sinusoidal", "SST_ds"))
gc()

# EDA -------------------------------------------------------------------------

# subsample for feasible EDA
set.seed(1997)
subind <- sample.int(N, round(N*0.005))

# computing the maximum distance
d.max <- sqrt((max(SST_train$lon) - min(SST_train$lon))^2 + 
                (max(SST_train$lat) - min(SST_train$lat))^2)
d.max # around 391.4029 (multiply by 111.139) ~> 43,520 KM

# check the variogram 
v.resp <- variog(coords = crd_S[subind, ], data = Y[subind], 
                 uvec = (seq(0, 0.8*d.max, length = 50))) # 30

vario.fit <- variofit(v.resp, cov.model="exponential")
summary(vario.fit)
print(vario.fit)

# plotting the variogram
# plotting results
width <- 360*2
height <- 360
pointsize <- 12
png("output/eda_univariate.png", width = width, height = height, pointsize = pointsize, family = "sans")

plot(v.resp, pch = 19, frame.plot = FALSE, axes = F, ylim = c(0, max(v.resp$v)), xlim = c(0, vario.fit$max.dist))
box(bty="l")
axis(2)
axis(1) 
lines(vario.fit, col = 4, lwd = 2, lty = 2)
abline(h = vario.fit$nugget, col = 2, lwd = 3)
text(70, 12.5, 
     bquote(tau^2 ~ " = " ~ .(round(vario.fit$nugget, 2))),
     col = 2, cex = 2)
abline(h = vario.fit$cov.pars[1], col = "green4", lwd = 2, lty = 2)
text(25, 110, 
     bquote(sigma^2 ~ " = " ~ .(round(vario.fit$cov.pars[1], 2))), 
     col = "green4", cex = 1.5)
abline(v = vario.fit$practicalRange, col = "green4", lwd = 2, lty = 2)
text(120, 60, 
     bquote(rho[0] ~ " = " ~ .(round(vario.fit$practicalRange, 2))), 
     col = "green4", cex = 1.5)

dev.off()

# free memory
gc()

# SubSubsample for model fitting ----------------------------------------------

set.seed(1997)
n <- 1000000
subsample <- sample(1:N, n)
y <- Y[subsample]
x <- cbind(1, matrix(as.matrix(X[subsample, ]), ncol = 2))
crd_s <- matrix(as.matrix(crd_S[subsample, ]), ncol = 2)

u <- 2500
subsampleu <- sample(1:U, u)
y_u <- Y_U[subsampleu]
x_u <- cbind(1, matrix(as.matrix(X_U[subsampleu, ]), ncol = 2))
crd_u <- matrix(as.matrix(crd_U[subsampleu, ]), ncol = 2)


# Fit Bayesian linear model----------------------------------------------------

# dimension
n <- nrow(x)
p <- ncol(x)

set.seed(1234)
bLM <- bayesLMConjugate(y~x-1, n.samples = 20000, 
                        beta.prior.mean = rep(0, times = p),
                        beta.prior.precision = matrix(0, nrow=p, ncol=p),
                        prior.shape = 2, prior.rate = 2)

round(summary(bLM$p.beta.tauSq.samples)$statistics, 2)
round(summary(bLM$p.beta.tauSq.samples)$quantiles, 2)

# posterior predictive
bLM.pred <- spPredict(bLM, pred.covars = x_u, pred.coords = crd_u,
                      start = 10000)

y.bLM <- apply(bLM.pred$p.y.predictive.samples, 1, mean)

(RMSPE_bLM <- sqrt(mean((y_u - y.bLM)^2)))
(RMSPE_LM <- sqrt(mean((y_u - x_u %*% lm.fit(x = x,y = y)$coef)^2)))


#####################################################################################################################################################
# Subset posterior models -----------------------------------------------------

# chioce the hyperparameter values by looking at: summary(vario.fit)
delta_seq <- ifelse(vario.fit$nugget/vario.fit$cov.pars[1]==0, 1e-6, vario.fit$nugget/vario.fit$cov.pars[1])
phi_vario <- 1 / vario.fit$cov.pars[2]
phi_step <- (1 / vario.fit$cov.pars[2])/3
phi_seq <- seq(-1*phi_step, 1*phi_step, phi_step) + phi_vario

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
                            X_u = x_u, crd_u = crd_u,
                            priors = list(mu_b = matrix(rep(0, p)),
                                          V_b = diag(10, p),
                                          a = 2,
                                          b = 2),
                            hyperpar = list(delta = delta_seq, phi = phi_seq),
                            W = Ws, R = 1)
  
  return(result)
}

# Subsetting data ----------------------------------------------------------

subset_size <- 250
K <- n/subset_size
data_part <- subset_data(data = list(Y = matrix(y), X = x, crd = crd_s), K = K)

# spBPS parallel fit -------------------------------------------------------

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
pred_mat_W <- sapply(1:R, function(r){predictions[[r]][[1]]})
post_mean_W <- rowMeans(pred_mat_W)
post_var_W <- apply(pred_mat_W, 1, sd)
post_qnt_W <- apply(pred_mat_W, 1, quantile, c(0.025, 0.975))

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){predictions[[r]][[2]]})
post_mean_Y <- rowMeans(pred_mat_Y)
post_var_Y <- apply(pred_mat_Y, 1, sd)
post_qnt_Y <- apply(pred_mat_Y, 1, quantile, c(0.025, 0.975))

# Empirical coverage for Y
coverage_Y <- mean(y_u >= post_qnt_Y[1,] & y_u <= post_qnt_Y[2,])
cat("Empirical coverage for Response:", round(coverage_Y, 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,]-post_qnt_Y[1,]))

# Root Mean Square Prediction Error
(rmspe_Y <- sqrt( mean( (y_u - post_mean_Y)^2 ) ))

# Naive linear model RMSPE
naive <- x_u %*% lm.fit(y = y, x = x)$coef
(rmspe_naive <- sqrt( mean( (y_u - naive)^2 ) ))


# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                   "Combination" = as.numeric(comb_time$toc-comb_time$tic),
                   "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                   "Total time" = as.numeric(tot_time$toc-tot_time$tic))

cat("minutes elapsed for fully model-based uncertainty quantification : \n"); round(elapsed_times/60, 2)


# Posterior inference ---------------------------------------------------------

# collecting posterior sample
smp <- sapply(1:R, function(r){c(predictions[[r]][[4]], predictions[[r]][[3]][,1:p])})
post_mean_smp <- rowMeans(smp)
post_var_smp <- apply(smp, 1, sd)
post_qnt_smp <- apply(smp, 1, quantile, c(0.025, 0.5, 0.975))
post_mean_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps
post_var_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(delta_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)  

posterior_bps <- t(smp)
colnames(posterior_bps) <- c("sigma^2", "beta[1]", "beta[2]", "beta[3]")

# Plotting data --------------------------------------------------

# plot data on worldmap
newmap <- getMap(resolution = "low")

# world surface interpolation
h <- 12
surf.raw.train <- mba.surf(cbind(crd_s, y), no.X = 300, 
                           no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est
surf.raw.test <- mba.surf(cbind(crd_u, y_u), no.X = 300, 
                          no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

surf.brks <- classIntervals(surf.raw.train$z, 50, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11, 'RdBu')[1:11])
xlim <- range(crd_u[, 1])
zlim <- range(surf.raw.train[["z"]][which(!is.na(surf.raw.train[["z"]]))])

# ggplot version plot
surf_df <- as.data.frame(surf.raw.train)
(combined_train <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = rev(col.pal(length(surf.brks)-1)), limits = zlim) +
    geom_sf(data = st_as_sf(newmap), fill = "#E1F1F1") +
    coord_sf(xlim = c(-175, 175), ylim = c(-55, 85), expand = FALSE) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

surf_df <- as.data.frame(surf.raw.test)
(combined_test <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = rev(col.pal(length(surf.brks)-1)), limits = zlim) +
    geom_sf(data = st_as_sf(newmap), fill = "#E5F5FF") +
    coord_sf(xlim = c(-175, 175), ylim = c(-60, 85), expand = FALSE) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))


# Plotting results ------------------------------------------------------------

h <- 12
surf.raw.pred <- mba.surf(cbind(crd_u, post_mean_Y), no.X = 300, 
                          no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est

# ggplot version plot
surf_df <- as.data.frame(surf.raw.pred)
(combined_pred <- ggplot() +
    geom_tile(data = na.omit(surf_df), aes(x = x, y = y, fill = z), na.rm = TRUE) +
    scale_fill_gradientn(colours = rev(col.pal(length(surf.brks)-1)), limits = zlim) +
    geom_sf(data = st_as_sf(newmap), fill = "#E5F5FF") +
    coord_sf(xlim = c(-175, 175), ylim = c(-60, 85), expand = FALSE) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()))

# graphical UC for Y (ordered)
ord_y <- order(y_u)
set.seed(1995)
plot_ind <- sample(1:u, 250)
df <- data.frame(
  x_ax = (1:u)[plot_ind],
  Yu_ord = y_u[ord_y][plot_ind],
  CI_Y_lower = post_qnt_Y[, ord_y][1,][plot_ind],
  CI_Y_upper = post_qnt_Y[, ord_y][2,][plot_ind],
  Ymap_ord = post_mean_Y[ord_y][plot_ind])
# Create the ggplot
(uc_Y <- ggplot(df, aes(x = x_ax, y = Yu_ord)) +
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

# plotting results
width <- 360*5
height <- 360*2
pointsize <- 12
png("output/dataanalysis_univariate250.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(combined_train, combined_test, uc_Y, combined_pred, nrow = 2, ncol = 2)
dev.off()

# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("time"    = elapsed_times,
                "comb"    = comb_bps,
                "pred"    = predictions,
                "metrics" = c("RMSPE" = rmspe_Y, "naive" = rmspe_naive, "CI_len" = CI_avlen_bps, "ECoverage" = coverage_Y))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "output/dataanalysis_univariate250.RData")

# to read results
# load("output/dataanalysis_univariate250.RData")
# results$time; cat("minutes elapsed for fully model-based uncertainty quantification : \n"); round(results$time/60, 2)
# results$metrics
# results$comb$W
# results$pred
