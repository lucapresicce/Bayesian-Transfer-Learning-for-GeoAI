#####################################################################################################################################################
## Multivariate DATA ANALYSIS - univariate model #########################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-for-GeoAI")

# Packages --------------------------------------------------------------------
library(spBPS)
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
library(sf)
library(geoR)
library(spBayes)
library(bayesplot)
library(raster)
library(corrplot)
library(mapproj)
library(rworldmap)
library(ggmap)
library(maps)
library(gridExtra)
library(rnaturalearth)
library(rnaturalearthdata)

# Data loading ----------------------------------------------------------------

# load preprocessed RData
load("data/NDVI_data_2024_05_12.RData")
full_data <- NDVIdata
rm(list = c("NDVIdata"))
names(full_data)

# take a look to data structure and responses variables
head(full_data)
dim(full_data)

# train data and test data
set.seed(1997)
test_ind <- sample.int(nrow(full_data), floor(0.15 * nrow(full_data)) )
train_data <- full_data[-test_ind, ]
test_data<- full_data[test_ind, ]

# select responses and predictors sets
response_set <- c("lnNDVI", "lnRedRefl")
q <- length(response_set)
predictor_set <- c("lnZenith")
p <- length(predictor_set)+1

# define train dimensions
N <- nrow(train_data)
# crd_S <- matrix(as.matrix(train_data[, c("scaled_x","scaled_y")]), ncol = 2)
crd_S <- matrix(as.matrix(train_data[, c("lon","lat")]), ncol = 2)
Y_S   <- matrix(as.matrix(train_data[, response_set]), ncol = q)
X_S   <- cbind(1, matrix(as.matrix(train_data[, predictor_set]), ncol = (p-1)))

# define test dimensions
U <- nrow(test_data)
crd_U <- matrix(as.matrix(test_data[, c("lon","lat")]), ncol = 2)
Y_U   <- matrix(as.matrix(test_data[, response_set]), ncol = q)
X_U   <- cbind(1, matrix(as.matrix(test_data[, predictor_set]), ncol = (p-1)))

# remove full dataset and free memory
rm(list = c("full_data", "train_data", "test_data"))
gc()


# EDA -------------------------------------------------------------------------

# subsample for feasible EDA
set.seed(1997)
eda_ind <- sample.int(N, round(N*0.025))

# computing the maximum distance
d.max <- sqrt((max(crd_S[,1]) - min(crd_S[,1]))^2 +
                (max(crd_S[,2]) - min(crd_S[,2]))^2)
d.max # around 386.0908 (multiply by 111.139) ~> 42,909 KM

# check the variogram for the first response
v.res_y1 <- variog(coords = crd_S[eda_ind, ], data = Y_S[eda_ind, 1],
                   uvec = (seq(0, 0.8*d.max, length = 50)))

vario.fit_y1 <- variofit(v.res_y1, cov.model="exponential")
summary(vario.fit_y1)
print(vario.fit_y1)

# check the variogram for the second response
v.res_y2 <- variog(coords = crd_S[eda_ind, ], data = Y_S[eda_ind, 2],
                   uvec = (seq(0, 0.8*d.max, length = 50)))

vario.fit_y2 <- variofit(v.res_y2, cov.model="exponential")
summary(vario.fit_y2)
print(vario.fit_y2)

# plotting the variogram
# plotting results
width <- 360*3
height <- 360
pointsize <- 12
png("output/eda_multivariate.png", width = width, height = height, pointsize = pointsize, family = "sans")

par(mfrow = c(1, 2))
plot(v.res_y1, pch = 19, frame.plot = FALSE, axes = F, ylim = c(0, max(v.res_y1$v)), xlim = c(0, vario.fit_y1$max.dist))
box(bty="l")
axis(2)
axis(1) 
lines(vario.fit_y1, col = 4, lwd = 2, lty = 2)
abline(h = vario.fit_y1$nugget, col = 2, lwd = 3)
text(170, 0.05, 
     bquote(tau^2 ~ " = " ~ .(round(vario.fit_y1$nugget, 2))),
     col = 2, cex = 2)
abline(h = vario.fit_y1$cov.pars[1], col = "green4", lwd = 2, lty = 2)
text(200, 0.25, 
     bquote(sigma^2 ~ " = " ~ .(round(vario.fit_y1$cov.pars[1], 2))), 
     col = "green4", cex = 1.5)
abline(v = vario.fit_y1$practicalRange, col = "green4", lwd = 2, lty = 2)
text(90, 0.175, 
     bquote(rho[0] ~ " = " ~ .(round(vario.fit_y1$practicalRange, 2))), 
     col = "green4", cex = 1.5)

plot(v.res_y2, pch = 19, frame.plot = FALSE, axes = F, ylim = c(0, max(v.res_y2$v)), xlim = c(0, vario.fit_y2$max.dist))
box(bty="l")
axis(2)
axis(1)
lines(vario.fit_y2, col = 4, lwd = 2, lty = 2)
abline(h = vario.fit_y2$nugget, col = 2, lwd = 3)
text(170, 0.06, 
     bquote(tau^2 ~ " = " ~ .(round(vario.fit_y2$nugget, 2))),
     col = 2, cex = 2)
abline(h = vario.fit_y2$cov.pars[1], col = "green4", lwd = 2, lty = 2)
text(200, 0.175, 
     bquote(sigma^2 ~ " = " ~ .(round(vario.fit_y2$cov.pars[1], 2))), 
     col = "green4", cex = 1.5)
abline(v = vario.fit_y2$practicalRange, col = "green4", lwd = 2, lty = 2)
text(120, 0.125, 
     bquote(rho[0] ~ " = " ~ .(round(vario.fit_y2$practicalRange, 2))), 
     col = "green4", cex = 1.5)
par(mfrow = c(1, 1))

dev.off()

# Response correlation matrix
cor(Y_S)

# free memory
gc()

# SubSubsample for model fitting ----------------------------------------------

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

# Prior parameters
mu_B <- matrix(0, p, q)
V_B <- diag(10, p)
nu <- 3
Psi <- diag(q)

# Run the Gibbs sampler
n_iter <- 1000
burn_in <- 500
tic()
set.seed(1234)
samples <- spBPS::bayesMvLMconjugate(y, x, mu_B, V_B, nu, Psi, n_iter, burn_in)
fit_conj_time <- toc()

# Extract posterior samples
B_samples <- samples$B_samples
Sigma_samples <- samples$Sigma_samples

# Summary statistics
B_mean <- apply(B_samples, 1:2, mean)
apply(B_samples, 1:2, quantile, c(0.025, 0.5, 0.975))
Sigma_mean <- apply(Sigma_samples, 1:2, mean)
apply(Sigma_samples, 1:2, quantile, c(0.025, 0.5, 0.975))

tic()
naive <- spBPS::pred_bayesMvLMconjugate(x_u, B_samples, Sigma_samples)
pred_conj_time <- toc()

(rmspe_naive <- sqrt(colMeans((y_u - naive$Y_pred)^2))); mean(rmspe_naive)
(mape_naive <- colMeans( abs(y_u - naive$Y_pred) ) ); mean(mape_naive)

post_cor_Y <- apply(naive$Y_pred_samples, 3, cor)[2,]
post_cor_Y |> quantile(c(0.025, 0.5, 0.975))

tot_conj_time <- as.numeric(fit_conj_time$toc-fit_conj_time$tic) + as.numeric(pred_conj_time$toc-pred_conj_time$tic)


#####################################################################################################################################################

# Subset posterior models -----------------------------------------------------

# hyperparameters values by looking at variograms
variofitphi1 <- 1 / vario.fit_y1$cov.pars[2]; variofitphi1
variofitalpha1 <- vario.fit_y1$cov.pars[1] / (vario.fit_y1$nugget+vario.fit_y1$cov.pars[1]); variofitalpha1
variofitphi2 <- 1 / vario.fit_y2$cov.pars[2]; variofitphi2
variofitalpha2 <- vario.fit_y2$cov.pars[1] / (vario.fit_y2$nugget+vario.fit_y2$cov.pars[1]); variofitalpha2

# define hyperparameters sets
(alfa_seq <- sort(c(variofitalpha1, variofitalpha2)))
(phi_seq <- sort(2*c(variofitphi1, variofitphi2)))

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
                                X_u = x_u, crd_u = crd_u,
                                priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                              V_r = diag(10, p),
                                              Psi = diag(1, q),
                                              nu = 3),
                                hyperpar = list(alpha = alfa_seq, phi = phi_seq),
                                W = Ws, R = 1)
  
  return(result)
}

# subsetting data -----------------------------------------------------------

subset_size <- 250
K <- n/subset_size
data_part <- subset_data(data = list(Y = y, X = x, crd = crd_s), K = K)

# spBPS parallel fit -------------------------------------------------------

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
pred_mat_W <- sapply(1:R, function(r){matrix(predictions[[r]]$Pred[[1]]$Wu, nrow = u, ncol = q)}, simplify = "array")
post_mean_W <- apply(pred_mat_W, c(1,2), mean)
post_var_W <- apply(pred_mat_W, c(1,2), sd)
post_qnt_W <- apply(pred_mat_W, c(1,2), quantile, c(0.025, 0.975))

# statistics computations Y
pred_mat_Y <- sapply(1:R, function(r){matrix(predictions[[r]]$Pred[[1]]$Yu, nrow = u, ncol = q)}, simplify = "array")
post_mean_Y <- apply(pred_mat_Y, c(1,2), mean)
post_var_Y <- apply(pred_mat_Y, c(1,2), sd)
post_qnt_Y <- apply(pred_mat_Y, c(1,2), quantile, c(0.025, 0.975))

# Empirical correlation
post_cor_Y <- apply(pred_mat_Y, 3, cor)[2,]
post_cor_Y |> quantile(c(0.025, 0.5, 0.975))

# Empirical coverage for Y
coverage_Y <- c(mean(y_u[,1] >= post_qnt_Y[1,,1] & y_u[,1] <= post_qnt_Y[2,,1]),
                mean(y_u[,2] >= post_qnt_Y[1,,2] & y_u[,2] <= post_qnt_Y[2,,2]))
# mean(Y_u >= post_qnt_Y[1,,] & Y_u <= post_qnt_Y[2,,])
cat("Empirical average coverage for Spatial process:", round(mean(coverage_Y), 3))
(CI_avlen_bps <- mean(post_qnt_Y[2,,]-post_qnt_Y[1,,]))

# Root Mean Square Prediction Error
(rmspe_Y <- sqrt( colMeans( (y_u - post_mean_Y)^2 ) )); mean(rmspe_Y)
(mape_Y <- colMeans( abs(y_u - post_mean_Y) ) ); mean(mape_Y)

# Posterior inference -----------------------------------------------------

beta_smp <- sapply(1:R, function(r){predictions[[r]]$Post[[1]]$beta[1:p,]}, simplify = "array")
post_mean_beta <-  apply(beta_smp, c(1,2), mean)
post_var_beta <- apply(beta_smp, c(1,2), sd)
post_qnt_beta <- apply(beta_smp, c(1,2), quantile, c(0.025, 0.5, 0.975))

sigma_smp <- sapply(1:R, function(r){predictions[[r]]$Post[[1]]$sigma}, simplify = "array")
post_mean_sigma <- apply(sigma_smp, c(1,2), mean)
post_var_sigma <- apply(sigma_smp, c(1,2), sd)
post_qnt_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.025, 0.5, 0.975))

post_mean_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]]) %*% Wbps
post_var_hyp <- sapply(1:K, function(k)t(expand_grid_cpp(alfa_seq, phi_seq)) %*% W_list[[k]])^2 %*% Wbps - (post_mean_hyp^2)

posterior_bps <- cbind(t(sapply(1:R, function(r)matrix(beta_smp[,,r]))),
                       t(sapply(1:R, function(r)matrix(sigma_smp[,,r])))[,-3])
colnames(posterior_bps) <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]", "Sigma[1,1]", "Sigma[2,1]", "Sigma[2,2]")

# Save timing result ----------------------------------------------------------

elapsed_times <- c("Fitting" = as.numeric(fit_time$toc-fit_time$tic),
                   "Combination" = as.numeric(comb_time$toc-comb_time$tic),
                   "Prediction" = as.numeric(prd_time$toc-prd_time$tic),
                   "Total time" = as.numeric(tot_time$toc-tot_time$tic))

cat("minutes elapsed for fully model-based uncertainty quantification : \n"); round(elapsed_times/60, 2)


# Plotting result ------------------------------------------------

# Crea un data frame con i dati
df_train <- data.frame(lon = crd_s[,1], lat = crd_s[,2],
                       NDVI = y[,1], RedRefl = y[,2])
df_test <- data.frame(lon = crd_u[,1], lat = crd_u[,2],
                      NDVI = y_u[,1], RedRefl = y_u[,2], hatNDVI = post_mean_Y[,1], hatRedRefl = post_mean_Y[,2])

# Funzione per creare la superficie interpolata
create_interpolated_surface <- function(data, variable) {
  h <- 12
  mba_result <- MBA::mba.surf(data[, c("lon", "lat", variable)], no.X = 300, no.Y = 300, exten = F, sp = TRUE, h = h)$xyz.est
  return(mba_result)
}

# Crea le superfici interpolate
NDVI_surface      <- create_interpolated_surface(df_train, "NDVI")
RR_surface        <- create_interpolated_surface(df_train, "RedRefl")
NDVI_test_surface <- create_interpolated_surface(df_test, "NDVI")
RR_test_surface   <- create_interpolated_surface(df_test, "RedRefl")
NDVI_hat_surface  <- create_interpolated_surface(df_test, "hatNDVI")
RR_hat_surface    <- create_interpolated_surface(df_test, "hatRedRefl")

# Funzione per convertire la superficie in data frame per ggplot
surface_to_df <- function(surface) {
  df <- data.frame(lon = surface$x, lat = surface$y, z = as.vector(surface$z))
  return(df)
}

# Converte le superfici in data frame
NDVI_df      <- surface_to_df(NDVI_surface)
RR_df        <- surface_to_df(RR_surface) 
NDVI_test_df <- surface_to_df(NDVI_test_surface) 
RR_test_df   <- surface_to_df(RR_test_surface)
NDVI_hat_df  <- surface_to_df(NDVI_hat_surface)
RR_hat_df    <- surface_to_df(RR_hat_surface)

# Ottenere i contorni dei continenti
world <- ne_countries(scale = "medium", returnclass = "sf")

# Ottenere i contorni dei continenti e convertirli in raster per la maschera
world <- ne_countries(scale = "medium", returnclass = "sf")
world_raster <- raster(extent(-180, 180, -90, 90), resolution = 1)
world_raster <- rasterize(world, world_raster, field = 1, background = NA)

# Funzione per applicare la maschera dei continenti ai dati interpolati
apply_continent_mask <- function(df, raster_mask) {
  coords <- data.frame(lon = df$lon, lat = df$lat)
  mask_values <- raster::extract(raster_mask, coords)
  df[!is.na(mask_values), ]
}

# Applica la maschera dei continenti ai dati interpolati
NDVI_df      <- apply_continent_mask(NDVI_df, world_raster)
RR_df        <- apply_continent_mask(RR_df, world_raster)
NDVI_test_df <- apply_continent_mask(NDVI_test_df, world_raster)
RR_test_df   <- apply_continent_mask(RR_test_df, world_raster)
NDVI_hat_df  <- apply_continent_mask(NDVI_hat_df, world_raster)
RR_hat_df    <- apply_continent_mask(RR_hat_df, world_raster)

# Crea le mappe con ggplot2
plot_map_col <- function(data, title, fill_label, low_color, high_color) {
  ggplot(data, aes(x = lon, y = lat, fill = z)) +
    geom_raster(interpolate = TRUE) +
    borders("world", colour = "black", fill = NA) +
    scale_fill_gradient(low = low_color, high = high_color, na.value = NA, name = fill_label) +
    coord_fixed(ratio = 1.3, ylim = c(-55, 80), xlim = c(-180, 180)) +
    # theme_void() +
    # labs(title = title)
    labs(x = "", y = "") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line()) +
    scale_x_continuous(breaks=c(-120, -60, 0, 60, 120),
                       labels=c("120°W", "60°W", "0°", "60°E", "120°E" )) +
    scale_y_continuous(breaks=c(-40, -20, 0, 20, 40, 60, 80),
                       labels=c("40°S", "20°S", "0°", "20°N", "40°N", "60°N", "80°N" ))
}
plot_map <- function(data, title, fill_label, Pal) {
  # Calcola i breakpoints e la palette di colori
  breaks <- classIntervals(data$z, n = 50, style = "pretty")$brks
  col.pal <- rev(colorRampPalette(brewer.pal(9, Pal)[1:9])(length(breaks) - 1))
  # Miglioramo la legenda
  legend_labels <- legend_labels <- pretty(data$z, n = 5)
  breaks <- classIntervals(data$z, n = 5, style = "pretty")$brks
  # Plot
  ggplot(data, aes(x = lon, y = lat, fill = z)) +
    geom_raster(interpolate = TRUE) +
    borders("world", colour = "black", fill = NA) +
    scale_fill_gradientn(colors = col.pal, breaks = breaks, labels = legend_labels, na.value = NA, name = fill_label) +
    coord_fixed(ratio = 1.3, ylim = c(-55, 80), xlim = c(-180, 180)) +
    # theme_void() +
    # labs(title = title)
    labs(x = "", y = "") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(), plot.margin = unit(c(-.05, -.05, -.05, -.05), "cm")) +
    scale_x_continuous(breaks=c(-120, -60, 0, 60, 120),
                       labels=c("120°W", "60°W", "0°", "60°E", "120°E" )) +
    scale_y_continuous(breaks=c(-40, -20, 0, 20, 40, 60, 80),
                       labels=c("40°S", "20°S", "0°", "20°N", "40°N", "60°N", "80°N" ))
}

# Crea i quattro grafici
NDVI_plot      <- plot_map(NDVI_df, "NDVI - train set", "NDVI"            , Pal = "YlGn")
RR_plot        <- plot_map(RR_df, "Red Reflectance - train set", "RedRefl", Pal = "RdBu")
NDVI_test_plot <- plot_map(NDVI_test_df, "NDVI - train set", "NDVI"            , Pal = "YlGn")
RR_test_plot   <- plot_map(RR_test_df, "Red Reflectance - train set", "RedRefl", Pal = "RdBu")
NDVI_hat_plot  <- plot_map(NDVI_hat_df, "NDVI - train set", "NDVI"            , Pal = "YlGn")
RR_hat_plot    <- plot_map(RR_hat_df, "Red Reflectance - train set", "RedRefl", Pal = "RdBu")

# graphical UC for Y1 (ordered)
ord_y <- order(y_u[,1])
set.seed(1997)
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
    labs(x = "Ordered locations", y = "NDVI") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"), plot.margin = unit(c(.75, .75, .75, .75), "cm")))


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
    labs(x = "Ordered locations", y = "Red Reflectance") +
    geom_point(aes(y = Ymap_ord), pch = 15, size = 1.5, col = "#D41159") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"), plot.margin = unit(c(.75, .75, .75, .75), "cm")))

# Combina i grafici in un'unica figura (per ogni risposta)
# NDVI
width <- 360*5
height <- 360*2
pointsize <- 12
png("output/dataanalysis_multivariate_NDVI250.png", width = width, height = height, pointsize = pointsize, family = "sans")
grid.arrange(NDVI_plot, NDVI_test_plot, uc_Y1, NDVI_hat_plot,  nrow = 2, ncol = 2)
dev.off()
# RR
width <- 360*5
height <- 360*2
pointsize <- 12
png("output/dataanalysis_multivariate_RR250.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(RR_plot, RR_test_plot, uc_Y2, RR_hat_plot, nrow = 2, ncol = 2)
dev.off()



# Save results ----------------------------------------------------------------

gc()
# Save the entire environment
results <- list("time"    = elapsed_times,
                "comb"    = comb_bps,
                "pred"    = predictions,
                "metrics" = c("RMSPE_bps" = rmspe_Y, "RMPSE_naive" = rmspe_naive, "CI_len" = CI_avlen_bps, "ECoverage" = coverage_Y, "MAPE_bps" = mape_Y, "MAPE_naive" = mape_naive))

rm(list = ls()[which(!(ls() %in% c("results")))])
save.image(file = "output/dataanalysis_multivariate250.RData")

# to read results
# load("output/dataanalysis_multivariate.RData")
# results$time; cat("minutes elapsed for fully model-based uncertainty quantification : \n"); round(results$time/60, 2)
# results$metrics
# results$comb
# results$pred
