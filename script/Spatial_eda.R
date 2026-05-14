#####################################################################################################################################################
## Multivariate SPATIAL EXPLORATORY DATA ANALYSIS #########################################
rm(list = ls())
gc()
setwd(".../Bayesian-Transfer-Learning-for-GeoAI")

# re-install damaged packages
remove.packages("gstat")
install.packages("gstat", dependencies = T)
library(gstat)

# Packages --------------------------------------------------------------------
library(spBPS)
library(Rcpp)
library(RcppArmadillo)
library(mniw)
library(ggplot2)
library(patchwork)
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
library(corrplot)
library(mapproj)
library(ggmap)
library(maps)
library(gridExtra)
library(rnaturalearth)
library(rnaturalearthdata)
library(h2o)
library(dplyr)


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
crd_S <- matrix(as.matrix(train_data[, c("lon","lat")]), ncol = 2)
Y_S   <- matrix(as.matrix(train_data[, response_set]), ncol = q)
X_S   <- cbind(1, matrix(as.matrix(train_data[, predictor_set]), ncol = (p-1)))

# define test dimensions
U <- nrow(test_data)
crd_U <- matrix(as.matrix(test_data[, c("lon","lat")]), ncol = 2)
Y_U   <- matrix(as.matrix(test_data[, response_set]), ncol = q)
X_U   <- cbind(1, matrix(as.matrix(test_data[, predictor_set]), ncol = (p-1)))


# EDA -------------------------------------------------------------------------

# total number of locations
print(nrow(full_data))

# computing the maximum inter-site distance (train data)
d.max <- sqrt((max(crd_S[,1]) - min(crd_S[,1]))^2 +
                (max(crd_S[,2]) - min(crd_S[,2]))^2)
d.max # around 386.0908 (multiply by 111.139) ~> 42,909 KM

# statistical summary/feature of response variables
summary(full_data)[,3:4]
apply(full_data[,3:4], 2, sd)

# lnNDVI            lnRedRefl    
# Mean     :8.593   Mean     :8.563
# Std.Dev. :0.517   Std.Dev. :0.447
# Min.     :6.909   Min.     :8.007  
# Max.     :9.469   Max.     :9.472

# plot histogram for tables
png(filename = "hist_lNDVI.png", width = 360, height = 120, bg = "transparent")
par(mar = c(0, 0, 0, 0))
hist(full_data[,3], col = 1, axes = F, main = "", xlab = "", ylab = "")
dev.off()

png(filename = "boxplot_lNDVI.png", width = 360, height = 120, bg = "transparent")
par(mar = c(0, 0, 0, 0))
boxplot(full_data[,3], col = 0, axes = F, main = "", xlab = "", ylab = "", horizontal = T)
dev.off()

png(filename = "hist_lRR.png", width = 360, height = 120, bg = "transparent")
par(mar = c(0, 0, 0, 0))
hist(full_data[,4], col = 1, axes = F, main = "", xlab = "", ylab = "")
dev.off()

png(filename = "boxplot_lRR.png", width = 360, height = 120, bg = "transparent")
par(mar = c(0, 0, 0, 0))
boxplot(full_data[,4], col = 0, axes = F, main = "", xlab = "", ylab = "", horizontal = T)
dev.off()


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
Sigma_samples <- samples$Sigma_samples


# NON-SPATIAL EDA  --------------------------------------------------------

# Summary statistics on non-spatial association
Sigma_mean <- apply(Sigma_samples, 1:2, mean)
apply(Sigma_samples, 1:2, quantile, c(0.025, 0.5, 0.975))
Cor_samples <- sapply(1:dim(Sigma_samples)[3], function(i) {Sigma_samples[1, 2, i]/sqrt(Sigma_samples[1, 1, i]*Sigma_samples[2, 2, i])})
quantile(Cor_samples, c(0.025, 0.5, 0.975))

# Response empirical correlation matrix
cor(Y_S)


# SPATIAL EDA -------------------------------------------------------------

# subsample for feasible variograms computation
set.seed(1997)
eda_ind <- sample.int(N, round(N*0.025))

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

# free memory
gc()

# convert data to sp object (gstat works with sp)
dat <- cbind(crd_S, Y_S)[eda_ind,]
coord_cols <- 1:2
colnames(dat) <- c("lon", "lat", "lnNDVI", "lnRR")
dat <- data.frame(dat)
head(dat)

# prepare data for cross-variogram
coordinates(dat) <- ~lon+lat
gridded(dat) <- TRUE
dat <- as(dat, "SpatialPointsDataFrame") # to full grid
# str(dat)

# create gstat object for cross variogram
g <- gstat(NULL, id = "ndvi", formula = as.formula(paste("lnNDVI", "~1")), data = dat)
g <- gstat(g, id = "rr", formula = as.formula(paste("lnRR", "~1")), data = dat)

# cross-variogram
vgm_cross <- gstat::variogram(g, cross = T, boundaries = seq(0, 0.8*d.max, length.out = 50))
g_cross <- gstat(g, model = vgm(psill = 1, model = "Exp", range = 10, nugget = 0.5), fill.all = TRUE)
g.fit <- fit.lmc(vgm_cross, g_cross, fit.ranges = T)


# plotting eda ------------------------------------------------------------

vario_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text  = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

fitted_vario_df <- function(vfit, n = 300) {
  h <- seq(0, vfit$max.dist, length.out = n)
  data.frame(
    dist  = h,
    gamma = vfit$nugget + vfit$cov.pars[1] * (1 - exp(-h / vfit$cov.pars[2]))
  )
}

plot_vario <- function(v_res, vfit, title,
                       lbl_nugget_x, lbl_nugget_y,
                       lbl_sill_x,   lbl_sill_y,
                       lbl_range_x,  lbl_range_y) {
  df_pts <- data.frame(dist = v_res$u, gamma = v_res$v)
  df_fit <- fitted_vario_df(vfit)
  
  ggplot(df_pts, aes(x = dist, y = gamma)) +
    geom_point(size = 1.5) +
    geom_line(data = df_fit, aes(x = dist, y = gamma),
              color = "steelblue", linewidth = 1, linetype = "dashed") +
    geom_hline(yintercept = vfit$nugget,
               color = "red", linewidth = 1) +
    annotate("text", x = lbl_nugget_x, y = lbl_nugget_y, parse = TRUE,
             label = paste0("tau^2 == ", round(vfit$nugget, 2)),
             color = "red", size = 6, hjust = 0) +
    geom_hline(yintercept = vfit$cov.pars[1],
               color = "green4", linewidth = 0.9, linetype = "dashed") +
    annotate("text", x = lbl_sill_x, y = lbl_sill_y, parse = TRUE,
             label = paste0("sigma^2 == ", round(vfit$cov.pars[1], 2)),
             color = "green4", size = 6, hjust = 0) +
    geom_vline(xintercept = vfit$practicalRange,
               color = "green4", linewidth = 0.9, linetype = "dashed") +
    annotate("text", x = lbl_range_x, y = lbl_range_y, parse = TRUE,
             label = paste0("rho[0] == ", round(vfit$practicalRange, 2)),
             color = "green4", size = 6, hjust = 0) +
    coord_cartesian(ylim = c(0, max(df_pts$gamma) * 1.15),
                    xlim = c(0, vfit$max.dist)) +
    labs(x = "Distance", y = "Semivariance", title = title) +
    vario_theme
}

p_v1 <- plot_vario(v.res_y1, vario.fit_y1,
                   title        = "", #expression(Y[1] ~ "(NDVI)"),
                   lbl_nugget_x = 180, lbl_nugget_y = vario.fit_y1$nugget + 0.02,
                   lbl_sill_x   = 220, lbl_sill_y   = vario.fit_y1$cov.pars[1] - 0.02,
                   lbl_range_x  = vario.fit_y1$practicalRange + 3,
                   lbl_range_y  = max(v.res_y1$v) * 0.45)

p_v2 <- plot_vario(v.res_y2, vario.fit_y2,
                   title        = "", #expression(Y[2] ~ "(RR)"),
                   lbl_nugget_x = 180, lbl_nugget_y = vario.fit_y2$nugget + 0.02,
                   lbl_sill_x   = 220, lbl_sill_y   = vario.fit_y2$cov.pars[1] - 0.02,
                   lbl_range_x  = vario.fit_y2$practicalRange + 3,
                   lbl_range_y  = max(v.res_y2$v) * 0.45)

png("output/variogram.png", width = 12, height = 4, units = "in", res = 300, family = "sans")
print(p_v1 | p_v2)
dev.off()

n_lags   <- 50
df_cv    <- data.frame(dist  = vgm_cross$dist[1:n_lags],
                       gamma = vgm_cross$gamma[1:n_lags])

vl_cross <- gstat::variogramLine(
  g.fit$model$ndvi.rr,
  dist_vector = seq(0, max(df_cv$dist), length.out = 300)
)

cross_sill  <- g.fit$model$ndvi.rr$psill[2]
cross_range <- 3 * g.fit$model$ndvi.rr$range[2]
y_rng       <- range(df_cv$gamma)

p_cross <- ggplot(df_cv, aes(x = dist, y = gamma)) +
  geom_point(size = 2.0) +
  geom_line(data = vl_cross, aes(x = dist, y = gamma),
            color = "steelblue", linewidth = 1.2, linetype = "dashed") +
  geom_hline(yintercept = cross_sill,
             color = "green4", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = max(df_cv$dist) * 0.75,
           y = cross_sill + diff(y_rng) * 0.079, parse = TRUE,
           label = paste0("sigma[list(NDVI,RR)] == ", round(cross_sill, 2)),
           color = "green4", size = 6, hjust = 0) +
  geom_vline(xintercept = cross_range,
             color = "green4", linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = cross_range + max(df_cv$dist) * 0.02,
           y = y_rng[1] + diff(y_rng) * 0.7, parse = TRUE,
           label = paste0("rho[0] == ", round(cross_range, 2)),
           color = "green4", size = 6, hjust = 0) +
  coord_cartesian(ylim = c(y_rng[1] * 1.1, max(df_cv$gamma) * 1.15),
                  xlim = c(0, max(df_cv$dist))) +
  labs(x = "Distance", y = "Semivariance",
       title = "") + #expression("Cross-variogram: " ~ Y[1] ~ "×" ~ Y[2])) +
  vario_theme

png("output/crossvariogram.png", width = 10, height = 5, units = "in", res = 300, family = "sans")
print(p_cross)
dev.off()