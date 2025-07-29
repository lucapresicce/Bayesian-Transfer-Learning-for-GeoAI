# Load necessary libraries
library(keras)
library(spBPS)
library(akima)
library(ggplot2)

# Define dimensions
n <- 500      # Number of rows in input matrix
q <- 2
p <- 2
qq <- q+p     # Number of columns in input matrix
pp <- (p*q)+((q*q)-1)
nq <- n*q
m  <- pp + nq # Number of rows in output matrix
t  <- 3       # Number of columns in output matrix

# Fixed coordinates and spatial process (used for all data)
set.seed(42)
fixed_crd <- matrix(runif(n * 2), ncol = 2)
sigma2 <- matrix(c(1, -0.3, -0.3, 1), 2, 2)
phi <- 4
Rphi_val <- exp(-phi * arma_dist(fixed_crd))
fixed_W <- mniw::rMNorm(1, Lambda = matrix(0, n, q), SigmaR = Rphi_val, SigmaC = sigma2)


# Data generation ---------------------------------------------------------

data_gen <- function(n, crd, W) {
  
  p <- 2
  q <- 2
  B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
  sigma2 <- matrix(c(1, -0.3, -0.3, 1), q, q)
  alfa <- 0.8
  
  X <- cbind(rep(1, n), matrix(runif((p - 1) * n), ncol = p - 1))
  Y <- X %*% B + W + mniw::rMNorm(1, Lambda = matrix(0, n, q),
                                  SigmaR = diag((1 / alfa) - 1, n),
                                  SigmaC = sigma2)
  return(cbind(Y, X))
}


# Posterior generation ----------------------------------------------------

calculate_post <- function(Z, crd) {
  
  p <- 2
  q <- 2
  pp <- p + q + 2
  
  alfa_seq <- c(0.7, 0.8, 0.9)
  phi_seq <- c(3, 4, 5)
  
  K <- 5
  Y <- Z[, 1:q]
  X <- Z[, -(1:q)]
  data_part <- subset_data(data = list(Y = Y, X = X, crd = crd), K = K)
  
  fit_loop <- function(i) {
    Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
    bps <- spBPS::BPS_weights_MvT(data = list(Y = Yi, X = Xi),
                                  priors = list(mu_B = matrix(0, p, q),
                                                V_r = diag(10, p),
                                                Psi = diag(1, q),
                                                nu = 3),
                                  coords = crd_i,
                                  hyperpar = list(alpha = alfa_seq, phi = phi_seq), K = 5)
    list(bps$epd, bps$W)
  }
  
  pred_loop <- function(r) {
    ind_s <- subset_ind[r]
    Ys <- data_part$Y_list[[ind_s]]; Xs <- data_part$X_list[[ind_s]]
    crds <- data_part$crd_list[[ind_s]]; Ws <- W_list[[ind_s]]
    spBPS::BPS_post_MvT(data = list(Y = Ys, X = Xs), coords = crds,
                        X_u = X, crd_u = crd,
                        priors = list(mu_B = matrix(0, p, q),
                                      V_r = diag(10, p),
                                      Psi = diag(1, q),
                                      nu = 3),
                        hyperpar = list(alpha = alfa_seq, phi = phi_seq),
                        W = Ws, R = 1)
  }
  
  cl <- parallel::makeCluster(parallel::detectCores(logical = FALSE) - 1)
  doParallel::registerDoParallel(cl)
  
  obj_fit <- foreach::foreach(i = 1:K) %dopar% { fit_loop(i) }
  comb_bps <- spBPS::BPS_combine(obj_fit, K, 1)
  Wbps <- comb_bps$W
  W_list <- comb_bps$W_list
  R <- 250
  subset_ind <- sample(1:K, R, TRUE, Wbps)
  predictions <- foreach::foreach(r = 1:R) %dopar% { pred_loop(r) }
  
  stopCluster(cl)
  gc()
  
  beta_smp <- sapply(1:R, function(r){ predictions[[r]]$Post[[1]]$beta[1:p,] }, simplify = "array")
  post_qnt_beta <- apply(beta_smp, c(1,2), quantile, c(0.5, 0.025, 0.975))
  vbet <- matrix(post_qnt_beta, nrow = 3)
  
  sigma_smp <- sapply(1:R, function(r){ predictions[[r]]$Post[[1]]$sigma }, simplify = "array")
  post_qnt_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.5, 0.025, 0.975))
  vsigm <- matrix(cbind(post_qnt_sigma[,,1], post_qnt_sigma[,,2])[, -p], nrow = 3)
  
  Omega_smp <- sapply(1:R, function(r){ predictions[[r]]$Pred[[1]]$Wu }, simplify = "array")
  post_qnt_Omega <- apply(Omega_smp, c(1,2), quantile, c(0.5, 0.025, 0.975))
  vomega <- matrix(cbind(post_qnt_Omega[,,1], post_qnt_Omega[,,2]), nrow = 3)
  
  return(t(cbind(vsigm, vbet, vomega)))
}


# Synthetic dataset -------------------------------------------------------

# Generate synthetic data
set.seed(42)  
num_samples <- 100

# Input data: (num_samples, n, q)
Z <- sapply(1:num_samples, function(i) {data_gen(n, crd = fixed_crd, W = fixed_W)}, simplify = "array")
Z <- aperm(Z, c(3,1,2))

# Output data: (num_samples, m, p)
Y <- sapply(1:num_samples, function(i) {calculate_post(Z = Z[i,,], crd = fixed_crd)}, simplify = "array")
Y <- aperm(Y, c(3,1,2))

# Function to create a residual block
residual_block <- function(x, units) {
  # Apply two dense layers
  x_main <- layer_dense(x, units = units, activation = 'relu')
  x_main <- layer_dense(x_main, units = units)
  
  # Create a shortcut connection that matches the output shape
  shortcut <- layer_dense(x, units = units)
  
  # Add the shortcut connection and apply activation
  x <- layer_add(list(x_main, shortcut))
  x <- layer_activation(x, activation = 'relu')
  
  return(x)
}

# Build the model
input_tensor <- layer_input(shape = c(n, qq))

# Create the ResNet architecture
x <- layer_flatten(input_tensor)
x <- residual_block(x, 128)
x <- residual_block(x, 256)
x <- residual_block(x, 512)
x <- layer_dense(x, units = (m * t), activation = 'linear')
output_tensor <- layer_reshape(x, target_shape = c(m, t))

# Create the model
model <- keras_model(inputs = input_tensor, outputs = output_tensor)

# Compile the model
model$compile(
  loss = 'mean_squared_error',
  optimizer = 'adam',
  metrics = list('mae')
)

# Train the model
model$fit(Z, Y, epochs = 50L, batch_size = 24L, validation_split = 0.2)

# Evaluate the model
score <- model$evaluate(Z, Y)
print(score)

# Generate synthetic input data for predictions
num_pred_samples <- 1  # Number of samples for prediction
X_new <- sapply(1:num_pred_samples, function(i) {data_gen(n, crd = fixed_crd, W = fixed_W)}, simplify = "array")
X_new <- aperm(X_new, c(3,1,2))

# Perform predictions - BPS
predictions_bps <- calculate_post(X_new[1,,], crd = fixed_crd)
print(dim(predictions_bps))

omega_pred_bps <- predictions_bps[-(1:7),]
omega_pred_med_bps <- matrix(omega_pred_bps[,1], ncol = q)
omega_pred_low_bps <- matrix(omega_pred_bps[,2], ncol = q)
omega_pred_upp_bps <- matrix(omega_pred_bps[,3], ncol = q)

rmpse <- (sqrt(colMeans( ( omega_pred_med_bps - fixed_W )^2 )))
cover <- (colMeans( (omega_pred_low_bps < fixed_W) & (omega_pred_upp_bps > fixed_W) ))

# Perform predictions
predictions <- model$predict(X_new)

# Manage predictions
predictions <- predictions[1,,]
pred_par <- predictions[  1:7,]
pred_spa <- predictions[-(1:7),]

# Spatial process prediction
omega_pred_med <- matrix(pred_spa[,1], ncol = q)
omega_pred_low <- matrix(pred_spa[,2], ncol = q)
omega_pred_upp <- matrix(pred_spa[,3], ncol = q)

(rmspe_spat <- sqrt(colMeans( ( omega_pred_med - fixed_W )^2 )))
(cov_spat <- colMeans( (omega_pred_low < fixed_W) & (omega_pred_upp > fixed_W) ))

# Plot surface function
plot_surface_interp <- function(mat, title, component = 1, grid_res = 100) {
  library(RColorBrewer)
  vals <- mat[, component]
  interp_result <- with(data.frame(x = fixed_crd[, 1], y = fixed_crd[, 2], z = vals),
                        interp(x, y, z, nx = grid_res, ny = grid_res, duplicate = "mean"))
  df_grid <- data.frame(
    x = rep(interp_result$x, times = length(interp_result$y)),
    y = rep(interp_result$y, each = length(interp_result$x)),
    value = as.vector(interp_result$z)
  )
  
  ggplot(df_grid, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(colors = brewer.pal(n = 11, name = "RdBu")) + 
    coord_equal() +
    ggtitle(title) +
    theme_minimal()
}

# heatmap plot
p0.1 <- plot_surface_interp(fixed_W,        "Spatial Process - 1st Response"                   , component = 1)
pb.1 <- plot_surface_interp(omega_pred_med_bps, "Double BPS Prediction (Median) - 1st Response", component = 1)
p1.1 <- plot_surface_interp(omega_pred_med, "Amortized Prediction (Median) - 1st Response"     , component = 1)
p2.1 <- plot_surface_interp(omega_pred_low, "Amortized Prediction (2.5th) - 1st Response"      , component = 1)
p3.1 <- plot_surface_interp(omega_pred_upp, "Amortized Prediction (97.5th) - 1st Response"     , component = 1)
p0.2 <- plot_surface_interp(fixed_W,        "Spatial Process - 2nd Response"                   , component = 2)
pb.2 <- plot_surface_interp(omega_pred_med_bps, "Double BPS Prediction (Median) - 2nd Response", component = 2)
p1.2 <- plot_surface_interp(omega_pred_med, "Amortized Prediction (Median) - 2nd Response"     , component = 2)
p2.2 <- plot_surface_interp(omega_pred_low, "Amortized Prediction (2.5th) - 2nd Response"      , component = 2)
p3.2 <- plot_surface_interp(omega_pred_upp, "Amortized Prediction (97.5th) - 2nd Response"     , component = 2)


# Size for the mapping
width <- 360*6
height <- 360*3
pointsize <- 16
png("heatmap-amortized.png", width = width, height = height, pointsize = pointsize, family = "sans")
gridExtra::grid.arrange(p0.1, pb.1, p1.1, p2.1, p3.1, p0.2, pb.2, p1.2, p2.2, p3.2, ncol = 5)
dev.off()

# parameter prediction
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- matrix(c(1, -0.3, -0.3, 1), q, q)
par_true <- c(c(sigma2)[c(1,3,4)], B)


# parameter names
greek_labels <- c("bold(Sigma)[list(1,1)]", "bold(Sigma)[list(2,2)]", "bold(Sigma)[list(1,2)]", "bold(beta)[list(0,1)]", "bold(beta)[list(1,1)]", "bold(beta)[list(0,2)]", "bold(beta)[list(1,2)]")
df <- data.frame(
  parameter = factor(1:7, 1:7),
  true = par_true,
  median = pred_par[, 1],
  lower = pred_par[,2],
  upper = pred_par[,3]
)

# CI plot
ci_plot <- ggplot(df, aes(x = parameter)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, color = "darkblue", linewidth = 2) +
    geom_point(aes(y = median, shape = "Estimated", fill = "Estimated"), size = 4, alpha = 0.85, color = "black") +
    geom_point(aes(y = true, shape = "True", fill = "True"), size = 4, alpha = 0.85, color = "black") +
    scale_x_discrete(labels = function(x) parse(text = greek_labels)) +
    scale_shape_manual(name = NULL, values = c("Estimated" = 21, "True" = 24)) +
    scale_fill_manual(name = NULL, values = c("Estimated" = "darkblue", "True" = "yellow")) +
    guides(shape = guide_legend(override.aes = list(fill = c("darkblue", "yellow"), color = "black", alpha = 0.85))) +
    labs(
      title = "Credible intervals - Amortized Inference",
      subtitle = "with posterior medians, true values, and 95% predicted credible intervals",
      y = "", x = ""
    ) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 16), #, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    )


# Size for the mapping
width <- 360*2
height <- 360
pointsize <- 16
png("parameters-amortized.png", width = width, height = height, pointsize = pointsize, family = "sans")
ci_plot
dev.off()
