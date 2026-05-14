# Load necessary libraries
library(keras)
library(spBPS)
library(akima)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

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
Rphi_val <- exp(-phi * spBPS:::arma_dist(fixed_crd))
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
  hyperpar <- list(alpha = alfa_seq, phi = phi_seq)
  
  K <- 5
  Y <- Z[, 1:q]
  X <- Z[, -(1:q)]
  
  invisible(
    capture.output(
      out <- spBPS(data = list(Y = Y, X = X),
               priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                             V_r = diag(10, p),
                             Psi = diag(1, q),
                             nu = 3),
               coords = crd,
               hyperpar = hyperpar,
               K = K,
               combine_method = "bps",
               draws = 200,
               newdata = list(X = X, coords = crd),
               n_cores = n.core,
               pred_batch_size = 50)
    ))
  
  beta_smp <- out$posterior$beta
  post_qnt_beta <- apply(beta_smp, c(1,2), quantile, c(0.5, 0.025, 0.975))
  vbet <- matrix(post_qnt_beta, nrow = 3)
  
  sigma_smp <- out$posterior$sigma
  post_qnt_sigma <- apply(sigma_smp, c(1,2), quantile, c(0.5, 0.025, 0.975))
  vsigm <- matrix(cbind(post_qnt_sigma[,,1], post_qnt_sigma[,,2])[, -p], nrow = 3)
  
  Omega_smp <- out$predictive$Wu
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


# plotting ----------------------------------------------------------------

# parameter prediction
B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
sigma2 <- matrix(c(1, -0.3, -0.3, 1), q, q)
par_true <- c(c(sigma2)[c(1,3,4)], B)

# parameter names
greek_labels <- c(
  "bold(Sigma)[list(1,1)]", "bold(Sigma)[list(2,2)]", "bold(Sigma)[list(1,2)]",
  "bold(beta)[list(0,1)]",  "bold(beta)[list(1,1)]",  "bold(beta)[list(0,2)]",  "bold(beta)[list(1,2)]"
)

df <- data.frame(
  parameter = factor(1:7, 1:7),
  true      = par_true,
  median    = pred_par[, 1],
  lower     = pred_par[, 2],
  upper     = pred_par[, 3]
)

ci_plot <- ggplot(df, aes(x = parameter)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, color = "darkblue", linewidth = 2) +
  geom_point(aes(y = median, shape = "Estimated", fill = "Estimated"), size = 5, alpha = 0.85, color = "black") +
  geom_point(aes(y = true,   shape = "True",      fill = "True"),      size = 5, alpha = 0.85, color = "black") +
  scale_x_discrete(labels = function(x) parse(text = greek_labels)) +
  scale_shape_manual(name = NULL, values = c("Estimated" = 21, "True" = 24)) +
  scale_fill_manual( name = NULL, values = c("Estimated" = "darkblue", "True" = "yellow")) +
  guides(shape = guide_legend(override.aes = list(
    fill = c("darkblue", "yellow"), color = "black", alpha = 0.85, size = 5
  ))) +
  labs(y = "", x = "") +   # title e subtitle rimossi
  theme_classic(base_size = 16) +
  theme(
    axis.text.y     = element_text(size = 18, face = "bold"),
    axis.text.x     = element_text(size = 20, face = "bold", angle = 45),
    axis.title      = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 15),
    legend.title    = element_blank(),
    legend.position = "right"
  )

png("output/parameters-amortized.png", width = 10, height = 5, units = "in", res = 300, family = "sans")
print(ci_plot)
dev.off()

titles <- list(
  y1 = list(
    expression(bold(Omega[1]) ~ "– True"),
    expression(bold(hat(Omega)[1]) ~ "– DBPS" ~ (p[50]  ) ),
    expression(bold(hat(Omega)[1]) ~ "– ABI"  ~ (p[50]  ) ),
    expression(bold(hat(Omega)[1]) ~ "– ABI"  ~ (p[2.5] ) ),
    expression(bold(hat(Omega)[1]) ~ "– ABI"  ~ (p[97.5]) )
  ),
  y2 = list(
    expression(bold(Omega[2]) ~ "– True"),
    expression(bold(hat(Omega)[2]) ~ "– DBPS" ~ (p[50]  ) ),
    expression(bold(hat(Omega)[2]) ~ "– ABI"  ~ (p[50]  ) ),
    expression(bold(hat(Omega)[2]) ~ "– ABI"  ~ (p[2.5] ) ),
    expression(bold(hat(Omega)[2]) ~ "– ABI"  ~ (p[97.5]) )
  )
)

get_vals <- function(mat, component, grid_res = 100) {
  vals <- mat[, component]
  interp_result <- with(data.frame(x = fixed_crd[, 1], y = fixed_crd[, 2], z = vals),
                        interp(x, y, z, nx = grid_res, ny = grid_res, duplicate = "mean"))
  as.vector(interp_result$z)
}

mat_list <- list(fixed_W, omega_pred_med_bps, omega_pred_med, omega_pred_low, omega_pred_upp)

lims_y1 <- range(unlist(lapply(mat_list, get_vals, component = 1)), na.rm = TRUE)
lims_y2 <- range(unlist(lapply(mat_list, get_vals, component = 2)), na.rm = TRUE)

plot_surface_interp <- function(mat, title, component = 1, grid_res = 100, scale_limits = NULL) {
  vals <- mat[, component]
  interp_result <- with(data.frame(x = fixed_crd[, 1], y = fixed_crd[, 2], z = vals),
                        interp(x, y, z, nx = grid_res, ny = grid_res, duplicate = "mean"))
  df_grid <- data.frame(
    x     = rep(interp_result$x, times = length(interp_result$y)),
    y     = rep(interp_result$y, each  = length(interp_result$x)),
    value = as.vector(interp_result$z)
  )
  
  ggplot(df_grid, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(
      colors = brewer.pal(n = 11, name = "RdBu"),
      limits = scale_limits,
      name   = "Value"
    ) +
    coord_equal() +
    ggtitle(title) +
    labs(x = "Easting", y = "Northing") +
    theme_minimal(base_size = 15) +
    theme(
      plot.title        = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.text         = element_text(size = 12, face = "bold"),
      axis.title        = element_text(size = 14, face = "bold"),
      legend.title      = element_text(size = 14, face = "bold"),
      legend.text       = element_text(size = 12, face = "bold"),
      legend.key.height = unit(1.5, "cm")
    )
}

panels_y1 <- mapply(function(mat, ttl) {
  plot_surface_interp(mat, ttl, component = 1, scale_limits = lims_y1)
}, mat_list, titles$y1, SIMPLIFY = FALSE)

panels_y2 <- mapply(function(mat, ttl) {
  plot_surface_interp(mat, ttl, component = 2, scale_limits = lims_y2)
}, mat_list, titles$y2, SIMPLIFY = FALSE)

row1 <- wrap_plots(panels_y1, nrow = 1) + plot_layout(guides = "collect")
row2 <- wrap_plots(panels_y2, nrow = 1) + plot_layout(guides = "collect")

final_heatmap <- row1 / row2

png("output/heatmap-amortized.png", width = 18, height = 9, units = "in", res = 300, family = "sans")
print(final_heatmap)
dev.off()