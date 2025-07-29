#### asymptotic upper bound - Monte Carlo simulations ####

# packages ----------------------------------------------------------------
library(spBPS)
library(Rcpp)
library(RcppArmadillo)
library(tictoc)


# Working functions -------------------------------------------------------

# loading fitting function from spBPS
Sys.setenv(PKG_CXXFLAGS = "-Ofast"); sourceCpp("code/src/code.cpp")

# Helper function to find most balanced factor pair
get_balanced_factors <- function(J) {
  if (J %% 2 > 1e-10) stop("J must be even.")
  
  factors <- which(J %% 1:J < 1e-10)
  pairs <- lapply(factors, function(f) c(f, J / f))
  best <- pairs[[which.min(sapply(pairs, function(p) abs(p[1] - p[2])))]]
  return(best)
}

# data generation
data_gen <- function(n) {
  
  # dimensions
  q <- 2
  p <- 2
  
  # parameters
  B <- matrix(c(-0.75, 0.90, 1.85, -1.1), p, q)
  sigma2 <- matrix(c(1, -0.3, -0.3, 1), q, q)
  
  # hyperparameters
  alfa <- 0.8
  phi <- 4
  
  # fixed coordinates
  crd <- matrix(runif((n) * 2), ncol = 2)

  # distance matrix
  D <- arma_dist(crd)
  Rphi <- exp(-phi * D)
  
  # generate sintethic data
  X <- cbind(rep(1, n), matrix(runif((p-1)*(n)), ncol = (p-1)))
  W <- matrix(0, n, q) + mniw::rMNorm(1, Lambda = matrix(0, n, q), SigmaR = Rphi, SigmaC = sigma2)
  Y <- X %*% B + W + mniw::rMNorm(1, Lambda = matrix(0, n, q), SigmaR = diag((1/alfa)-1, n), SigmaC = sigma2)
  
  # return( cbind(crd, Y, X) )
  
  return( list("crd" = crd, "Y" = Y, "X" = X) )
  
}

# posterior computation
calculate_post <- function(Z, K, sqJ) {
  
  #  dimensions
  p <- 2
  q <- 2
  
  # hyperparameters values
  # J <- 22
  # sqJ <- sqrt(J)
  J <- ceiling(sqJ^2 - 1e-6)
  if (sqJ == floor(sqJ)) {
    # J is a perfect square
    alfa_seq <- seq(from  = 0.7, to = 0.9, length.out = sqJ)
    phi_seq  <- seq(from  = 3,   to = 5,   length.out = sqJ)
  } else {

    pair <- get_balanced_factors(J)

    if (runif(1) > 0.5) {
      n1 <- pair[1]
      n2 <- pair[2]
    } else {
      n1 <- pair[2]
      n2 <- pair[1]
    }

    alfa_seq <- seq(from  = 0.7, to = 0.9, length.out = as.integer(n1))
    phi_seq  <- seq(from  = 3,   to = 5,   length.out = as.integer(n2))
  }
  # length(alfa_seq)
  # length(phi_seq)
  # nrow(expand_grid_cpp(alfa_seq, phi_seq))
  
  
  
  # subsetting data
  K <- K
  crd <- Z$crd
  Y   <- Z$Y
  X <- Z$X
  data_part <- subset_data(data = list(Y = Y, X = X, crd = crd), K = K)
  
  # function for the loops
  fit_loop <- function(i,K) {
    
    Yi <- data_part$Y_list[[i]]; Xi <- data_part$X_list[[i]]; crd_i <- data_part$crd_list[[i]]
    p <- ncol(Xi); q <- ncol(Yi)
    bps <- spBPS::BPS_weights_MvT(data = list(Y = Yi, X = Xi),
                                  priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                                V_r = diag(10, p),
                                                Psi = diag(1, q),
                                                nu = 3), coords = crd_i,
                                  hyperpar = list(alpha = alfa_seq, phi = phi_seq), K = K)
    w_hat <- bps$W
    epd <- bps$epd
    
    result <- list(epd, w_hat)
    return(result)
    
  }
  funs_fit <- lsf.str()[which(lsf.str() != "fit_loop")]
  pred_loop <- function(r) {
    
    ind_s <- subset_ind[r]
    Ys <- data_part$Y_list[[ind_s]]; Xs <- data_part$X_list[[ind_s]]; crds <- data_part$crd_list[[ind_s]]; Ws <- W_list[[ind_s]]
    result <- spBPS::BPS_post_MvT(data = list(Y = Ys, X = Xs), coords = crds,
                                  X_u = Xs[1:2,], crd_u = crds[1:2,],
                                  priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                                V_r = diag(10, p),
                                                Psi = diag(1, q),
                                                nu = 3),
                                  hyperpar = list(alpha = alfa_seq, phi = phi_seq),
                                  W = Ws, R = 1)
    
    return(result)
  }
  funs_pred <- lsf.str()[which(lsf.str() != "pred_loop")]
  
  # number of clusters for parallel implementation
  n.core <- parallel::detectCores(logical=F)-1
  
  # starting cluster
  cl <- parallel::makeCluster(n.core)
  doParallel::registerDoParallel(cl)
  
  # parallelized subset computation of GP in different cores
  obj_fit <- foreach::foreach(i = 1:K, .noexport = funs_fit) %dopar% { fit_loop(i,K = K) }
  comb_bps <- spBPS::BPS_combine(obj_fit, K, 1)
  Wbps <- comb_bps$W
  W_list <- comb_bps$W_list
  # R <- 250
  # subset_ind <- sample(1:K, R, T, Wbps)
  # predictions <- foreach::foreach(r = 1:R, .noexport = funs_pred) %dopar% { pred_loop(r) }
  
  # closing cluster
  stopCluster(cl)
  gc()
  
  # return(list("post" = obj_fit,
  #             "wts"  = Wbps,
  #             "pred" = predictions))
  
  return(list("post" = obj_fit,
              "wts"  = Wbps,
              "data" = data_part,
              "alpha" = alfa_seq,
              "phi"   = phi_seq))
  
}

# expectiation MC approximation
ExpVal_jk <- function(k, j, J, L, z_jk, data, alpha_seq, phi_seq, true_data) {
  
    # subsetting data
    K <- k
    crd <- data$crd_list[[k]]
    Y   <- data$Y_list[[k]]
    X   <- data$X_list[[k]]
    
    # prior inizialization
    p <- ncol(X)
    n <- nrow(X)
    q <- ncol(Y)
    priors_list <- list(mu_B = matrix(0, nrow = p, ncol = q),
                        V_r = diag(10, p),
                        Psi = diag(1, q),
                        nu = 3)
    
    # setting the hyperparameters
    par_grid <- expand_grid_cpp(alpha_seq, phi_seq)
    hpar <- list(alpha = par_grid[j, 1], phi = par_grid[j, 2])
    
    # Posterior parameters calculation function
    out <- fit_cpp_MvT(data = list(Y = Y, X = X),
                       priors = priors_list,
                       coords = crd,
                       hyperpar = hpar)
    
    # MC predictive sampling
    X_ul <- cbind(1, matrix(runif((p-1)), ncol = (p-1)))
    crd_ul <- matrix(runif(2), ncol = 2)
    d_ul  <- arma_dist(crd_ul)
    d_uls <- arma_dist(rbind(crd_ul, crd))
    smp <- r_pred_marg_MvT(data = list(Y = Y, X = X),
                           hyperpar = hpar,
                           poster = out, X_u = X_ul, d_u = d_ul, d_us = d_uls, R = L)$Yu |> 
      unlist() |> matrix(ncol = q, byrow = T)
    
    p_jkL <- matrix(0, L, J)
    for (i in 1:J) {
    
    # MC predictive evaluation
    hpar_j <- list(alpha = par_grid[i, 1], phi = par_grid[i, 2])
    p_jkL[, i] <- matrix(sapply(1:L, function(l) d_pred_cpp_MvT(data = list(Y = Y, X = X),
                                                           Y_u = matrix(smp[l,], ncol = q), hyperpar = hpar_j,
                                                           poster = out, X_u = X_ul, d_u = d_ul, d_us = d_uls) ))
    
  }
  
  # true predictive
  crd_t <- true_data$crd
  Y_t   <- true_data$Y
  X_t   <- true_data$X
  true_par <- list(alpha = 0.8, phi = 4)
  out_t <- fit_cpp_MvT(data = list(Y = Y_t, X = X_t),
                       priors = priors_list,
                       coords = crd_t,
                       hyperpar = true_par)
  
  p_tL <- matrix(sapply(1:L, function(a) d_pred_cpp_MvT(data = list(Y = Y_t, X = X_t), Y_u = matrix(smp[a,], ncol = q),
                                                        hyperpar = list(alpha = 0.8, phi = 4),
                                                        poster = out, X_u = X_ul, d_u = d_ul, d_us = d_uls)))
  
  # MC ratios
  s_p_jkL <- p_jkL %*% matrix(z_jk[,k])
  e_jk <- mean(s_p_jkL / p_tL)
  
  return(e_jk)
  
  }

# KL upper bound computation - M data/replications
KL_ub <- function(n, K, J, M = 1) {
  
  # problem dimensions
  q <- 2
  p <- 2
  
  res <- numeric(M)
  pb <- txtProgressBar(0, M)
  for (m in 1:M) {
    
  tic()
  # data generation
  df <- data_gen(n = n)
  post_weights <- calculate_post(Z = df, K = K, sqJ = sqrt(J))
  
  # double BPS weights
  z_jk <- sapply(post_weights$post, function(a) a[[2]])
  w_k <- post_weights$wts
  
  # original data
  crd_t <- df$crd
  Y_t <- df$Y
  X_t <- df$X
  
  # extract post data
  data <- post_weights$data
  alpha_seq <- post_weights$alpha
  phi_seq <- post_weights$phi
  
  # compute MC expectations
  E_JK <- matrix(0, J, K)
  for (k in 1:K) {
    
    for (j in 1:J) {
      
      E_JK[j, k] <- ExpVal_jk(k = k, j = j, J = J, L = 100, z_jk = z_jk,
                              data = data, alpha_seq = alpha_seq, phi_seq = phi_seq, true_data = df)
      
      cat(sprintf("\rComputation done for: j = %d of: k = %d", j, k))
      # flush.console()
      
    }
    
    cat("\n")
    
  }
  
  # final result
  res[m] <- log(prod((sapply(1:K, function(k) t(E_JK) %*% matrix(z_jk[,k]) ) %*% w_k)^(w_k)))
  
  cat(as.numeric(toc()$elapsed), "\n")
  setTxtProgressBar(pb, m)
  cat("\n")
  
  }
  
  cat("\n", M, " replications done for J : ", J, " and K : ", K, "\n")
  return(res)
  
}


# Simulations -------------------------------------------------------------

# Let K vary
K_set <- seq(from = 5, to = 100, by = 5)
kk <- length(K_set)
K_values <- numeric(kk)
pb_k <- txtProgressBar(0, kk, style = 3)
tic()
for (kkk in K_set) {
  k_idx <- which(K_set == kkk)
  K_values[k_idx] <- mean(KL_ub(n = 1000, K = kkk, J = 4, M = 10))
  setTxtProgressBar(pb_k, k_idx)
}
K_time <- toc()
df_K <- data.frame(x = K_set, y = K_values)
K_plot <- ggplot(df_K, aes(x, y)) +
  geom_smooth(method = "loess", se = T, color = "firebrick", linewidth = 1.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Upper bound - variable K",
       x = "K", y = "KL divergence")

# Let J vary
J_set <- seq(from = 2, to = 40, by = 2)
jj <- length(J_set)
J_values <- numeric(jj)
pb_j <- txtProgressBar(0, jj, style = 3)
tic()
for (jjj in J_set) {
  j_idx <- which(J_set == jjj)
  J_values[j_idx] <- mean(KL_ub(n = 1000, K = 5, J = jjj, M = 10))
  setTxtProgressBar(pb_j, j_idx)
}
J_time <- toc()
df_J <- data.frame(x = J_set, y = J_values)
J_plot <- ggplot(df_J, aes(x, y)) +
  geom_smooth(method = "loess", se = T, color = "firebrick", linewidth = 1.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Upper bound - variable J",
       x = "J", y = "KL divergence")


# save results ------------------------------------------------------------

# Size for the mapping
width <- 360*3
height <- 360*1
pointsize <- 12
# Plotting
png("upperbound_simN.png", width = width, height = height, pointsize = pointsize, family = "sans")
cowplot::plot_grid(K_plot, J_plot, nrow = 1)
dev.off()

gc()
# Save the entire environment
results_ub <- list("J" = list("val" = J_values, "set" = J_set, "time" = J_time),
                   "K" = list("val" = K_values, "set" = K_set, "time" = K_time))

rm(list = ls()[which(!(ls() %in% c("results_ub")))])
save.image(file = "upperbound.RData")
# load("upperbound.RData")

# K_time   <- results_ub$K$time
# K_values <- results_ub$K$val
# K_set    <- results_ub$K$set

# J_time   <- results_ub$J$time
# J_values <- results_ub$J$val
# J_set    <- results_ub$J$set
