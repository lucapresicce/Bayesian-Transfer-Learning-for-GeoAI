# Compute the BPS stacking weights for univariate latent spatial regression model
#
# Implement Bayesian predictive stacking for the univariate latent spatial regression models.
#
# @param Y [matrix] \eqn{N \times 1} of sample response variable
# @param X [matrix] \eqn{N \times P} matrix of sample covarites
# @param crd_s [matrix] \eqn{N \times 2} matrix of sample coordinates
# @param Delta [vector] (univariate models) candidate values for hyperparameter \eqn{\delta}
# @param Alfa [vector] (mulrivariate models) candidate values for hyperparameter \eqn{\alpha}
# @param Fi [vector] candidate values for hyperparameter \eqn{\phi}
# @param KCV [boolean] flag to use K-fold cross validation instead of LOOCV (\code{default = FALSE})
# @param K [integer] if \code{KCV = TRUE}, represent the number of desired K-fold
#
# @return A list A list with the following components:
# \item{Grid}{[matrix] models configuration and the associated weights}
# \item{W}{[matrix] stacking weights}
#
# @importFrom CVXR Variable Maximize Problem solve
#
#
#' PredictiveStackingWeights_cpp <- function(Y, X, crd_s, Delta = NULL, Alfa = NULL, Fi, KCV = F, K = 10) {
#'
#'   ## mandatory packages --------------------------------------------------------
#'   # library(CVXR, quietly = T)
#'
#'   ## evaluate the loocv predictive density -------------------------------------
#'   q <- ncol(Y)
#'   if (q > 1) {
#'
#'     ## control on Alfa
#'     if (is.null(Alfa)) stop("Set of values for alpha (Alfa) is missing")
#'
#'     ## dimensions for priors
#'     q <- ncol(Y)
#'     p <- ncol(X)
#'
#'     ## evaluate the loocv/kcv predictive density
#'     out <- models_dens_latent(data = list(Y = Y, X = X), coords = crd_s,
#'                               priors = list(mu_B = matrix(0, nrow = p, ncol = q),
#'                                             V_r = diag(10, p),
#'                                             Psi = diag(1, q),
#'                                             nu = 3),
#'                               hyperpar = list(alpha = Alfa,
#'                                               phi = Fi),
#'                               useKCV = KCV, K = K)
#'   } else {
#'
#'     ## control on Delta
#'     if (is.null(Delta)) stop("Set of values for delta (Delta) is missing")
#'
#'     ## dimensions for priors
#'     q <- ncol(Y)
#'     p <- ncol(X)
#'
#'     ## evaluate the loocv/kcv predictive density
#'     out <- models_dens(data = list(Y = Y, X = X), coords = crd_s,
#'                        priors = list(mu_b = matrix(rep(0, p)),
#'                                      V_b = diag(10, p),
#'                                      a = 2,
#'                                      b = 2),
#'                        hyperpar = list(delta = Delta,
#'                                        phi = Fi),
#'                        useKCV = KCV, K = K)
#'   }
#'
#'   ## solve the convex optimization problem -------------------------------------
#'
#'   # declare variable
#'   scores <- out
#'   weights <- Variable( ncol(scores) )
#'
#'   # set up minimization problem and solve it
#'   constraints <- list(weights >= 0, sum(weights) == 1)
#'   # the constraint for sum up to 1 with positive weights
#'   f <- Maximize( mean( log( scores %*% weights ) ) )
#'   problem <- Problem(f, constraints)
#'   result <- solve(problem, solver = "ECOS_BB") # ECOS, SCS, OSQP
#'
#'   # define the stacking predictive distribution obtained
#'   w_hat <- result$getValue(weights)
#'
#'   ## function return
#'   if (q > 1) {
#'     grid <- expand_grid_cpp(Alfa, Fi)
#'     dfw <- as.matrix(data.frame(
#'       "W" = round(w_hat, 3), "alfa" = grid[, 1], "phi" = grid[, 2]
#'     ))
#'   } else {
#'     grid <- expand_grid_cpp(Delta, Fi)
#'     dfw <- as.matrix(data.frame(
#'       "W" = round(w_hat, 3), "delta" = grid[, 1], "phi" = grid[, 2]
#'     ))
#'   }
#'
#'   return(list("Grid" = dfw, "W" = w_hat))
#'
#' }
#'
#'
# Compute the BPS stacking weights for univariate latent spatial regression model
#
# @param Y [matrix] \eqn{N \times 1} of sample response variable
# @param X [matrix] \eqn{N \times P} matrix of sample covarites
# @param crd_s [matrix] \eqn{N \times 2} matrix of sample coordinates
# @param Delta [vector] (univariate models) candidate values for hyperparameter \eqn{\delta}
# @param Fi [vector] candidate values for hyperparameter \eqn{\phi}
# @param K [integer] if \code{KCV = TRUE}, represent the number of desired K-fold
#
# @return A list with the following components:
# \item{Grid}{[matrix] models configuration and the associated weights}
# \item{W}{[matrix] stacking weights}
#
# @importFrom CVXR Variable Maximize Problem solve
#
#
#' # Compute the stacking weights for univariate spatial regression (latent model)
#' PredictiveStackingWeights_cpp2 <- function(Y, X, crd_s, Delta, Fi, K = 10) {
#'
#'   ## mandatory packages --------------------------------------------------------
#'   # library(CVXR, quietly = T)
#'
#'   ## evaluate the loocv predictive density
#'   p <- ncol(X)
#'   out <- models_dens2(data = list(Y = Y, X = X), coords = crd_s,
#'                       priors = list(mu_b = matrix(rep(0, p)),
#'                                     V_b = diag(10, p),
#'                                     a = 2,
#'                                     b = 2),
#'                       hyperpar = list(delta = Delta,
#'                                       phi = Fi), K = K)
#'
#'   ## solve the convex optimization problem -------------------------------------
#'
#'   # declare variable
#'   scores <- out$out
#'   weights <- Variable( ncol(scores) )
#'
#'   # set up minimization problem and solve it
#'   constraints <- list(weights >= 0, sum(weights) == 1)
#'   # the constraint for sum up to 1 with positive weights
#'   f <- Maximize( mean( log( scores %*% weights ) ) )
#'   problem <- Problem(f, constraints)
#'   result <- solve(problem, solver = "ECOS_BB") # ECOS, SCS, OSQP
#'
#'   # define the stacking predictive distribution obtained
#'   w_hat <- result$getValue(weights)
#'
#'   ## function return
#'   grid <- expand_grid_cpp(Delta, Fi)
#'   dfw <- as.matrix(data.frame(
#'     "W" = round(w_hat, 3), "delta" = grid[, 1], "phi" = grid[, 2]
#'   ))
#'   return(list("Grid" = dfw, "W" = w_hat, "beta" = out$beta, "sigma" = out$sigma))
#'
#' }
#'
#'
# Using Bayesian Predictive Stacking to combine subset models
#
# @param fit_list [list] of fitted models for all the subsets
#
# @return A list with the following components:
# \item{W}{[matrix] stacking weights of all model configurations}
# \item{W_lsit}{[list] stackign weights for each model configurations}
#
#' BPScombine <- function(fit_list) {
#'
#'   # number of subsets
#'   K <- length(fit_list)
#'
#'   # weights list
#'   W_list <- list()
#'   for(i in 1:K){
#'     W_list[[i]] <- fit_list[[i]][[2]]
#'     attr(W_list[[i]], "class") <- "matrix"
#'   }
#'
#'   # epd list
#'   out3 <- lapply(1:K, function(i){fit_list[[i]][[1]]})
#'   out3 <- do.call(rbind, out3)
#'   # ind_0 <- which(apply(out3, 1, sum) !=0 )
#'   # out3 <- out3[ind_0,]
#'   out4 <- sapply(1:K, function(k){
#'     out3 %*% W_list[[k]]
#'   })
#'
#'   # convex optimization
#'   Wbps <- conv_opt(scores = out4)
#'
#'   # remove small weights - (3.1) of "Scalable and Robust Bayesian Inference via the Median Posterior"
#'   ind_W <- which(Wbps > 1/(2*length(Wbps)))
#'   Wbps[-ind_W,] <- 0
#'
#'   return(list("W" = Wbps,
#'               "W_list" = W_list))
#' }
