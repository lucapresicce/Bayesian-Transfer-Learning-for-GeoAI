#' Solver for Bayesian Predictive Stacking of Predictive densities convex optimization problem
#'
#' @param scores [matrix] \eqn{N \times K} of expected predictive density evaluations for the K models considered
#'
#' @return W [matrix] of Bayesian Predictive Stacking weights for the K models considered
#'
#' @importFrom CVXR Variable Maximize Problem solve
#'
#' @export
conv_opt <- function(scores) {
  # library(CVXR, quietly = T)

  # set up minimization problem and solve it
  weights <- Variable( ncol(scores) )
  constraints <- list(weights >= 0, sum(weights) == 1)

  # the constraint for sum up to 1 with positive weights
  f <- Maximize( mean( log( scores %*% weights ) ) )
  problem <- Problem(f, constraints)
  result <- solve(problem, solver = "ECOS_BB") # ECOS, SCS, OSQP

  # return the weights
  W <- if(result$status == "solver_error") {
    matrix(rep(1/ncol(scores), ncol(scores)))
  } else {
    result$getValue(weights)
  }
  return(W)
}
