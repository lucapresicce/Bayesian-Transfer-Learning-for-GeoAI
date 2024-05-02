// #include <RcppArmadillo.h>
// #include "code.h"
// #include "utilsC.h"
//
// // [[Rcpp::depends(RcppArmadillo)]]
//
// using namespace Rcpp;
// using namespace arma;
//
//
// // ##################################################################################################################################################
// // BACKUP FUNTCIONS - UNIVARIATE ####################################################################################################################
// // ##################################################################################################################################################
//
// //' Compute the parameters for the posteriors distribution of \eqn{\beta} and \eqn{\Sigma} (i.e. updated parameters)
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //'
// //' @return [list] posterior update parameters
// //'
// // [[Rcpp::export]]
// List fit_cpp_optimized(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar) {
//
//   // Unpack data and priors
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   arma::vec mu_b = as<arma::vec>(priors["mu_b"]);
//   arma::mat V_b = as<arma::mat>(priors["V_b"]);
//   double b = as<double>(priors["b"]);
//   double a = as<double>(priors["a"]);
//   double delta = as<double>(hyperpar["delta"]);
//   double phi = as<double>(hyperpar["phi"]);
//
//   int n = Y.n_rows;
//
//   arma::mat d_s = arma_dist(coords);
//   arma::mat Rphi_s = exp(-phi * d_s);
//
//   // Precompute some reusable values
//   double d = (1 / delta);
//   arma::mat tX = trans(X);
//   arma::mat iV_b = arma::inv(V_b);
//   arma::mat iR_s = arma::inv(Rphi_s);
//
//   // Compute posterior updating
//   arma::mat iM_B = d * tX * X + iV_b;
//   arma::mat iM_BW = d * tX;
//   arma::mat iM_WB = trans(iM_BW);
//   arma::mat iM_W = iR_s + (d * eye<arma::mat>(n, n));
//
//   arma::mat iM_star1 = join_horiz(iM_B, iM_BW);
//   arma::mat iM_star2 = join_horiz(iM_WB, iM_W);
//   arma::mat iM_star = join_vert( iM_star1, iM_star2);
//   arma::mat M_star = arma::inv(iM_star);
//
//   arma::vec M = join_vert( (d * tX * Y) + (iV_b * mu_b) , d * Y );
//   arma::vec gamma_hat = M_star * M;
//
//   double dYY = as_scalar(d * trans(Y) * Y);
//   double mbVbmb = as_scalar(trans(mu_b) * iV_b * mu_b);
//   double bb = dYY + mbVbmb + as_scalar(trans(gamma_hat) * iM_star * gamma_hat) - as_scalar(2 * trans(gamma_hat) * M);
//
//   double b_star = b + 0.5 * as_scalar(bb);
//   double a_star = a + (n/2);
//
//   // Return results as an R list
//   return List::create(Named("M_star") = M_star,
//                       Named("gamma_hat") = gamma_hat,
//                       Named("b_star") = b_star,
//                       Named("a_star") = a_star,
//                       Named("iRphi_s") = iR_s);
// }
//
// //' Draw from the conditional posterior predictive for a set of unobserved covariates
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param iRphi_s [matrix] inverse of the sample correlation matrix
// //' @param d_u [matrix] unobserved instances distance matrix
// //' @param d_us [matrix] cross-distance between unobserved and observed instances matrix
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param poster [list] output from \code{fit_cpp} function
// //'
// //' @return [list] posterior predictive samples
// //'
// // [[Rcpp::export]]
// List r_pred_MC(const List& data, const arma::mat& X_u, const arma::mat& iRphi_s, const arma::mat& d_u, const arma::mat& d_us, const List& hyperpar, const List& poster) {
//
//   // Unpack data, posterior sample and hyperparameters
//   arma::mat beta = as<arma::mat>(poster["Betas"]);
//   arma::vec sigma = as<arma::vec>(poster["Sigmas"]);
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   double delta = as<double>(hyperpar["delta"]);
//   double phi = as<double>(hyperpar["phi"]);
//
//   // extract info from data
//   int m = X_u.n_rows;
//   int p = X_u.n_cols;
//   int n = d_us.n_rows-m;
//   int R = beta.n_rows;
//
//   // covariance matrices
//   // arma::mat iR_s = arma::inv(Rphi_s);
//   arma::mat Rphi_u = exp(-phi * d_u);
//   arma::mat Rphi_us = exp(-phi * d_us.submat(0, m, m-1, m+n-1));
//
//   // environment
//   Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
//   Rcpp::Function rmNorm_R = mniw["rmNorm"];
//   //
//   // initialize return objects
//   arma::mat Z_u(m, R);
//   arma::mat Y_u(m, R);
//
//   // compute reusable
//   arma::mat RiR = Rphi_us * iRphi_s;
//   arma::mat V_z = Rphi_u - RiR * trans(Rphi_us);
//
//   for (int r = 0; r < R; ++r) {
//
//     // unpack posterior sample
//     arma::vec gamma_hat_r = trans(beta.row(r));
//     arma::vec b = gamma_hat_r.subvec(0, p - 1);
//     arma::vec gamma_r = gamma_hat_r.subvec(p, gamma_hat_r.n_elem - 1);
//     double s = sigma(r);
//
//
//     // predictive conjugate parameters
//     // arma::mat mu_z = Rphi_us * iR_s * gamma_r;
//     // arma::mat V_z = Rphi_u - Rphi_us * iR_s * trans(Rphi_us);
//     arma::mat mu_z = RiR * gamma_r;
//     Z_u.col(r) = as<arma::vec>(rmNorm_R(Named("n", 1), Named("mu", trans(mu_z)), Named("Sigma", s * V_z)));
//
//     arma::mat mu_y = X_u * b + Z_u.col(r);
//     arma::mat V_y = (s * delta) * eye<arma::mat>(m, m);
//     Y_u.col(r) = as<arma::vec>(rmNorm_R(Named("n", 1), Named("mu", trans(mu_y)), Named("Sigma", V_y)));
//
//   }
//
//   return List::create(Named("Z_u") = Z_u,
//                       Named("Y_u") = Y_u);
// }
//
//
// //' Evaluate the density of a set of unobserved response with respect to the conditional posterior predictive
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param Y_u [matrix] unobserved instances response matrix
// //' @param iRphi_s [matrix] inverse of the sample correlation matrix
// //' @param d_u [matrix] unobserved instances distance matrix
// //' @param d_us [matrix] cross-distance between unobserved and observed instances matrix
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param poster [list] output from \code{fit_cpp} function
// //'
// //' @return [vector] posterior predictive density evaluations
// //'
// // [[Rcpp::export]]
// arma::vec d_pred_MC(const List& data, const arma::mat& X_u, const arma::vec& Y_u, const arma::mat& iRphi_s, const arma::mat& d_u, const arma::mat& d_us, const List& hyperpar, const List& poster) {
//
//   // Unpack data, posterior sample and hyperparameters
//   arma::mat beta = as<arma::mat>(poster["Betas"]);
//   arma::vec sigma = as<arma::vec>(poster["Sigmas"]);
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   double delta = as<double>(hyperpar["delta"]);
//   double phi = as<double>(hyperpar["phi"]);
//
//   // extract info from data
//   int m = X_u.n_rows;
//   int p = X_u.n_cols;
//   int n = d_us.n_rows-m;
//   int R = beta.n_rows;
//
//   // covariance matrices
//   // arma::mat iR_s = arma::inv(Rphi_s);
//   arma::mat Rphi_u = exp(-phi * d_u);
//   arma::mat Rphi_us = exp(-phi * d_us.submat(0, m, m-1, m+n-1));
//
//   // environment
//   Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
//   Rcpp::Function rmNorm_R = mniw["rmNorm"];
//   Rcpp::Function dmNorm_R = mniw["dmNorm"];
//
//
//   // initialize return objects
//   arma::vec P_u(R);
//
//   // compute reusable
//   arma::mat RiR = Rphi_us * iRphi_s;
//   arma::mat V_z = Rphi_u - RiR * trans(Rphi_us);
//
//   for (int r = 0; r < R; ++r) {
//
//     // unpack posterior sample
//     arma::vec gamma_hat_r = trans(beta.row(r));
//     arma::vec b = gamma_hat_r.subvec(0, p - 1);
//     arma::vec gamma_r = gamma_hat_r.subvec(p, gamma_hat_r.n_elem - 1);
//     double s = sigma(r);
//
//     // predictive conjugate parameters
//     // arma::mat mu_z = Rphi_us * iR_s * gamma_r;
//     // arma::mat V_z = Rphi_u - Rphi_us * iR_s * trans(Rphi_us);
//     arma::mat mu_z = RiR * gamma_r;
//     arma::vec Z_u = as<arma::vec>(rmNorm_R(Named("n", 1), Named("mu", trans(mu_z)), Named("Sigma", s * V_z)));
//
//     arma::mat mu_y = X_u * b + Z_u;
//     arma::mat V_y = (s * delta) * eye<arma::mat>(m, m);
//     P_u(r) = as<double>(dmNorm_R(Named("x", trans(Y_u)), Named("mu", trans(mu_y)), Named("Sigma", V_y)));
//
//   }
//
//   return P_u;
// }
//
//
// //' Compute the LOOCV of the density evaluations for fixed values of the hyperparameters
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //'
// //' @return [vector] posterior predictive density evaluations
// //'
// // [[Rcpp::export]]
// arma::vec dens_loocv_MC(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar) {
//
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   int n = Y.n_rows;
//   arma::vec predictions(n);
//
//   for (int i = 0; i < n; i++) {
//
//     // Rcout << "i: " << i << std::endl;
//
//     // extract the training data (excluding the i-th row) (create a mask to exclude the i-th row)
//     arma::mat Ytr, Xtr, coords_i;
//     if (i == 0) {
//       // Exclude the 0-th row
//       Ytr = Y.rows(1, n - 1);
//       Xtr = X.rows(1, n - 1);
//       coords_i = coords.rows(1, n - 1);
//     } else if (i == 1) {
//       // Exclude the 1-th row
//       Ytr = join_vert(Y.row(0), Y.rows(2, n - 1));
//       Xtr = join_vert(X.row(0), X.rows(2, n - 1));
//       coords_i = join_vert(coords.row(0), coords.rows(2, n - 1));
//     } else if (i == n - 1) {
//       // Exclude the last row
//       Ytr = Y.rows(0, n - 2);
//       Xtr = X.rows(0, n - 2);
//       coords_i = coords.rows(0, n - 2);
//     } else {
//       // Exclude the i-th row for i > 1
//       Ytr = join_vert(Y.rows(0, i - 1), Y.rows(i + 1, n - 1));
//       Xtr = join_vert(X.rows(0, i - 1), X.rows(i + 1, n - 1));
//       coords_i = join_vert(coords.rows(0, i - 1), coords.rows(i + 1, n - 1));
//     }
//     List data_i = List::create(
//       Named("Y") = Ytr,
//       Named("X") = Xtr);
//
//     // extract the test data
//     arma::mat crd_i = coords.row(i);
//     arma::mat X_i = X.row(i);
//     arma::mat Y_i = Y.row(i);
//
//     // Fit your model on the training data
//     List poster_i = fit_cpp(data_i, priors, coords_i, hyperpar);
//
//     // posterior draws
//     List post_i = post_draws(poster_i, 1);
//
//     // evaluate predictive density
//     arma::mat Rphi_s = as<arma::mat>(poster_i["Rphi_s"]);
//     arma::mat d_i = arma_dist(crd_i);
//     arma::mat crd_is = join_cols(crd_i, coords_i);
//     arma::mat d_is = arma_dist(crd_is);
//
//     double dens = as_scalar(d_pred_cpp(data_i, X_i, Y_i, Rphi_s, d_i, d_is, hyperpar, post_i));
//
//     // Store the prediction
//     predictions[i] = dens;
//   }
//
//   return predictions;
//
// }
//
//
// //' Compute the KCV of the density evaluations for fixed values of the hyperparameters
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param K [integer] number of folds
// //' @param g [integer] number of draws for MC density approximations
// //'
// //' @return [vector] posterior predictive density evaluations
// //'
// // [[Rcpp::export]]
// arma::vec dens_kcv_MC(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, const int& K, const int& g) { //, const int& seed) {
//
//   // unpack data
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   int n = Y.n_rows;
//   arma::vec predictions(n, arma::fill::zeros);
//
//   // Create a random permutation of indices from 1 to K
//   arma::vec p = arma::vec(K, arma::fill::ones) / K;
//   arma::uvec foldIndices = sample_index(K, n, p);
//
//   for (int k = 0; k < K; k++) {
//
//     // Define the indices for the current fold
//     arma::uvec testSet = find(foldIndices == (k + 1)); // Find indices that match the current fold number
//
//     // Create a boolean vector to identify the training set for this fold
//     arma::uvec trainSet = find(foldIndices != (k + 1)); // Find indices that do not match the current fold number
//
//     // Extract the training data for the current fold
//     arma::mat Ytr = Y.rows(trainSet);
//     arma::mat Xtr = X.rows(trainSet);
//     arma::mat coords_tr = coords.rows(trainSet);
//
//     List data_tr = List::create(
//       Named("Y") = Ytr,
//       Named("X") = Xtr);
//
//     // Extract the test data for the current fold
//     arma::mat Y_test = Y.rows(testSet);
//     arma::mat X_test = X.rows(testSet);
//     arma::mat coords_test = coords.rows(testSet);
//
//     // Fit your model on the training data
//     List poster_k = fit_cpp(data_tr, priors, coords_tr, hyperpar);
//     arma::mat iRphi_k = as<arma::mat>(poster_k["iRphi_s"]);
//
//     // posterior draws
//     List post_k = post_draws(poster_k, g);
//
//     // evaluate predictive density for the test set of the current fold
//     for (uword i = 0; i < testSet.n_elem; i++) {
//
//       arma::mat crd_i = coords_test.row(i);
//       arma::mat X_i = X_test.row(i);
//       arma::mat Y_i = Y_test.row(i);
//       arma::mat d_i = arma_dist(crd_i);
//
//       // Rcout << crd_i << std::endl;
//       // Rcout << coords_tr << std::endl;
//
//       arma::mat crd_is = join_cols(crd_i, coords_tr);
//       arma::mat d_is = arma_dist(crd_is);
//
//       double dens = mean(d_pred_cpp(data_tr, X_i, Y_i, iRphi_k, d_i, d_is, hyperpar, post_k));
//       predictions(testSet(i)) = dens;
//     }
//   }
//
//   return predictions;
//
// }
//
//
// //' Return the CV predictive density evaluations for all the model combinations
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param useKCV [boolean] if \code{TRUE} K-fold cross validation is used instead of LOOCV (no \code{default})
// //' @param K [integer] number of folds
// //' @param g [integer] number of draws for MC density approximations
// //'
// //' @return [matrix] posterior predictive density evaluations (each columns represent a different model)
// //'
// // [[Rcpp::export]]
// arma::mat models_dens_MC(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, const bool& useKCV, const int& K, const int& g) {
//
//   // build the grid of hyperparameters
//   arma::vec Delta = hyperpar["delta"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Delta, Fi);
//   int k = Grid.n_rows;
//
//   arma::mat out;
//
//   for(int j = 0; j < k; j++) {
//
//     // identify the model
//     arma::rowvec hpar = Grid.row(j);
//     double delt = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("delta") = delt,
//       Named("phi") = fi);
//
//     // Call the appropriate function based on the 'useKCV' argument
//     arma::vec out_j;
//     if (useKCV) {
//       out_j = dens_kcv(data, priors, coords, hmod, K, g);
//     } else {
//       out_j = dens_loocv(data, priors, coords, hmod);
//     }
//
//     out =  join_horiz(out, out_j);
//
//   }
//
//   return out;
// }
//
//
// //' Compute the KCV of the density evaluations for fixed values of the hyperparameters
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param K [integer] number of folds
// //'
// //' @return [list] posterior predictive density evaluations
// //'
// // [[Rcpp::export]]
// List dens_kcv2(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, const int& K) { //, const int& seed) {
//
//   // unpack data
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   int n = Y.n_rows;
//   int p = X.n_cols;
//   arma::vec predictions(n, arma::fill::zeros);
//   arma::mat all_beta(K, p, arma::fill::zeros);
//   arma::vec all_sigma(K, arma::fill::zeros);
//
//   // Create a random permutation of indices from 1 to K
//   arma::vec pr = arma::vec(K, arma::fill::ones) / K;
//   arma::uvec foldIndices = sample_index(K, n, pr);
//
//   for (int k = 0; k < K; k++) {
//
//     // Define the indices for the current fold
//     arma::uvec testSet = find(foldIndices == (k + 1)); // Find indices that match the current fold number
//
//     // Create a boolean vector to identify the training set for this fold
//     arma::uvec trainSet = find(foldIndices != (k + 1)); // Find indices that do not match the current fold number
//
//     // Extract the training data for the current fold
//     arma::mat Ytr = Y.rows(trainSet);
//     arma::mat Xtr = X.rows(trainSet);
//     arma::mat coords_tr = coords.rows(trainSet);
//
//     List data_tr = List::create(
//       Named("Y") = Ytr,
//       Named("X") = Xtr);
//
//     // Extract the test data for the current fold
//     arma::mat Y_test = Y.rows(testSet);
//     arma::mat X_test = X.rows(testSet);
//     arma::mat coords_test = coords.rows(testSet);
//
//     // Fit your model on the training data
//     List poster_k = fit_cpp(data_tr, priors, coords_tr, hyperpar);
//     arma::mat Rphi_k = as<arma::mat>(poster_k["Rphi_s"]);
//
//     // posterior draws
//     List post_k = post_draws(poster_k, 1);
//
//     // Unpack posterior sample
//     arma::rowvec beta = as<arma::rowvec>(post_k["Betas"]);
//     arma::vec sigma = as<arma::vec>(post_k["Sigmas"]);
//
//     // Store the results for the current fold
//     all_beta.row(k) = beta.subvec(0, p - 1);
//     all_sigma(k) = sigma(0);
//
//     // evaluate predictive density for the test set of the current fold
//     for (uword i = 0; i < testSet.n_elem; i++) {
//
//       arma::mat crd_i = coords_test.row(i);
//       arma::mat X_i = X_test.row(i);
//       arma::mat Y_i = Y_test.row(i);
//       arma::mat d_i = arma_dist(crd_i);
//
//       // Rcout << crd_i << std::endl;
//       // Rcout << coords_tr << std::endl;
//
//       arma::mat crd_is = join_cols(crd_i, coords_tr);
//       arma::mat d_is = arma_dist(crd_is);
//
//       double dens = as_scalar(d_pred_cpp(data_tr, X_i, Y_i, Rphi_k, d_i, d_is, hyperpar, post_k));
//       predictions(testSet(i)) = dens;
//     }
//   }
//
//   // Create a list to hold both predictions and post_k outputs
//   List results;
//   results["predictions"] = predictions;
//   results["beta"] = all_beta;
//   results["sigma"] = all_sigma;
//
//   return results;
// }
//
//
// //' Return the CV predictive density evaluations for all the model combinations
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param K [integer] number of folds
// //'
// //' @return [matrix] posterior predictive density evaluations (each columns represent a different model)
// //'
// // [[Rcpp::export]]
// List models_dens2(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, const int& K) {
//
//   // build the grid of hyperparameters
//   arma::vec Delta = hyperpar["delta"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Delta, Fi);
//   int k = Grid.n_rows;
//
//   arma::mat out;
//   List beta(k);
//   List sigma(k);
//
//   for(int j = 0; j < k; j++) {
//
//     // identify the model
//     arma::rowvec hpar = Grid.row(j);
//     double delt = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("delta") = delt,
//       Named("phi") = fi);
//
//     // Call the KCV function based on K folds
//     List res = dens_kcv2(data, priors, coords, hmod, K);
//     arma::vec out_j = res["predictions"];
//     arma::mat beta_j = res["beta"];
//     arma::vec sigma_j = res["sigma"];
//
//     out =  join_horiz(out, out_j);
//     beta(j) =  beta_j;
//     sigma(j) =  sigma_j;
//
//   }
//
//   List results;
//   results["out"] = out;
//   results["beta"] = beta;
//   results["sigma"] = sigma;
//
//   return results;
// }
//
//
// //' Compute the BPS weights by convex optimization
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param K [integer] number of folds
// //'
// //' @return [matrix] posterior predictive density evaluations (each columns represent a different model)
// //'
// // [[Rcpp::export]]
// arma::mat BPSweights_cpp2(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, int K) {
//
//   // compute predictive density evaluations
//   arma::mat out = models_dens(data, priors, coords, hyperpar, true, K);
//
//   // compute the weights
//   arma::mat weights = as<arma::mat>(CVXR_opt(out));
//
//   // return the list (BPS weights over model configurations)
//   arma::vec Delta = hyperpar["delta"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Delta, Fi);
//   arma::mat res = join_horiz(weights, Grid);
//
//   return res;
// }
//
//
// //' Compute the BPS spatial prediction given a set of stacking weights
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //'
// //' @return [list] BPS posterior predictive samples
// //'
// // [[Rcpp::export]]
// List BPS_post(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R) {
//
//   arma::mat Z_pred;
//   arma::mat Y_pred;
//   arma::mat Betas;
//   arma::vec Sigmas(R);
//
//   // compute distance matrices
//   arma::mat d_u = arma_dist(crd_u);
//   arma::mat crd_us = join_cols(crd_u, coords);
//   arma::mat d_us = arma_dist(crd_us);
//
//   for(int r = 0; r < R; r++) {
//
//     // build the grid of hyperparameters
//     arma::vec Delta = hyperpar["delta"];
//     arma::vec Fi = hyperpar["phi"];
//     arma::mat Grid = expand_grid_cpp(Delta, Fi);
//     int k = Grid.n_rows;
//
//     // sample the model
//     arma::uvec kmod = sample_index(k, 1, W);
//     arma::uword k_mod = kmod(0);
//
//     // identify the k-th model
//     arma::rowvec hpar = Grid.row(k_mod);
//     double delt = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("delta") = delt,
//       Named("phi") = fi);
//
//     // fit your model on the training data
//     List poster = fit_cpp(data, priors, coords, hmod);
//
//     // posterior draws
//     List post = post_draws(poster, 1);
//
//     // save posterior samples
//     arma::mat beta = as<arma::mat>(post["Betas"]);
//     arma::vec sigma = as<arma::vec>(post["Sigmas"]);
//     Betas =  join_vert(Betas, beta);
//     Sigmas(r) =  sigma(0);
//
//     // draw from conditional posterior predictive
//     arma::mat iRphi_s = as<arma::mat>(poster["iRphi_s"]);
//     List pred_R = r_pred_cpp(data, X_u, iRphi_s, d_u, d_us, hmod, post);
//
//     arma::vec Z_pred_r = as<arma::mat>(pred_R["Z_u"]);
//     Z_pred =  join_horiz(Z_pred, Z_pred_r);
//
//     arma::vec Y_pred_r = as<arma::mat>(pred_R["Y_u"]);
//     Y_pred =  join_horiz(Y_pred, Y_pred_r);
//
//   }
//
//   // return pred;
//   return List::create(Named("Z_hat") = Z_pred,
//                       Named("Y_hat") = Y_pred,
//                       Named("Betas") = Betas,
//                       Named("Sigmas") = Sigmas);
//
// }
//
//
// //' Compute the BPS spatial prediction given a set of stacking weights
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //'
// //' @return [list] BPS posterior predictive samples
// //'
// // [[Rcpp::export]]
// List BPS_SpatialPrediction_cpp2(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R) {
//
//   // build the grid of hyperparameters
//   arma::vec Delta = hyperpar["delta"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Delta, Fi);
//   int k = Grid.n_rows;
//
//   List Z_pred(k);
//   List Y_pred(k);
//
//   for(int i = 0; i < k; i++) {
//
//     // identify the i-th set of hyperparameters
//     arma::rowvec hpar = Grid.row(i);
//     double delt = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("delta") = delt,
//       Named("phi") = fi);
//
//     // fit your model on the training data
//     List poster = fit_cpp(data, priors, coords, hmod);
//
//     // posterior draws
//     // List post = post_draws(poster, 1);
//     List post = post_draws(poster, R);
//
//     // compute distance matrices
//     arma::mat iRphi_s = as<arma::mat>(poster["iRphi_s"]);
//     arma::mat d_u = arma_dist(crd_u);
//     arma::mat crd_us = join_cols(crd_u, coords);
//     arma::mat d_us = arma_dist(crd_us);
//
//     // draw from conditional posterior predictive
//     List pred_R = r_pred_cpp(data, X_u, iRphi_s, d_u, d_us, hmod, post);
//
//     // Extract and assign for the given set of indices
//     arma::mat Z_pred_r = as<arma::mat>(pred_R["Z_u"]);
//     Z_pred(i) = Z_pred_r;
//
//     // vec Y_pred_r = as<mat>(pred_R["Y_u"]);
//     arma::mat Y_pred_r = as<arma::mat>(pred_R["Y_u"]);
//     Y_pred(i) = Y_pred_r;
//   }
//
//   // sample the models outside the loop
//   arma::uvec kmod = sample_index(k, R, W);
//
//
//   arma::mat Z_hat(X_u.n_rows, R, arma::fill::zeros);
//   arma::mat Y_hat(X_u.n_rows, R, arma::fill::zeros);
//
//   // Imputation loop
//   for(int r = 0; r < R; r++) {
//     arma::uword k_mod = kmod(r);
//
//     arma::mat Z_mod = Z_pred(k_mod);
//     Z_hat.col(r) = Z_mod.col(r);
//
//     arma::mat Y_mod = Y_pred(k_mod);
//     Y_hat.col(r) = Y_mod.col(r);
//   }
//
//   // return as a List
//   return List::create(Named("Z_hat") = Z_hat,
//                       Named("Y_hat") = Y_hat);
// }
//
//
// //' Compute the BPS spatial prediction given a set of stacking weights
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //'
// //' @return [list] BPS posterior predictive samples
// //'
// // [[Rcpp::export]]
// List fast_BPSpred2(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R) {
//
//   // build the grid of hyperparameters
//   arma::vec Delta = hyperpar["delta"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Delta, Fi);
//   int k = Grid.n_rows;
//
//   List Z_pred(k);
//   List Y_pred(k);
//
//   for(int i = 0; i < k; i++) {
//
//     // identify the i-th set of hyperparameters
//     arma::rowvec hpar = Grid.row(i);
//     double delt = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("delta") = delt,
//       Named("phi") = fi);
//
//     // fit your model on the training data
//     List poster = fit_cpp(data, priors, coords, hmod);
//
//     // posterior draws
//     // List post = post_draws(poster, 1);
//     List post = post_draws(poster, R);
//
//     // compute distance matrices
//     arma::mat iRphi_s = as<mat>(poster["iRphi_s"]);
//     arma::mat d_u = arma_dist(crd_u);
//     arma::mat crd_us = join_cols(crd_u, coords);
//     arma::mat d_us = arma_dist(crd_us);
//
//     // draw from conditional posterior predictive
//     List pred_R = r_pred_cpp(data, X_u, iRphi_s, d_u, d_us, hmod, post);
//
//     // Extract and assign for the given set of indices
//     arma::mat Z_pred_r = as<arma::mat>(pred_R["Z_u"]);
//     Z_pred(i) = Z_pred_r;
//
//     arma::mat Y_pred_r = as<arma::mat>(pred_R["Y_u"]);
//     Y_pred(i) = Y_pred_r;
//   }
//
//   // sample the models outside the loop
//   arma::uvec kmod = sample_index(k, R, W);
//
//
//   arma::mat Z_hat(X_u.n_rows, R, arma::fill::zeros);
//   arma::mat Y_hat(X_u.n_rows, R, arma::fill::zeros);
//
//   // Imputation loop
//   for(int r = 0; r < R; r++) {
//     arma::uword k_mod = kmod(r);
//
//     arma::mat Z_mod = Z_pred(k_mod);
//     Z_hat.col(r) = Z_mod.col(r);
//
//     arma::mat Y_mod = Y_pred(k_mod);
//     Y_hat.col(r) = Y_mod.col(r);
//   }
//
//   // return as a List
//   return List::create(Named("Z_hat") = Z_hat,
//                       Named("Y_hat") = Y_hat);
// }
//
//
// //' Compute the BPS spatial prediction given a set of stacking weights
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //'
// //' @return [list] BPS posterior predictive samples
// //'
// // [[Rcpp::export]]
// List fast_BPSpred3(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R) {
//
//   // build the grid of hyperparameters
//   arma::vec Delta = hyperpar["delta"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Delta, Fi);
//   int k = Grid.n_rows;
//
//   // sample the models from which predict
//   arma::uvec kmod = sample_index(k, R, W);
//   arma::uvec uniqmod = arma::unique(kmod);
//   int l = uniqmod.size();
//
//   List Z_pred(l);
//   List Y_pred(l);
//
//   for(int i = 0; i < l; i++) {
//
//     // select the model
//     arma::uword k_mod = uniqmod(i);
//
//     // identify the i-th set of hyperparameters
//     arma::rowvec hpar = Grid.row(k_mod);
//     double delt = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("delta") = delt,
//       Named("phi") = fi);
//
//     // fit your model on the training data
//     List poster = fit_cpp(data, priors, coords, hmod);
//
//     // posterior draws
//     // List post = post_draws(poster, 1);
//     List post = post_draws(poster, R);
//
//     // compute distance matrices
//     arma::mat iRphi_s = as<mat>(poster["iRphi_s"]);
//     arma::mat d_u = arma_dist(crd_u);
//     arma::mat crd_us = join_cols(crd_u, coords);
//     arma::mat d_us = arma_dist(crd_us);
//
//     // draw from conditional posterior predictive
//     List pred_R = r_pred_cpp(data, X_u, iRphi_s, d_u, d_us, hmod, post);
//
//     // Extract and assign for the given set of indices
//     arma::mat Z_pred_r = as<arma::mat>(pred_R["Z_u"]);
//     Z_pred(i) = Z_pred_r;
//
//     // vec Y_pred_r = as<mat>(pred_R["Y_u"]);
//     arma::mat Y_pred_r = as<arma::mat>(pred_R["Y_u"]);
//     Y_pred(i) = Y_pred_r;
//   }
//
//   // initialize result containers
//   arma::mat Z_hat(X_u.n_rows, R, arma::fill::zeros);
//   arma::mat Y_hat(X_u.n_rows, R, arma::fill::zeros);
//
//   // Imputation loop
//   for(int r = 0; r < R; r++) {
//
//     // sample the models
//     arma::uword r_mod = kmod(r);
//     arma::uvec j_mod = arma::find(uniqmod == r_mod);
//
//     // assign prediction
//     arma::mat Z_mod = Z_pred(j_mod(0));
//     Z_hat.col(r) = Z_mod.col(r);
//
//     // assign prediction
//     arma::mat Y_mod = Y_pred(j_mod(0));
//     Y_hat.col(r) = Y_mod.col(r);
//   }
//
//   // return as a List
//   return List::create(Named("Z_hat") = Z_hat,
//                       Named("Y_hat") = Y_hat);
// }
//
//
// //' Perform prediction for ASMK models - loop over prediction set
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //' @param J [integer] number of desired partition of prediction set
// //'
// //' @return [list] BPS posterior predictive samples
// //'
// // [[Rcpp::export]]
// List spPredict_ASMK2(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R, const int& J) {
//
//   // Subset number for prediction data
//   int u = X_u.n_rows;
//   int set_size_u = floor(u / J);
//   arma::vec sets_u = arma::regspace(0, set_size_u, set_size_u * (J - 1));
//   sets_u = join_vert(sets_u, arma::vec(1, arma::fill::ones) * u);
//
//   // initialize objects
//   arma::mat out_rZ;
//   arma::mat out_rY;
//
//   for (size_t l = 0; l < sets_u.n_elem - 1; ++l) {
//
//     // Define the current set
//     arma::uvec set_j = arma::regspace<arma::uvec>(sets_u(l), sets_u(l + 1) - 1);
//
//     // Extract the current subset from X
//     arma::mat sub_X_u = X_u.rows(set_j);
//     arma::mat sub_crd_u = crd_u.rows(set_j);
//
//     // Call the user-defined function FUN
//     List out_j = fast_BPSpred3(data, sub_X_u, priors, coords, sub_crd_u, hyperpar, W, R);
//
//     // Append the results to the output list
//     out_rZ = join_vert(out_rZ, as<arma::mat>(out_j["Z_hat"]));
//     out_rY = join_vert(out_rY, as<arma::mat>(out_j["Y_hat"]));
//   }
//
//   return List::create(Named("Z_hat") = out_rZ,
//                       Named("Y_hat") = out_rY);
//
// }
//
//
// // ##################################################################################################################################################
// // BACKUP FUNTCIONS - MULTIVARIATE ##################################################################################################################
// // ##################################################################################################################################################
//
//
// //' The same fit_latent_cpp, but take as argument the distance matrix directly (does not compute it by itself from coords)
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_B},\eqn{V_r},\eqn{\Psi},\eqn{\nu}
// //' @param dist [matrix] sample distance matrix
// //' @param hyperpar [list] two elemets: first named \eqn{\alpha}, second named \eqn{\phi}
// //'
// //' @return [list] posterior update parameters
// //'
// // [[Rcpp::export]]
// List fit_latent_cpp2(const List& data, const List& priors, const arma::mat& dist, const List& hyperpar) {
//
//   // Unpack data and priors
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   arma::mat mu_B = as<arma::mat>(priors["mu_B"]);
//   arma::mat V_r = as<arma::mat>(priors["V_r"]);
//   arma::mat Psi = as<arma::mat>(priors["Psi"]);
//   double nu = as<double>(priors["nu"]);
//   double alpha = as<double>(hyperpar["alpha"]);
//   double phi = as<double>(hyperpar["phi"]);
//
//   int n = Y.n_rows;
//   arma::mat Rphi_s = exp(-phi * dist);
//
//   // Precompute some reusable values
//   double a = (alpha / (1 - alpha));
//   arma::mat tX = trans(X);
//   arma::mat iV_r = arma::inv(V_r);
//   arma::mat iR_s = arma::inv(Rphi_s);
//
//   // Compute posterior updating
//   arma::mat V_B = a * tX * X + iV_r;
//   arma::mat V_BW = a * tX;
//   arma::mat V_WB = trans(V_BW);
//   arma::mat V_W = iR_s + (a * eye<arma::mat>(n, n));
//
//   arma::mat V_star1 = join_horiz(V_B, V_BW);
//   arma::mat V_star2 = join_horiz(V_WB, V_W);
//   arma::mat iV_star = join_vert( V_star1, V_star2);
//   arma::mat V_star = arma::inv(iV_star);
//
//   arma::mat M = join_vert( (a * tX * Y) + (iV_r * mu_B) , a * Y );
//   arma::mat mu_star = V_star * M;
//
//   arma::mat aYY = (a * trans(Y) * Y);
//   arma::mat mbVrmb = (trans(mu_B) * iV_r * mu_B);
//   arma::mat msVsms = (trans(mu_star) * iV_star * mu_star);
//   arma::mat Psi_star = Psi + aYY  + mbVrmb - msVsms;
//   double nu_star = nu + n;
//
//   // Return results as an R list
//   return List::create(Named("V_star") = V_star,
//                       Named("mu_star") = mu_star,
//                       Named("Psi_star") = Psi_star,
//                       Named("nu_star") = nu_star,
//                       Named("iRphi_s") = iR_s);
// }
//
//
// //' Draw from the conditional posterior predictive for a set of unobserved covariates
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param d_u [matrix] unobserved instances distance matrix
// //' @param d_us [matrix] cross-distance between unobserved and observed instances matrix
// //' @param hyperpar [list] two elemets: first named \eqn{\alpha}, second named \eqn{\phi}
// //' @param poster [list] output from \code{fit_cpp} function
// //' @param beta [matrix] posterior sample for \eqn{\beta}
// //' @param sigma [matrix] posterior sample for \eqn{\Sigma}
// //'
// //' @return [list] posterior predictive samples
// //'
// // [[Rcpp::export]]
// List r_pred_latent_MC(const List& data, const arma::mat& X_u, const arma::mat& d_u, const arma::mat& d_us, const List& hyperpar, const List& poster, const arma::mat& beta, const arma::mat& sigma) {
//
//   // Unpack data, posterior sample and hyperparameters
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   arma::mat iR_s = as<arma::mat>(poster["iRphi_s"]);
//   double alpha = as<double>(hyperpar["alpha"]);
//   double phi = as<double>(hyperpar["phi"]);
//
//   // extract info from data
//   int m = X_u.n_rows;
//   int n = d_us.n_rows-m;
//   int p = X_u.n_cols;
//
//   // covariance matrices
//   arma::mat Rphi_u = exp(-phi * d_u);
//   arma::mat Rphi_us = exp(-phi * d_us.submat(0, m, m-1, m+n-1));
//
//   // sampling environment
//   Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
//   Rcpp::Function rMNorm_R = mniw["rMNorm"];
//
//   // predictive conjugate parameters
//   arma::mat b = beta.rows(0, p - 1);
//   arma::mat w = beta.rows(p, (n+p) - 1);
//
//   // prediction W_u
//   arma::mat M_u = Rphi_us * iR_s;
//   arma::mat mu_u = M_u * w;
//   arma::mat V_u = Rphi_u - Rphi_us * iR_s * trans(Rphi_us);
//   arma::mat resultW = as<arma::mat>(rMNorm_R(Named("n", 1), Named("Lambda", mu_u), Named("SigmaR", V_u), Named("SigmaC", sigma)));
//
//   // prediction Y_u
//   arma::mat mu_y = X_u * b + resultW;
//   double a = alpha / (1 - alpha);
//   arma::mat V_y = a * eye<arma::mat>(m, m);
//   arma::mat resultY = as<arma::mat>(rMNorm_R(Named("n", 1), Named("Lambda", mu_y), Named("SigmaR", V_y), Named("SigmaC", sigma)));
//
//   return List::create(Named("Wu") = resultW,
//                       Named("Yu") = resultY);
//
// }
//
//
// //' Evaluate the density of a set of unobserved response with respect to the conditional posterior predictive
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param Y_u [matrix] unobserved instances response matrix
// //' @param d_u [matrix] unobserved instances distance matrix
// //' @param d_us [matrix] cross-distance between unobserved and observed instances matrix
// //' @param hyperpar [list] two elemets: first named \eqn{\alpha}, second named \eqn{\phi}
// //' @param poster [list] output from \code{fit_cpp} function
// //' @param beta [matrix] posterior sample for \eqn{\beta}
// //' @param sigma [matrix] posterior sample for \eqn{\Sigma}
// //'
// //' @return [double] posterior predictive density evaluation
// //'
// // [[Rcpp::export]]
// double d_pred_latent_MC(const List& data, const arma::mat& X_u, const arma::mat& Y_u, const arma::mat& d_u, const arma::mat& d_us, const List& hyperpar, const List& poster, const arma::mat& beta, const arma::mat& sigma) {
//
//   // Unpack data, posterior sample and hyperparameters
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   arma::mat iR_s = as<arma::mat>(poster["iRphi_s"]);
//   double alpha = as<double>(hyperpar["alpha"]);
//   double phi = as<double>(hyperpar["phi"]);
//
//   // extract info from data
//   int m = X_u.n_rows;
//   int n = d_us.n_rows-m;
//   int p = X_u.n_cols;
//
//   // covariance matrices
//   arma::mat Rphi_u = exp(-phi * d_u);
//   arma::mat Rphi_us = exp(-phi * d_us.submat(0, m, m-1, m+n-1));
//
//   // sampling environment
//   Rcpp::Environment mniw = Rcpp::Environment::namespace_env("mniw");
//   Rcpp::Function rMNorm_R = mniw["rMNorm"];
//   Rcpp::Function dMNorm_R = mniw["dMNorm"];
//
//   // predictive conjugate parameters
//   arma::mat b = beta.rows(0, p - 1);
//   arma::mat w = beta.rows(p, (n+p) - 1);
//
//   // prediction W_u
//   arma::mat M_u = Rphi_us * iR_s;
//   arma::mat mu_u = M_u * w;
//   arma::mat V_u = Rphi_u - Rphi_us * iR_s * trans(Rphi_us);
//   arma::mat resultW = as<arma::mat>(rMNorm_R(Named("n", 1), Named("Lambda", mu_u), Named("SigmaR", V_u), Named("SigmaC", sigma)));
//
//   // prediction Y_u
//   arma::mat mu_y = X_u * b + resultW;
//   double a = alpha / (1 - alpha);
//   arma::mat V_y = a * eye<arma::mat>(m, m);
//   double resultY = as<double>(dMNorm_R(Named("X", Y_u), Named("Lambda", mu_y), Named("SigmaR", V_y), Named("SigmaC", sigma)));
//
//   return resultY;
// }
//
//
// //' Compute the LOOCV of the density evaluations for fixed values of the hyperparameters
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_B},\eqn{V_r},\eqn{\Psi},\eqn{\nu}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\alpha}, second named \eqn{\phi}
// //'
// //' @return [vector] posterior predictive density evaluations
// //'
// // [[Rcpp::export]]
// arma::vec dens_loocv_latent_MC(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar) {
//
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   int n = Y.n_rows;
//   arma::vec predictions(n);
//
//   for (int i = 0; i < n; i++) {
//
//     // Rcout << "i: " << i << std::endl;
//
//     // extract the training data (excluding the i-th row) (create a mask to exclude the i-th row)
//     arma::mat Ytr, Xtr, coords_i;
//     if (i == 0) {
//       // Exclude the 0-th row
//       Ytr = Y.rows(1, n - 1);
//       Xtr = X.rows(1, n - 1);
//       coords_i = coords.rows(1, n - 1);
//     } else if (i == 1) {
//       // Exclude the 1-th row
//       Ytr = join_vert(Y.row(0), Y.rows(2, n - 1));
//       Xtr = join_vert(X.row(0), X.rows(2, n - 1));
//       coords_i = join_vert(coords.row(0), coords.rows(2, n - 1));
//     } else if (i == n - 1) {
//       // Exclude the last row
//       Ytr = Y.rows(0, n - 2);
//       Xtr = X.rows(0, n - 2);
//       coords_i = coords.rows(0, n - 2);
//     } else {
//       // Exclude the i-th row for i > 1
//       Ytr = join_vert(Y.rows(0, i - 1), Y.rows(i + 1, n - 1));
//       Xtr = join_vert(X.rows(0, i - 1), X.rows(i + 1, n - 1));
//       coords_i = join_vert(coords.rows(0, i - 1), coords.rows(i + 1, n - 1));
//     }
//     List data_i = List::create(
//       Named("Y") = Ytr,
//       Named("X") = Xtr);
//
//     // extract the test data
//     arma::mat crd_i = coords.row(i);
//     arma::mat X_i = X.row(i);
//     arma::mat Y_i = Y.row(i);
//
//     // Fit your model on the training data
//     List poster_i = fit_latent_cpp(data_i, priors, coords_i, hyperpar);
//
//     // evaluate predictive density
//     arma::mat d_i = arma_dist(crd_i);
//     arma::mat crd_is = join_cols(crd_i, coords_i);
//     arma::mat d_is = arma_dist(crd_is);
//
//     // posterior draws
//     List post = post_draws_latent(poster_i, 1);
//     List drw = as<List>(post(0));
//     arma::mat b = as<arma::mat>(drw["beta"]);
//     arma::mat s = as<arma::mat>(drw["sigma"]);
//
//     double dens = d_pred_latent_cpp(data_i, X_i, Y_i, d_i, d_is, hyperpar, poster_i, b, s);
//
//     // Store the prediction
//     predictions[i] = dens;
//   }
//
//   return predictions;
//
// }
//
//
// //' Compute the KCV of the density evaluations for fixed values of the hyperparameters
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_B},\eqn{V_r},\eqn{\Psi},\eqn{\nu}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\alpha}, second named \eqn{\phi}
// //' @param K [integer] number of folds
// //'
// //' @return [vector] posterior predictive density evaluations
// //'
// // [[Rcpp::export]]
// arma::vec dens_kcv_latent_MC(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, const int& K) { //, const int& seed) {
//
//   // unpack data
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   int n = Y.n_rows;
//   arma::vec predictions(n, arma::fill::zeros);
//   // Rcout << "n " << n << std::endl;
//
//   // Create a random permutation of indices from 1 to K
//   arma::vec p = arma::vec(K, arma::fill::ones) / K;
//   arma::uvec foldIndices = sample_index(K, n, p);
//
//   for (int k = 0; k < K; k++) {
//
//     // Define the indices for the current fold
//     arma::uvec testSet = find(foldIndices == (k + 1)); // Find indices that match the current fold number
//
//     // Create a boolean vector to identify the training set for this fold
//     arma::uvec trainSet = find(foldIndices != (k + 1)); // Find indices that do not match the current fold number
//
//     // Extract the training data for the current fold
//     arma::mat Ytr = Y.rows(trainSet);
//     arma::mat Xtr = X.rows(trainSet);
//     arma::mat coords_tr = coords.rows(trainSet);
//
//     List data_tr = List::create(
//       Named("Y") = Ytr,
//       Named("X") = Xtr);
//
//     // Extract the test data for the current fold
//     arma::mat Y_test = Y.rows(testSet);
//     arma::mat X_test = X.rows(testSet);
//     arma::mat coords_test = coords.rows(testSet);
//
//     // Fit your model on the training data
//     List poster_k = fit_latent_cpp(data_tr, priors, coords_tr, hyperpar);
//
//     // posterior draws
//     List post = post_draws_latent(poster_k, 1);
//     List drw = as<List>(post(0));
//     arma::mat b = as<arma::mat>(drw["beta"]);
//     arma::mat s = as<arma::mat>(drw["sigma"]);
//
//     // evaluate predictive density for the test set of the current fold
//     for (uword i = 0; i < testSet.n_elem; i++) {
//
//       arma::mat crd_i = coords_test.row(i);
//       arma::mat X_i = X_test.row(i);
//       arma::mat Y_i = Y_test.row(i);
//       arma::mat d_i = arma_dist(crd_i);
//       arma::mat crd_is = join_cols(crd_i, coords_tr);
//       arma::mat d_is = arma_dist(crd_is);
//
//       double dens = d_pred_latent_cpp(data_tr, X_i, Y_i, d_i, d_is, hyperpar, poster_k, b, s);
//       predictions(testSet(i)) = dens;
//     }
//   }
//
//   return predictions;
//
// }
//
//
// //' Return the CV predictive density evaluations for all the model combinations
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_B},\eqn{V_r},\eqn{\Psi},\eqn{\nu}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\alpha}, second named \eqn{\phi}
// //' @param useKCV [boolean] if \code{TRUE} K-fold cross validation is used instead of LOOCV (no \code{default})
// //' @param K [integer] number of folds
// //'
// //' @return [matrix] posterior predictive density evaluations (each columns represent a different model)
// //'
// // [[Rcpp::export]]
// arma::mat models_dens_latent_MC(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, bool useKCV, int K) {
//
//   // build the grid of hyperparameters
//   arma::vec Alfa = hyperpar["alpha"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Alfa, Fi);
//   int k = Grid.n_rows;
//
//   arma::mat out;
//
//   for(int j = 0; j < k; j++) {
//
//     // identify the model
//     arma::rowvec hpar = Grid.row(j);
//     double alfa = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("alpha") = alfa,
//       Named("phi") = fi);
//
//     // Call the appropriate function based on the 'useKCV' argument
//     arma::vec out_j;
//     if (useKCV) {
//       out_j = dens_kcv_latent(data, priors, coords, hmod, K);
//     } else {
//       out_j = dens_loocv_latent(data, priors, coords, hmod);
//     }
//
//     out =  join_horiz(out, out_j);
//
//   }
//
//   return out;
// }
//
//
// //' Compute the BPS spatial prediction given a set of stacking weights
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_B},\eqn{V_r},\eqn{\Psi},\eqn{\nu}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\alpha}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //'
// //' @return [list] BPS posterior predictive samples
// //'
// // [[Rcpp::export]]
// List BPS_latent_SpatialPrediction_cpp(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R) {
//
//   List pred(R);
//
//   // compute distance matrices
//   arma::mat d_u = arma_dist(crd_u);
//   arma::mat crd_us = join_cols(crd_u, coords);
//   arma::mat d_us = arma_dist(crd_us);
//
//   // build the grid of hyperparameters
//   arma::vec Alfa = hyperpar["alpha"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Alfa, Fi);
//   int k = Grid.n_rows;
//
//   // sample the model
//   arma::uvec kmod = sample_index(k, R, W);
//
//   for(int r = 0; r < R; r++) {
//
//     // identify the k-th model
//     arma::uword k_mod = kmod(r);
//     arma::rowvec hpar = Grid.row(k_mod);
//     double alfa = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("alpha") = alfa,
//       Named("phi") = fi);
//
//     // fit your model on the training data
//     List poster = fit_latent_cpp(data, priors, coords, hmod);
//
//     // posterior draws
//     List post = post_draws_latent(poster, 1);
//     List drw = post(0);
//     arma::mat b = as<arma::mat>(drw["beta"]);
//     arma::mat s = as<arma::mat>(drw["sigma"]);
//
//     // draw from conditional posterior predictive
//     pred(r) = r_pred_latent_cpp(data, X_u, d_u, d_us, hmod, poster, b, s);
//
//   }
//
//   return pred;
//
// }
//
//
// //' Compute the BPS weights by convex optimization
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param K [integer] number of folds
// //'
// //' @return [matrix] posterior predictive density evaluations (each columns represent a different model)
// //' @export
// // [[Rcpp::export]]
// List BPS_weights_MvT(const List& data, const List& priors, const arma::mat& coords, const List& hyperpar, int K) {
//
//   // compute predictive density evaluations
//   arma::mat out = models_dens_latent(data, priors, coords, hyperpar, true, K);
//
//   // compute the weights
//   arma::mat weights = as<arma::mat>(CVXR_opt(out));
//   weights.elem(find(weights <= 0)).zeros();
//   weights /= sum(weights, 0).eval()(0, 0);
//
//   // return the list (BPS weights, Grid, and predictive density evaluations)
//   arma::vec Alfa = hyperpar["alpha"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Alfa, Fi);
//   arma::mat res = join_horiz(weights, Grid);
//
//   List Res = List::create(
//     Named("Grid") = res,
//     Named("W") = weights,
//     Named("epd") = out
//   );
//
//   return Res;
// }
//
//
// //' Compute the BPS spatial prediction given a set of stacking weights
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //'
// //' @return [list] BPS posterior predictive samples
// //' @export
// // [[Rcpp::export]]
// List BPS_pred_MvT(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R) {
//
//   List pred(R);
//
//   // compute distance matrices
//   arma::mat d_u = arma_dist(crd_u);
//   arma::mat crd_us = join_cols(crd_u, coords);
//   arma::mat d_us = arma_dist(crd_us);
//
//   // build the grid of hyperparameters
//   arma::vec Alfa = hyperpar["alpha"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Alfa, Fi);
//   int k = Grid.n_rows;
//
//   // sample the model
//   arma::uvec kmod = sample_index(k, R, W);
//
//   for(int r = 0; r < R; r++) {
//
//     // identify the k-th model
//     arma::uword k_mod = kmod(r);
//     arma::rowvec hpar = Grid.row(k_mod);
//     double alfa = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("alpha") = alfa,
//       Named("phi") = fi);
//
//     // fit your model on the training data
//     List poster = fit_latent_cpp(data, priors, coords, hmod);
//
//     // posterior draws
//     List post = post_draws_latent(poster, 1);
//     List drw = post(0);
//     arma::mat b = as<arma::mat>(drw["beta"]);
//     arma::mat s = as<arma::mat>(drw["sigma"]);
//
//     // draw from conditional posterior predictive
//     pred(r) = r_pred_latent_cpp(data, X_u, d_u, d_us, hmod, poster, b, s);
//
//   }
//
//   return pred;
//
// }
//
//
// //' Compute the BPS spatial prediction given a set of stacking weights
// //'
// //' @param data [list] two elements: first named \eqn{Y}, second named \eqn{X}
// //' @param X_u [matrix] unobserved instances covariate matrix
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param coords [matrix] sample coordinates for X and Y
// //' @param crd_u [matrix] unboserved instances coordinates
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param W [matrix] set of stacking weights
// //' @param R [integer] number of desired samples
// //'
// //' @return [list] BPS posterior predictive samples
// //' @export
// // [[Rcpp::export]]
// List BPS_post_MvT(const List& data, const arma::mat& X_u, const List& priors, const arma::mat& coords, const arma::mat& crd_u, const List& hyperpar, const arma::vec& W, const int& R) {
//
//   List pred(R);
//   List post_smp(R);
//
//   // compute distance matrices
//   arma::mat d_u = arma_dist(crd_u);
//   arma::mat crd_us = join_cols(crd_u, coords);
//   arma::mat d_us = arma_dist(crd_us);
//
//   // build the grid of hyperparameters
//   arma::vec Alfa = hyperpar["alpha"];
//   arma::vec Fi = hyperpar["phi"];
//   arma::mat Grid = expand_grid_cpp(Alfa, Fi);
//   int k = Grid.n_rows;
//
//   // sample the model
//   arma::uvec kmod = sample_index(k, R, W);
//
//   for(int r = 0; r < R; r++) {
//
//     // identify the k-th model
//     arma::uword k_mod = kmod(r);
//     arma::rowvec hpar = Grid.row(k_mod);
//     double alfa = hpar[0];
//     double fi = hpar[1];
//     List hmod = List::create(
//       Named("alpha") = alfa,
//       Named("phi") = fi);
//
//     // fit your model on the training data
//     List poster = fit_latent_cpp(data, priors, coords, hmod);
//
//     // posterior draws
//     List post = post_draws_latent(poster, 1);
//     List drw = post(0);
//     arma::mat b = as<arma::mat>(drw["beta"]);
//     arma::mat s = as<arma::mat>(drw["sigma"]);
//
//     // draw from conditional posterior predictive
//     pred(r) = r_pred_latent_cpp(data, X_u, d_u, d_us, hmod, poster, b, s);
//
//     // save posterior draw
//     post_smp(r) = drw;
//   }
//
//   return List::create(Named("Pred") = pred,
//                       Named("Post") = post_smp);
//
// }
