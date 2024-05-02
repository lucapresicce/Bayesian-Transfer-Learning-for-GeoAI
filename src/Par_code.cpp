// #include <RcppArmadillo.h>
// #include "code.h"
// #include "utilsC.h"
//
// // #include <omp.h>
// // // [[Rcpp::plugins(openmp)]]
//
// #include <RcppParallel.h>
// #include <RcppThread.h>
// #include <RcppClock.h>
// #include <thread>
// //[[Rcpp::depends(RcppClock)]]
//
// using namespace Rcpp;
// using namespace arma;
//
//
// C++ general OpenMP code structure
//
// main ()  {
//
//   int var1, var2, var3;
//
//   Serial code
//     .
//     .
//     .
//
//   Beginning of parallel section. Fork a team of threads.
//   Specify variable scoping
//
// #pragma omp parallel private(var1, var2) shared(var3)
//
// {
//
//   Parallel section executed by all threads
//   .
//   Other OpenMP directives
//   .
//   Run-time Library calls
//   .
//   All threads join master thread and disband
//
// }
//
// Resume serial code
//   .
//   .
//   .
//
// }
//
//
// //' Perform Accelerated Spatial Meta Kriging (ASMK)
// //'
// //' @param data [list] three elements: first named \eqn{Y}, second named \eqn{X}, third named \eqn{crd}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param K [integer] number of subsets
// //' @param newdata [list] two elements: second named \eqn{X}, third named \eqn{crd}
// //' @param R [integer] number of poterior predictive sample
// //'
// //' @return [list] posterior update parameters
// //'
// // [[Rcpp::export]]
// List spASMK(const List& data, const List& priors, const List& hyperpar, const int& K, const List& newdata, const int& R = 250) {
//
//   // Profiling
//   Rcpp::Clock clock; // mind the headers: //[[Rcpp::depends(RcppClock)]] #include <RcppClock.h> #include <thread>
//
//   // Unpack data
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   arma::mat coords = as<arma::mat>(data["crd"]);
//   // int n = Y.n_rows;
//   // int q = Y.n_cols;
//   // int p = X.n_cols;
//
//   // // hyperparamters grid
//   // arma::vec Delta = hyperpar["delta"];
//   // arma::vec Phi = hyperpar["phi"];
//   // arma::mat Grid = expand_grid_cpp(Delta, Phi);
//   // int K = Grid.n_rows;
//
//   // Subset data
//   List subsets = subset_data(data, K);
//   List Y_list = subsets[0];
//   List X_list = subsets[1];
//   List crd_list = subsets[2];
//
//   // Fit subset BPS-GP
//   List fit_list(K);
//   double threshold = 1.0 / (2.0 * K);
//
//   clock.tick("Subset model fitting");
//   // for loop over K
//   for (int i = 0; i < K; ++i) {
//
//     // subset data
//     arma::mat Y_i = as<arma::mat>(Y_list[i]);
//     arma::mat X_i = as<arma::mat>(X_list[i]);
//     arma::mat crd_i = as<arma::mat>(crd_list[i]);
//     List data_i = List::create(Named("Y") = Y_i,
//                                Named("X") = X_i);
//
//     // fit subset model
//     List out_i = BPS_weights(data_i, priors, crd_i, hyperpar);
//
//     // extract results
//     arma::mat epd_i = as<arma::mat>(out_i["epd"]);
//     arma::mat W_i = as<arma::mat>(out_i["W"]);
//     W_i.elem(find(W_i < threshold)).zeros();
//
//     // return
//     List fit_i = List::create(Named("epd") = epd_i,
//                               Named("W") = W_i);
//     fit_list[i] = fit_i;
//
//   }
//   clock.tock("Subset model fitting");
//
//   // Combine subset models
//   clock.tick("Subset models aggregation");
//   List comb_list = BPS_combine(fit_list, K);
//   clock.tock("Subset models aggregation");
//   arma::mat Wbps = as<arma::mat>(comb_list[0]);
//   List W_list = comb_list[1];
//
//   // Perform predictions
//   // List pred_list(R); // think about two matrices instead
//   arma::mat pred_Z;
//   arma::mat pred_Y;
//   arma::vec W_vec = arma::conv_to<arma::vec>::from(Wbps.col(0));
//   arma::uvec subset_ind = sample_index(K, R, W_vec);
//
//   clock.tick("Predictions");
//   // for loop over R
//   for (int r = 0; r < R; ++r) {
//
//     // model sample
//     int ind_r = subset_ind(r);
//     arma::mat Y_r = as<arma::mat>(Y_list[ind_r]);
//     arma::mat X_r = as<arma::mat>(X_list[ind_r]);
//     arma::mat crd_r = as<arma::mat>(crd_list[ind_r]);
//     arma::mat W_r = as<arma::mat>(W_list[ind_r]);
//     List data_r = List::create(Named("Y") = Y_r,
//                                Named("X") = X_r);
//
//     // newdata
//     arma::mat X_u = as<arma::mat>(newdata["X"]);
//     arma::mat crd_u = as<arma::mat>(newdata["crd"]);
//
//     // perform predictions
//     // int jj = X_u.n_rows/500;
//     List out_r = spPredict_ASMK(data_r, X_u, priors, crd_r, crd_u, hyperpar, W_r, 1, 1);
//     // pred_list[r] = out_r;
//     arma::mat out_rZ = out_r[0];
//     pred_Z = join_horiz(pred_Z, out_rZ);
//     arma::mat out_rY = out_r[1];
//     pred_Y = join_horiz(pred_Y, out_rY);
//
//   }
//   clock.tock("Predictions");
//
//   clock.stop("Timing");
//
//   // Return results
//   List pred_list = List::create(Named("Z_hat") = pred_Z,
//                                 Named("Y_hat") = pred_Y);
//
//   return List::create(Named("Predictions") = pred_list,
//                       Named("Comb_weights") = Wbps);
//
// }
//
//
// //  perform accelerated spatial meta kriging (ASMK)
// // [[Rcpp::export]]
// List spASMK(const List& data, const List& priors, const List& hyperpar, const int& K, const List& newdata, const int& R = 250, const int& num_cores = 1) {
//
//   // Profiling
//   Rcpp::Clock clock;
//
//   // Set the number of cores for OpenMP
//   omp_set_num_threads(num_cores);
//
//   // Unpack data
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   arma::mat coords = as<arma::mat>(data["crd"]);
//
//   // Subset data
//   List subsets = subsetAndSample(data, K);
//   List Y_list = subsets[0];
//   List X_list = subsets[1];
//   List crd_list = subsets[2];
//
//   // Fit subset BPS-GP
//   List fit_list(K);
//   double threshold = 1.0 / (2.0 * K);
//
//   clock.tick("Subset model fitting");
//   // Parallelize the loop over K using OpenMP
// #pragma omp parallel for
//   for (int i = 0; i < K; ++i) {
//
//     // thread-local storage for subset-specific variables
//     arma::mat Y_i_local = as<arma::mat>(Y_list[i]);
//     arma::mat X_i_local = as<arma::mat>(X_list[i]);
//     arma::mat crd_i_local = as<arma::mat>(crd_list[i]);
//     List data_i_local = List::create(Named("Y") = Y_i_local, Named("X") = X_i_local);
//
//     // fit subset model
//     List out_i_local = BPSweights_cpp2(data_i_local, priors, crd_i_local, hyperpar);
//
//     // extract results
//     arma::mat epd_i = as<arma::mat>(out_i_local["epd"]);
//     arma::mat W_i = as<arma::mat>(out_i_local["W"]);
//     W_i.elem(find(W_i < threshold)).zeros();
//
//     // return
//     List fit_i = List::create(Named("epd") = epd_i, Named("W") = W_i);
//
//     // Use a critical section to safely update shared fit_list
// #pragma omp critical
//     fit_list[i] = fit_i;
//   }
//   clock.tock("Subset model fitting");
//
//   clock.tick("Subset models aggregation");
//   // Combine subset models
//   List comb_list = BPS_combine(fit_list, K);
//   clock.tock("Subset models aggregation");
//   arma::mat Wbps = as<arma::mat>(comb_list[0]);
//   List W_list = comb_list[1];
//
//   // Perform predictions
//   arma::mat pred_Z;
//   arma::mat pred_Y;
//   arma::vec W_vec = arma::conv_to<arma::vec>::from(Wbps.col(0));
//   arma::uvec subset_ind = sample_index(K, R, W_vec);
//
//   clock.tick("Predictions");
//   // Parallelize the loop over R using OpenMP
// #pragma omp parallel for
//   for (int r = 0; r < R; ++r) {
//
//     // thread-local storage for prediction-specific variables
//     arma::mat pred_Z_local;
//     arma::mat pred_Y_local;
//
//     // model sample
//     int ind_r = subset_ind(r);
//     arma::mat Y_r_local = as<arma::mat>(Y_list[ind_r]);
//     arma::mat X_r_local = as<arma::mat>(X_list[ind_r]);
//     arma::mat crd_r_local = as<arma::mat>(crd_list[ind_r]);
//     arma::mat W_r_local = as<arma::mat>(W_list[ind_r]);
//     List data_r_local = List::create(Named("Y") = Y_r_local, Named("X") = X_r_local);
//
//     // newdata
//     arma::mat X_u = as<arma::mat>(newdata["X"]);
//     arma::mat crd_u = as<arma::mat>(newdata["crd"]);
//
//     // perform predictions
//     List out_r_local = spPredict_ASMK(data_r_local, X_u, priors, crd_r_local, crd_u, hyperpar, W_r_local, 1, 1);
//     arma::mat out_rZ_local = out_r_local[0];
//     arma::mat out_rY_local = out_r_local[1];
//
//     // Use a critical section to safely update shared pred_Z and pred_Y
// #pragma omp critical
// {
//   pred_Z_local = out_rZ_local;
//   pred_Y_local = out_rY_local;
//   pred_Z = join_horiz(pred_Z, pred_Z_local);
//   pred_Y = join_horiz(pred_Y, pred_Y_local);
// }
//   }
//
//   clock.tock("Predictions");
//
//   clock.stop("Timing");
//
//   // Return results
//   List pred_list = List::create(Named("Z_hat") = pred_Z, Named("Y_hat") = pred_Y);
//   return List::create(Named("Predictions") = pred_list, Named("Comb_weights") = Wbps, Named("subsetind") = subset_ind);
// }

// //' Perform Accelerated Spatial Meta Kriging (ASMK) PARALLEL
// //'
// //' @param data [list] three elements: first named \eqn{Y}, second named \eqn{X}, third named \eqn{crd}
// //' @param priors [list] priors: named \eqn{\mu_b},\eqn{V_b},\eqn{a},\eqn{b}
// //' @param hyperpar [list] two elemets: first named \eqn{\delta}, second named \eqn{\phi}
// //' @param K [integer] number of subsets
// //' @param newdata [list] two elements: second named \eqn{X}, third named \eqn{crd}
// //' @param R [integer] number of poterior predictive sample
// //'
// //' @return [list] posterior update parameters
// //' @export
// // [[Rcpp::export]]
// List PARspASMK(const List& data, const List& priors, const List& hyperpar, const int& K, const List& newdata, const int& R = 250, const int& num_cores = 1) {
//
//   // Profiling
//   Rcpp::Clock clock;
//
//   // Set the number of cores for RcppThread
//   RcppThread::ThreadPool pool(num_cores);
//
//   // Unpack data
//   arma::mat Y = as<arma::mat>(data["Y"]);
//   arma::mat X = as<arma::mat>(data["X"]);
//   arma::mat coords = as<arma::mat>(data["crd"]);
//
//   // Subset data
//   List subsets = subset_data(data, K);
//   List Y_list = subsets[0];
//   List X_list = subsets[1];
//   List crd_list = subsets[2];
//
//   // Fit subset BPS-GP
//   List fit_list(K);
//   double threshold = 1.0 / (2.0 * K);
//
//   clock.tick("Subset model fitting");
//
//   // Parallelize the loop over K using RcppThread
//   RcppThread::parallelFor(0, K, [&](int i) {
//     // thread-local storage for subset-specific variables
//     arma::mat Y_i_local = as<arma::mat>(Y_list[i]);
//     arma::mat X_i_local = as<arma::mat>(X_list[i]);
//     arma::mat crd_i_local = as<arma::mat>(crd_list[i]);
//     List data_i_local = List::create(Named("Y") = Y_i_local, Named("X") = X_i_local);
//
//     // fit subset model
//     List out_i_local = BPSweights_cpp2(data_i_local, priors, crd_i_local, hyperpar);
//
//     // extract results
//     arma::mat epd_i = as<arma::mat>(out_i_local["epd"]);
//     arma::mat W_i = as<arma::mat>(out_i_local["W"]);
//     W_i.elem(find(W_i < threshold)).zeros();
//
//     // return
//     List fit_i = List::create(Named("epd") = epd_i, Named("W") = W_i);
//
//     // Use a critical section to safely update shared fit_list
//     // RcppThread::Rcout lock;
//     fit_list[i] = fit_i;
//   });
//
//   clock.tock("Subset model fitting");
//
//   clock.tick("Subset models aggregation");
//   // Combine subset models
//   List comb_list = BPS_combine(fit_list, K);
//   clock.tock("Subset models aggregation");
//   arma::mat Wbps = as<arma::mat>(comb_list[0]);
//   List W_list = comb_list[1];
//
//   // Perform predictions
//   arma::mat pred_Z;
//   arma::mat pred_Y;
//   arma::vec W_vec = arma::conv_to<arma::vec>::from(Wbps.col(0));
//   arma::uvec subset_ind = sample_index(K, R, W_vec);
//
//   clock.tick("Predictions");
//
//   // Parallelize the loop over R using RcppThread
//   RcppThread::parallelFor(0, R, [&](int r) {
//     // thread-local storage for prediction-specific variables
//     arma::mat pred_Z_local;
//     arma::mat pred_Y_local;
//
//     // model sample
//     int ind_r = subset_ind(r);
//     arma::mat Y_r_local = as<arma::mat>(Y_list[ind_r]);
//     arma::mat X_r_local = as<arma::mat>(X_list[ind_r]);
//     arma::mat crd_r_local = as<arma::mat>(crd_list[ind_r]);
//     arma::mat W_r_local = as<arma::mat>(W_list[ind_r]);
//     List data_r_local = List::create(Named("Y") = Y_r_local, Named("X") = X_r_local);
//
//     // newdata
//     arma::mat X_u = as<arma::mat>(newdata["X"]);
//     arma::mat crd_u = as<arma::mat>(newdata["crd"]);
//
//     // perform predictions
//     List out_r_local = spPredict_ASMK(data_r_local, X_u, priors, crd_r_local, crd_u, hyperpar, W_r_local, 1, 1);
//     arma::mat out_rZ_local = out_r_local[0];
//     arma::mat out_rY_local = out_r_local[1];
//
//     // Use a critical section to safely update shared pred_Z and pred_Y
//     // RcppThread::Rcout lock;
//     pred_Z_local = out_rZ_local;
//     pred_Y_local = out_rY_local;
//     pred_Z = join_horiz(pred_Z, pred_Z_local);
//     pred_Y = join_horiz(pred_Y, pred_Y_local);
//   });
//
//   clock.tock("Predictions");
//
//   clock.stop("Timing");
//
//   // Return results
//   List pred_list = List::create(Named("Z_hat") = pred_Z, Named("Y_hat") = pred_Y);
//   return List::create(Named("Predictions") = pred_list, Named("Comb_weights") = Wbps, Named("subsetind") = subset_ind);
// }
