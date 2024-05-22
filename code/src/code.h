
#ifndef CODE_H
#define CODE_H

// Function declarations (conquer-step)

SEXP CVXR_opt(const arma::mat& scores);

Rcpp::List BPS_combine(const Rcpp::List& fit_list, const int& K, const double& rp = 1);

Rcpp::List BPS_PseudoBMA(const Rcpp::List& fit_list);

// Function declarations (univariate models)

Rcpp::List fit_cpp(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar);

Rcpp::List post_draws(const Rcpp::List& poster, const int& R = 50, const bool& par = false, const int& p = 1);

Rcpp::List r_pred_joint(const Rcpp::List& data, const arma::mat& X_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster, const int& R = 1);

Rcpp::List r_pred_marg(const Rcpp::List& data, const arma::mat& X_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster, const int& R = 1);

Rcpp::List r_pred_cond(const Rcpp::List& data, const arma::mat& X_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster, const Rcpp::List& post);

double d_pred_cpp(const Rcpp::List& data, const arma::mat& X_u, const arma::vec& Y_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster);

arma::vec dens_loocv(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar);

arma::vec dens_kcv(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, const int& K = 5);

arma::mat models_dens(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, const bool& useKCV = true, const int& K = 5);

Rcpp::List BPS_weights(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, int K = 5);

Rcpp::List BPS_pred(const Rcpp::List& data, const arma::mat& X_u, const Rcpp::List& priors, const arma::mat& coords, const arma::mat& crd_u, const Rcpp::List& hyperpar, const arma::vec& W, const int& R = 1);

arma::mat BPS_postdraws(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, const arma::vec& W, const int& R = 1);

Rcpp::List BPS_post(const Rcpp::List& data, const arma::mat& X_u, const Rcpp::List& priors, const arma::mat& coords, const arma::mat& crd_u, const Rcpp::List& hyperpar, const arma::vec& W, const int& R);

Rcpp::List spPredict_BPS(const Rcpp::List& data, const arma::mat& X_u, const Rcpp::List& priors, const arma::mat& coords, const arma::mat& crd_u, const Rcpp::List& hyperpar, const arma::vec& W, const int& R = 1, const int& J = 1);

// Function declarations (multivariate models)

Rcpp::List fit_cpp_MvT(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar);

Rcpp::List post_draws_MvT(const Rcpp::List& poster, const int& R = 1, const bool& par = false, const int& p = 2);

Rcpp::List r_pred_joint_MvT(const Rcpp::List& data, const arma::mat& X_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster, const int& R = 1);

Rcpp::List r_pred_marg_MvT(const Rcpp::List& data, const arma::mat& X_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster, const int& R = 1);

Rcpp::List r_pred_cond_MvT(const Rcpp::List& data, const arma::mat& X_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster, const Rcpp::List& post);

double d_pred_cpp_MvT(const Rcpp::List& data, const arma::mat& X_u, const arma::mat& Y_u, const arma::mat& d_u, const arma::mat& d_us, const Rcpp::List& hyperpar, const Rcpp::List& poster);

arma::vec dens_loocv_MvT(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar);

arma::vec dens_kcv_MvT(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, const int& K);

arma::mat models_dens_MvT(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, bool useKCV, int K = 10);

Rcpp::List BPS_weights_MvT(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, int K);

Rcpp::List BPS_pred_MvT(const Rcpp::List& data, const arma::mat& X_u, const Rcpp::List& priors, const arma::mat& coords, const arma::mat& crd_u, const Rcpp::List& hyperpar, const arma::vec& W, const int& R = 1);

Rcpp::List BPS_post_MvT(const Rcpp::List& data, const arma::mat& X_u, const Rcpp::List& priors, const arma::mat& coords, const arma::mat& crd_u, const Rcpp::List& hyperpar, const arma::vec& W, const int& R);

Rcpp::List BPS_postdraws_MvT(const Rcpp::List& data, const Rcpp::List& priors, const arma::mat& coords, const Rcpp::List& hyperpar, const arma::vec& W, const int& R = 1,  bool par = false);

#endif
