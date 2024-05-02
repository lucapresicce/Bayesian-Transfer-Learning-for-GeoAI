
#ifndef UTILSC_H
#define UTILSC_H

// Function declarations

arma::mat arma_dist(const arma::mat & X);

arma::mat expand_grid_cpp(const arma::vec& x, const arma::vec& y);

arma::uvec sample_index(const int& size, const int& length, const arma::vec& p);

Rcpp::List subset_data(const Rcpp::List& data, int K);

arma::mat forceSymmetry_cpp(const arma::mat& mat);

#endif // UTILSC_H
