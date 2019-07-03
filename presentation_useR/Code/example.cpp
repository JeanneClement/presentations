#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]

Rcpp::List alpha_sample(arma::mat X, arma::mat beta,
                        arma::mat Y, int Valpha, const int seed) {
  // Defining constants 
  const int NSITE = Y.n_rows;
  const int NSPECIES = Y.n_cols; 
  
  // Initialize random number generator 
  gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(s, seed);
  
  // Declaring new objects to store results 
  arma::vec alpha; alpha.zeros(NSITE);
  
  // Draw in the posterior distribution
  for (int i=0; i < NSITE ; i++) {
    double small_v = arma::sum(Y.row(i)-X.row(i)*beta);
    double big_V = 1/(1/Valpha + NSPECIES);
    alpha(i) = big_V*small_v + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V));
  }
  
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("alpha") = alpha);  
  return results;
}

