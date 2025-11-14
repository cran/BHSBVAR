#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Non-Central T-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_nonc_t(const double a1, const double p1, const double sigma1, const double nu, const double lam1) {
  double den = 0.0;
  if (lam1 < 0.0) {
    den = R::dnt((((-1.0) * (a1 - p1)) / sigma1), nu, std::abs(lam1), 0) / sigma1;
  } else {
    den = R::dnt(((a1 - p1) / sigma1), nu, lam1, 0) / sigma1;
  }
  return den;
}

// T-Distribution Truncated to be Negative.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_t_n(const double a1, const double p1, const double sigma1, const double nu) {
  double den = R::dt(((a1 - p1) / sigma1), nu, 0) / (sigma1 * (R::pt((((-1) * p1) / sigma1), nu, 1, 0)));
  return den;
}

// T-Distribution Truncated to be Positive.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_t_p(const double a1, const double p1, const double sigma1, const double nu) {
  double den = R::dt(((a1 - p1) / sigma1), nu, 0) / (sigma1 * (1 - R::pt((((-1) * p1) / sigma1), nu, 1, 0)));
  return den;
}

// T-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_t(const double a1, const double p1, const double sigma1, const double nu) {
  double den = R::dt(((a1 - p1) / sigma1), nu, 0) / sigma1;
  return den;
}

// Beta-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_beta(const double a1, const double sh1, const double sh2) {
  double den = R::dbeta(a1, sh1, sh2, 0);
  return den;
}

// Inverted Beta-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_ibeta(const double a1, const double sh1, const double sh2) {
  // returns NaN if a1 is less than 1.0
  double den = 0.0;
  if (a1 >= 1.0) {
    den = std::exp(((sh2 - 1.0) * std::log((a1 - 1.0))) + (((-1.0) * (sh2 + sh1)) * std::log((1.0 + (a1 - 1.0)))) - std::log(R::beta(sh2, sh1)));
  }
  return den;
}

// Sum the Log of Prior Densities.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double sum_log_prior_densities(const arma::mat& A_test, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH) {
  arma::uword nrow = pA.n_rows, ncol = pA.n_cols;
  arma::mat H_test = arma::inv(A_test);
  double sum_log_priors = 0.0, detA_test = arma::det(A_test);
  //compute sum of log priors
  for (arma::uword i = 0; i < nrow; ++i) {
    for (arma::uword j = 0; j < ncol; ++j) {
      if (std::isfinite(pA(i, j, 0))) {
        // if pA(i, j, 0) == 0, symmetric t-distribution
        if (pA(i, j, 0) == 0) {
          // if pA(i, j, 1) == 1,  positive sign restrictions
          if (pA(i, j, 1) == 1) {
            sum_log_priors += std::log(prior_t_p(A_test(i, j), pA(i, j, 2), pA(i, j, 3), pA(i, j, 4)));
          } else if (pA(i, j, 1) == (-1)) { // if pA(i, j, 1) == -1,  negative sign restrictions
            sum_log_priors += std::log(prior_t_n(A_test(i, j), pA(i, j, 2), pA(i, j, 3), pA(i, j, 4)));
          } else { // if pA(i, j, 1) == NA,  no sign restrictions
            sum_log_priors += std::log(prior_t(A_test(i, j), pA(i, j, 2), pA(i, j, 3), pA(i, j, 4)));
          }
        }
        // if pA(i, j, 0) == 1, non-central t-distribution
        if (pA(i, j, 0) == 1) {
          sum_log_priors += std::log(prior_nonc_t(A_test(i, j), pA(i, j, 2), pA(i, j, 3), pA(i, j, 4), pA(i, j, 5)));
        }
        // if pA(i, j, 0) == 2, inverted beta-distribution
        if (pA(i, j, 0) == 2) {
          sum_log_priors += std::log(prior_ibeta((pA(i, j, 1) * A_test(i, j)), pA(i, j, 3), pA(i, j, 4)));
        }
        // if pA(i, j, 0) == 3, beta-distribution
        if (pA(i, j, 0) == 3) {
          sum_log_priors += std::log(prior_beta((pA(i, j, 1) * A_test(i, j)), pA(i, j, 3), pA(i, j, 4)));
        }        
      }
      if (std::isfinite(pH(i, j, 0))) {
        // if pH(i, j, 0) == 0, symmetric t-distribution
        if (pH(i, j, 0) == 0) {
          // if pH(i, j, 1) == 1,  positive sign restrictions
          if (pH(i, j, 1) == 1) {
            sum_log_priors += std::log(prior_t_p(H_test(i, j), pH(i, j, 2), pH(i, j, 3), pH(i, j, 4)));
          } else if (pH(i, j, 1) == (-1)) { // if pH(i, j, 1) == -1,  negative sign restrictions
            sum_log_priors += std::log(prior_t_n(H_test(i, j), pH(i, j, 2), pH(i, j, 3), pH(i, j, 4)));
          } else { // if pH(i, j, 1) == NA,  no sign restrictions
            sum_log_priors += std::log(prior_t(H_test(i, j), pH(i, j, 2), pH(i, j, 3), pH(i, j, 4)));
          }
        }
        // if pH(i, j, 0) == 1, non-central t-distribution
        if (pH(i, j, 0) == 1) {
          sum_log_priors += std::log(prior_nonc_t(H_test(i, j), pH(i, j, 2), pH(i, j, 3), pH(i, j, 4), pH(i, j, 5)));
        }
        // if pH(i, j, 0) == 2, inverted beta-distribution
        if (pH(i, j, 0) == 2) {
          sum_log_priors += std::log(prior_ibeta((pH(i, j, 1) * H_test(i, j)), pH(i, j, 3), pH(i, j, 4)));
        }
        // if pH(i, j, 0) == 3, beta-distribution
        if (pH(i, j, 0) == 3) {
          sum_log_priors += std::log(prior_beta((pH(i, j, 1) * H_test(i, j)), pH(i, j, 3), pH(i, j, 4)));
        }
      }
    }
  }
  if (std::isfinite(pdetA(0, 0, 0))) {
    // if pdetA(0, 0, 0) == 0, symmetric t-distribution
    if (pdetA(0, 0, 0) == 0) {
      // if pdetA(0, 0, 0) == 1,  positive sign restrictions
      if (pdetA(0, 0, 1) == 1) {
        sum_log_priors += std::log(prior_t_p(detA_test, pdetA(0, 0, 2), pdetA(0, 0, 3), pdetA(0, 0, 4)));
      } else if (pdetA(0, 0, 1) == (-1)) { // if pdetA(0, 0, 1) == -1,  negative sign restrictions
        sum_log_priors += std::log(prior_t_n(detA_test, pdetA(0, 0, 2), pdetA(0, 0, 3), pdetA(0, 0, 4)));
      } else { // if pdetA(0, 0, 1) == NA, no sign restrictions
        sum_log_priors += std::log(prior_t(detA_test,pdetA(0, 0, 2), pdetA(0 ,0 ,3), pdetA(0, 0, 4)));
      }
    }
    // if pdetA(0, 0, 0) == 1, non-central t-distribution
    if (pdetA(0, 0, 0) == 1) {
      sum_log_priors += std::log(prior_nonc_t(detA_test, pdetA(0, 0, 2), pdetA(0, 0, 3), pdetA(0, 0, 4), pdetA(0, 0, 5)));
    }
    // if pdetA(0, 0, 0) == 2, inverted beta-distribution
    if (pdetA(0, 0, 0) == 2) {
      sum_log_priors += std::log(prior_ibeta((pdetA(0, 0, 1) * detA_test), pdetA(0, 0, 3), pdetA(0, 0, 4)));
    }
    // if pdetA(0, 0, 0) == 3, beta-distribution
    if (pdetA(0, 0, 0) == 3) {
      sum_log_priors += std::log(prior_beta((pdetA(0, 0, 1) * detA_test), pdetA(0, 0, 3), pdetA(0, 0, 4)));
    }
  }
  
  return sum_log_priors;
}

// Likelihood Function.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double log_likelihood_function(const arma::mat& A_test, const arma::mat& kappa1, const arma::mat& y1, const arma::mat& omega, const arma::mat& zeta_test, const arma::mat& somega) {
  arma::uword kappa_ncol = kappa1.n_cols;
  double ynrow = double (y1.n_rows);
  double lik_A_numerator = (ynrow / 2.0) * std::log(arma::det(A_test.t() * omega * A_test));
  arma::mat tau = arma::diagmat(kappa1) * arma::diagmat(A_test.t() * somega * A_test);
  for (arma::uword i = 0; i < kappa_ncol; ++i) {
    if (kappa1(0, i) > 0.0) {
      lik_A_numerator += kappa1(0, i) * std::log(tau(i, i));
    }
  }
  double lik_A_denomenator = arma::as_scalar((kappa1 + ((ynrow / 2.0) * arma::ones(1, kappa_ncol))) * arma::log(arma::diagvec((2.0 / ynrow) * (tau + (arma::diagmat(zeta_test) / 2.0)))));
  double lik_A = lik_A_numerator - lik_A_denomenator;
  return lik_A;
}

// Posterior Density Function.
double post_A_function(const arma::mat& A_test, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& kappa1, const arma::mat& y1, const arma::mat& omega, const arma::mat& zeta_test, const arma::mat& somega) {
  double priors = sum_log_prior_densities(A_test, pA, pdetA, pH);
  double likelihood = log_likelihood_function(A_test, kappa1, y1, omega, zeta_test, somega);
  double posterior = priors + likelihood;
  return posterior;
}

// Generate Proposal Values for A.
arma::mat proposal_function(const arma::mat& A_old, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& scale1) {
  
  const arma::uword nrow = pA.n_rows, ncol = pA.n_cols;
  
  arma::uword nH = 0;
  for (arma::uword i = 0; i < ncol; ++i) {
    for (arma::uword j = 0; j < nrow; ++j) {
      if ((std::isfinite(pH(j, i, 0))) && (pH(j, i, 0) == 0.0) && (std::isfinite(pH(j, i, 1)))) {
        nH += 1;
      }
    }
  }
  
  arma::mat A_test(nrow, ncol), H_test(nrow, ncol);
  std::fill(A_test.begin(), A_test.end(), Rcpp::NumericVector::get_na());
  std::fill(H_test.begin(), H_test.end(), Rcpp::NumericVector::get_na());
  
  double detA_test = 0.0;
  
  arma::uword sign_test = 1, aa = 0, bb = 0, cc = 0;
  
  while (sign_test != 0) {
    
    sign_test = 0;
    aa = 0;
    bb += 1;
    
    if (bb == 1000) {
      std::string message = "No draws match all sign restrictions:\n  Iterations: " + std::to_string(bb);
      Rcpp::stop(message);
    }
    
    for (arma::uword i = 0; i < ncol; ++i) {
      for (arma::uword j = 0; j < nrow; ++j) {
        if (std::isfinite(pA(j, i, 0))) {
          
          A_test(j, i) = A_old(j, i) + scale1(aa, aa) * (R::rnorm(0, 1) / std::sqrt(0.5 * (std::pow(R::rnorm(0, 1), 2) + std::pow(R::rnorm(0, 1), 2))));
          
          if (std::isfinite(pA(j, i, 1))) {
            if (pA(j, i, 0) == 0.0) {
              while (((pA(j, i, 1) > 0.0) && (A_test(j, i) < 0.0)) || ((pA(j, i, 1) < 0.0) && (A_test(j, i) > 0.0))) {
                A_test(j, i) = A_old(j, i) + scale1(aa, aa) * (R::rnorm(0, 1) / std::sqrt(0.5 * (std::pow(R::rnorm(0, 1), 2) + std::pow(R::rnorm(0, 1), 2))));
                cc += 1;
                if (cc == 1000) {
                  std::string message = "No draws match all sign restrictions:\n  Iterations: " + std::to_string(cc);
                  Rcpp::stop(message);
                }
              }
            }
            if (pA(j, i, 0) == 2.0){
              while (((pA(j, i, 1) == 1.0) && (A_test(j, i) < 1.0)) || ((pA(j, i, 1) == -1.0) && (A_test(j, i) > -1.0))) {
                A_test(j, i) = A_old(j, i) + scale1(aa, aa) * (R::rnorm(0, 1) / std::sqrt(0.5 * (std::pow(R::rnorm(0, 1), 2) + std::pow(R::rnorm(0, 1), 2))));
                cc += 1;
                if (cc == 1000) {
                  std::string message = "No draws match all sign restrictions:\n  Iterations: " + std::to_string(cc);
                  Rcpp::stop(message);
                }
              }
            }
            if (pA(j, i, 0) == 3.0){
              while (((pA(j, i, 1) == 1.0) && ((A_test(j, i) < 0.0) || (A_test(j, i) > 1.0))) || ((pA(j, i, 1) == -1.0) && ((A_test(j, i) > 0.0) || (A_test(j, i) < -1.0)))) {
                A_test(j, i) = A_old(j, i) + scale1(aa, aa) * (R::rnorm(0, 1) / std::sqrt(0.5 * (std::pow(R::rnorm(0, 1), 2) + std::pow(R::rnorm(0, 1), 2))));
                cc += 1;
                if (cc == 1000) {
                  std::string message = "No draws match all sign restrictions:\n  Iterations: " + std::to_string(cc);
                  Rcpp::stop(message);
                }
              }
            }
          }
          
          aa += 1;
          cc = 0;
          
        } else {
          
          A_test(j, i) = pA(j, i, 2);
          
        }
      }
    }
    
    if (nH > 0) {
      H_test = arma::inv(A_test);
      for (arma::uword i = 0; i < ncol; ++i) {
        for (arma::uword j = 0; j < nrow; ++j) {
          if ((std::isfinite(pH(j, i, 0))) && (std::isfinite(pH(j, i, 1)))) {
            if (pH(j, i, 0) == 0.0) {
              if (((pH(j, i, 1) > 0.0) && (H_test(j, i) < 0.0)) || ((pH(j, i, 1) < 0.0) && (H_test(j, i) > 0.0))) {
                sign_test += 1;
              }
            }
            if (pH(j, i, 0) == 2.0) {
              if (((pH(j, i, 1) == 1.0) && (H_test(j, i) < 1.0)) || ((pH(j, i, 1) == -1.0) && (H_test(j, i) > -1.0))) {
                sign_test += 1;
              }
            }
            if (pH(j, i, 0) == 3.0) {
              if (((pH(j, i, 1) == 1.0) && ((H_test(j, i) < 0.0) || (H_test(j, i) > 1.0))) || ((pH(j, i, 1) == -1.0) && ((H_test(j, i) > 0.0) || (H_test(j, i) < -1.0)))) {
                sign_test += 1;
              }
            }
          }
        }
      }
    }
    
    if ((sign_test == 0) && (std::isfinite(pdetA(0, 0, 0))) && (std::isfinite(pdetA(0, 0, 1)))) {
      detA_test = arma::det(A_test);
      if (pdetA(0, 0, 0) == 0.0) {
        if (((pdetA(0, 0, 1) > 0.0) && (detA_test < 0.0)) || ((pdetA(0, 0, 1) < 0.0) && (detA_test > 0.0))) {
          sign_test += 1;
        }
      }
      if (pdetA(0, 0, 0) == 2.0) {
        if (((pdetA(0, 0, 1) == 1.0) && (detA_test < 1.0)) || ((pdetA(0, 0, 1) == -1.0) && (detA_test > -1.0))) {
          sign_test += 1;
        }
      }
      if (pdetA(0, 0, 0) == 3.0) {
        if (((pdetA(0, 0, 1) == 1.0) && ((detA_test < 0.0) || (detA_test > 1.0))) || ((pdetA(0, 0, 1) == -1.0) && ((detA_test > 0.0) || (detA_test < -1.0)))) {
          sign_test += 1;
        }
      }
    }
    
  }
  
  return A_test;
  
}

// Start Random-Walk Metropolis-Hastings Algorithm.
arma::field<arma::cube> MH(const arma::mat& y1, const arma::mat& x1, const arma::uword nlags, const arma::mat& omega, const arma::mat& somega, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& pP, const arma::mat& pP_sig, const arma::cube& pR_sig, const arma::mat& kappa1, arma::mat A_old, const arma::uword itr, const arma::uword burn, const arma::mat& scale1) {
  const arma::uword pA_nrow = pA.n_rows, pA_ncol = pA.n_cols;
  
  A_old = proposal_function(A_old, pA, pdetA, pH, scale1);

  arma::mat B_old(((pA_ncol * nlags) + 1), pA_nrow, arma::fill::zeros), zeta_old(pA_nrow, pA_ncol, arma::fill::zeros);
  
  arma::mat A_test(pA_nrow, pA_ncol, arma::fill::zeros), B_test(((pA_ncol * nlags) + 1), pA_nrow, arma::fill::zeros), zeta_test(pA_nrow, pA_ncol, arma::fill::zeros);
  
  arma::cube A_chain(pA_nrow, pA_ncol, (itr - burn)), B_chain(((pA_ncol * nlags) + 1), pA_ncol, (itr - burn)), zeta_chain(pA_nrow, pA_ncol, (itr - burn));
  std::fill(A_chain.begin(), A_chain.end(), Rcpp::NumericVector::get_na());
  std::fill(B_chain.begin(), B_chain.end(), Rcpp::NumericVector::get_na());
  std::fill(zeta_chain.begin(), zeta_chain.end(), Rcpp::NumericVector::get_na());
  
  arma::cube pR(((pA_ncol * nlags) + 1), pA_ncol, pA_ncol, arma::fill::zeros);
  
  const arma::mat temp0 = x1.t() * x1 + pP_sig;
  const arma::mat temp1 = x1.t() * y1 + pP_sig * pP;
  const arma::mat temp2 = y1.t() * y1 + pP.t() * pP_sig * pP;
  const arma::mat Phi0 = arma::solve(temp0, temp1);
  const arma::mat temp4 = temp2 - Phi0.t() * temp1;
  arma::mat temp5(((pA_ncol * nlags) + 1), 1, arma::fill::zeros);
  
  arma::vec nR(pA_ncol, arma::fill::zeros);
  for (arma::uword i = 0; i < pA_ncol; ++i) {
    if (arma::any(arma::vectorise(pR_sig.slice(i)) > 0.0)) {
      nR(i) = 1.0;
    }
  }
  
  if (arma::all(nR == 0.0)) {
    B_old = Phi0 * A_old;
    zeta_old = arma::diagmat(A_old.t() * temp4 * A_old);
  } else {
    for (arma::uword i = 0; i < pA_ncol; ++i) {
      if (nR(i) == 0.0) {
        B_old.col(i) = Phi0 * A_old.col(i);
        zeta_old(i, i) = arma::as_scalar(arma::trans(A_old.col(i)) * temp4 * A_old.col(i));
      } else {
        for (arma::uword j = 0; j < pA_nrow; ++j) {
          if (std::isfinite(pA(j, i, 6))) {
            pR(j, i, i) = A_old(j, i);
          }
        }
        temp5 = pR_sig.slice(i) * pR.slice(i).col(i);
        B_old.col(i) = arma::inv((temp0 + pR_sig.slice(i))) * ((temp1 * A_old.col(i)) + temp5);
        zeta_old(i, i) = arma::as_scalar(((A_old.col(i).t() * temp2 * A_old.col(i)) + (pR.slice(i).col(i).t() * temp5)) - (arma::trans((temp1 * A_old.col(i)) + temp5) * arma::inv(temp0 + pR_sig.slice(i)) * ((temp1 * A_old.col(i)) + temp5)));
      }
    }
  }
  
  double post_A_old = post_A_function(A_old, pA, pdetA, pH, kappa1, y1, omega, zeta_old, somega), post_A_test = 0.0, accept = 0.0, naccept = 0.0, ru = 0.0;
  
  //start Metropolis-Hastings algorithm
  for (arma::uword c = 0; c < itr; ++c) {
    
    //check for interruptions
    if (c % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    //random draw proposals
    A_test = proposal_function(A_old, pA, pdetA, pH, scale1);
    
    if (arma::all(nR == 0.0)) {
      B_test = Phi0 * A_test;
      zeta_test = arma::diagmat(A_test.t() * temp4 * A_test);
    } else {
      for (arma::uword i = 0; i < pA_ncol; ++i) {
        if (nR(i) == 0.0) {
          B_test.col(i) = Phi0 * A_test.col(i);
          zeta_test(i, i) = arma::as_scalar(arma::trans(A_test.col(i)) * temp4 * A_test.col(i));
        } else {
          for (arma::uword j = 0; j < pA_nrow; ++j) {
            if (std::isfinite(pA(j, i, 6))) {
              pR(j, i, i) = A_test(j, i);
            }
          }
          temp5 = pR_sig.slice(i) * pR.slice(i).col(i);
          B_test.col(i) = arma::inv((temp0 + pR_sig.slice(i))) * ((temp1 * A_test.col(i)) + temp5);
          zeta_test(i, i) = arma::as_scalar(((A_test.col(i).t() * temp2 * A_test.col(i)) + (pR.slice(i).col(i).t() * temp5)) - (arma::trans((temp1 * A_test.col(i)) + temp5) * arma::inv(temp0 + pR_sig.slice(i)) * ((temp1 * A_test.col(i)) + temp5)));
        }
      }
    }
    
    //compute posterior density
    post_A_test = post_A_function(A_test, pA, pdetA, pH, kappa1, y1, omega, zeta_test, somega);
    //acceptance value
    accept = std::exp(post_A_test - post_A_old);
    //threshold to determine whether to accept proposals
    ru = R::runif(0.0, 1.0);
    
    //determine if the proposals will be kept
    if (ru <= accept) {
      A_old = A_test;
      B_old = B_test;
      zeta_old = zeta_test;
      post_A_old = post_A_test;
      if (c >= burn) {
        naccept += 1.0;
      }
    }
    
    if (c >= burn) {
      A_chain.slice((c - burn)) = A_old;
      B_chain.slice((c - burn)) = B_old;
      zeta_chain.slice((c - burn)) = zeta_old;
    }
    
  }
  
  //compute acceptance rate
  double total = double (itr - burn);
  arma::cube accept_rate(1, 1, 1);
  accept_rate(0, 0, 0) = (naccept / total);
  
  arma::field<arma::cube> list2(4);
  list2(0) = accept_rate;
  list2(1) = A_chain;
  list2(2) = B_chain;
  list2(3) = zeta_chain;
  
  return list2;
  
}

// Thin Markov Chains.
arma::cube thin_function(const arma::cube& chain, const arma::uword thin) {
  const arma::uword nrow = chain.n_rows, ncol = chain.n_cols;
  const double totals = double (chain.n_slices), totalt = double (thin);
  const arma::uword nsli = arma::uword (std::floor((totals / totalt)));
  arma::cube chain1(nrow, ncol, nsli);
  std::fill(chain1.begin(), chain1.end(), Rcpp::NumericVector::get_na());
  
  //store a fraction of the estimates
  for (arma::uword i = 0; i < nsli; ++i) {
    chain1.slice(i) = chain.slice((thin * (i + 1) - 1));
  }
  return chain1;
}

// Process Raw Results.
arma::cube results_function(const arma::cube& raw, const double ci) {
  const arma::uword nrow = raw.n_rows, ncol = raw.n_cols, nsli = raw.n_slices;
  
  arma::mat temp(nrow, nsli), temp1(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());
  
  arma::cube results1(nrow, ncol, 3);
  std::fill(results1.begin(), results1.end(), Rcpp::NumericVector::get_na());
  
  const double total = double (nsli);
  const arma::uword lb = arma::uword (std::round(total * (1.0 - ci)));
  const arma::uword ub = arma::uword (std::round((total * ci) - 1.0));
  
  for (arma::uword j = 0; j < ncol; ++j) {
    
    temp = raw(arma::span::all, arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    results1.slice(0).col(j) = temp1.col(lb);
    results1.slice(1).col(j) = arma::median(temp1, 1);
    results1.slice(2).col(j) = temp1.col(ub);
    
  }
  
  return results1;
  
}

// Estimate Density Coordinates.
Rcpp::List den_function(const arma::cube& raw, const arma::cube& priors) {
  
  double ll = 0.0, from0 = 0.0, to0 = 0.0, by0 = 0.0, from1 = 0.0, to1 = 0.0, by1 = 0.0, t1 = 0.0;
  const double crr_lb = 0.001, nbins_d = 300.0, nbins_d1 = 300.0, nsli_d = double (raw.n_slices);
  
  const arma::uword nrow = raw.n_rows, ncol = raw.n_cols, nsli = raw.n_slices, nbins = 300, nbins1 = 300;
  
  arma::uword raw_fi = 0, raw_ti = 0;
  
  arma::vec breaks0(nbins), counts0(nbins), counts1(nbins1), t2(nsli), raw_temp(nsli);
  std::fill(breaks0.begin(), breaks0.end(), Rcpp::NumericVector::get_na());
  std::fill(counts0.begin(), counts0.end(), Rcpp::NumericVector::get_na());
  std::fill(counts1.begin(), counts1.end(), Rcpp::NumericVector::get_na());
  std::fill(t2.begin(), t2.end(), Rcpp::NumericVector::get_na());
  std::fill(raw_temp.begin(), raw_temp.end(), Rcpp::NumericVector::get_na());
  
  arma::cube hori(nrow, ncol, (nbins1 + 4)), vert(nrow, ncol, (nbins1 + 4));
  std::fill(hori.begin(), hori.end(), Rcpp::NumericVector::get_na());
  std::fill(vert.begin(), vert.end(), Rcpp::NumericVector::get_na());
  
  for (arma::uword i = 0; i < nrow; ++i) {
    for (arma::uword j = 0; j < ncol; ++j) {
      if (std::isfinite(priors(i, j, 0))) {
        t1 = raw(i, j, 0);
        t2(arma::span::all) = arma::vectorise(raw(arma::span(i), arma::span(j), arma::span::all));
        if (arma::any(t2 != t1)) {
          
          raw_temp(arma::span::all) = arma::sort(arma::vectorise(raw(arma::span(i), arma::span(j), arma::span::all)));
          
          from0 = arma::min(raw_temp);
          to0 = arma::max(raw_temp);
          by0 = (to0 - from0) / nbins_d;
          breaks0(0) = from0 + (by0 / 2.0);
          for (arma::uword k = 1; k < nbins; ++k) {
            ll = double (k);
            breaks0(k) = breaks0(0) + (by0 * ll);
          }
          
          counts0(arma::span::all) = arma::conv_to<arma::vec>::from(arma::hist(raw_temp, nbins_d));
          
          from1 = raw_temp(arma::as_scalar(arma::find(raw_temp <= arma::min(breaks0.elem(arma::find(counts0 > (arma::max(counts0) * crr_lb)))), 1, "last")));
          to1 = raw_temp(arma::as_scalar(arma::find(raw_temp >= arma::max(breaks0.elem(arma::find(counts0 > (arma::max(counts0) * crr_lb)))), 1, "first")));
          by1 = (to1 - from1) / nbins_d1;
          
          hori(i, j, 2) = from1 + (by1 / 2.0);
          for (arma::uword k = 3; k < (nbins1 + 2); ++k) {
            ll = double (k - 2);
            hori(i, j, k) = hori(i, j, 2) + (by1 * ll);
          }
          hori(i, j, 0) = from0;
          hori(i, j, 1) = from0;
          hori(i, j, (nbins1 + 2)) = to0;
          hori(i, j, (nbins1 + 3)) = to0;
          
          raw_fi = arma::as_scalar(arma::find(raw_temp <= from1, 1, "last"));
          raw_ti = arma::as_scalar(arma::find(raw_temp >= to1, 1, "first"));
          
          counts1 = arma::conv_to<arma::vec>::from(arma::hist(raw_temp(arma::span(raw_fi, raw_ti)), nbins_d1));
          
          vert.subcube(i, j, 0, i, j, (nbins1 + 3)).fill(0.0);
          vert.subcube(i, j, 2, i, j, (nbins1 + 1)) = counts1 * (1.0 / (nsli_d * by1));
          vert(i, j, 1) = vert(i, j, 2);
          vert(i, j, (nbins1 + 2)) = vert(i, j, (nbins1 + 1));
          
        }
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("hori") = hori,
                            Rcpp::Named("vert") = vert
                            );
  
}

// Estimate Historical Decompositions and Process Raw Results.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
arma::cube hd_estimates(const arma::cube& A_chain, const arma::cube& B_chain, const arma::mat& y1, const arma::mat& x1, const arma::uword nlags, const double ci) {
  const arma::uword nvar = A_chain.n_cols, nsli = A_chain.n_slices;
  const arma::uword nrow = y1.n_rows, ncol = (nvar * nvar);

  arma::cube HD_chain(nrow, ncol, nsli);
  std::fill(HD_chain.begin(), HD_chain.end(), Rcpp::NumericVector::get_na());

  arma::mat HD_draw((nrow + nlags), (nvar * nvar), arma::fill::zeros);
  arma::mat u1(nrow, nvar), H_draw(nvar, nvar), Phi_draw(((nvar * nlags) + 1), nvar);
  std::fill(u1.begin(), u1.end(), Rcpp::NumericVector::get_na());
  std::fill(H_draw.begin(), H_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(Phi_draw.begin(), Phi_draw.end(), Rcpp::NumericVector::get_na());

  const arma::mat y1_temp = arma::reverse(y1);
  const arma::mat x1_temp = arma::reverse(x1);

  for (arma::uword t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //historical decomposition (HD_chain, U_chain, H_chain, Phi_chain)
    u1 = (y1_temp * A_chain.slice(t)) - (x1_temp * B_chain.slice(t));
    H_draw = arma::inv(A_chain.slice(t));
    Phi_draw = B_chain.slice(t) * H_draw;
    for (arma::uword i = 0; i < nvar; ++i) {
      HD_draw.submat(0, (nvar * i), (nrow - 1), ((nvar * (i + 1)) - 1)) = u1.col(i) * H_draw.row(i);
      for (arma::uword j = (nrow - 1); j >= 1; --j) {
        HD_draw.submat((j - 1), (nvar * i), (j - 1), ((nvar * (i + 1)) - 1)) += arma::vectorise(HD_draw.submat(j, (nvar * i), (j + nlags - 1), ((nvar * (i + 1)) - 1)).t()).t() * Phi_draw(arma::span(0, ((nvar * nlags) - 1)), arma::span::all);
      }
    }
    HD_chain.slice(t) = HD_draw(arma::span(0, (nrow - 1)), arma::span::all);
  }

  //process HD_chain
  arma::mat temp(1, nsli), temp1(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());

  arma::cube HD(nrow, ncol, 3);
  std::fill(HD.begin(), HD.end(), Rcpp::NumericVector::get_na());

  const double total = double (nsli);
  const arma::uword lb = arma::uword (std::round(total * (1.0 - ci)));
  const arma::uword ub = arma::uword (std::round((total * ci) - 1.0));

  for (arma::uword j = 0; j < ncol; ++j) {

    for (arma::uword i = 0; i < nrow; ++i) {
      temp = HD_chain.tube(i, j);
      temp1.row(i) = arma::sort(temp, "ascend", 1);
    }
    HD.slice(0).col(j) = arma::reverse(temp1.col(lb));
    HD.slice(1).col(j) = arma::reverse(arma::median(temp1, 1));
    HD.slice(2).col(j) = arma::reverse(temp1.col(ub));

  }

  return HD;
}

// Estimate Impulse Responses and Process Raw Results.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
arma::cube irf_estimates(const arma::cube& A_chain, const arma::cube& B_chain, const arma::uword nlags, const arma::uword h, const bool acc, const double ci) {
  const arma::uword nvar = A_chain.n_cols, nsli = A_chain.n_slices;
  const arma::uword nrow = (h + 1), ncol = (nvar * nvar);

  arma::cube IRF_chain(nrow, ncol, nsli);
  std::fill(IRF_chain.begin(), IRF_chain.end(), Rcpp::NumericVector::get_na());

  arma::mat H_draw(nvar, nvar), IRF_draw((h + nlags), (nvar * nvar)), Phi_draw(((nvar * nlags) + 1), nvar);
  std::fill(H_draw.begin(), H_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(IRF_draw.begin(), IRF_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(Phi_draw.begin(), Phi_draw.end(), Rcpp::NumericVector::get_na());

  for (arma::uword t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //impulse responses (H_chain, IRF_chain, Phi_chain)
    H_draw = arma::inv(A_chain.slice(t));
    Phi_draw = B_chain.slice(t) * H_draw;
    IRF_draw.fill(0.0);
    for (arma::uword i = 0; i < nvar; ++i) {
      IRF_draw.submat((nrow - 1), (nvar * i), (nrow - 1), ((nvar * (i + 1)) - 1)) = H_draw.row(i);
      for (arma::uword j = (nrow - 1); j >= 1; --j) {
        IRF_draw.submat((j - 1), (nvar * i), (j - 1), ((nvar * (i + 1)) - 1)) = arma::vectorise(IRF_draw.submat(j, (nvar * i), (j + nlags - 1), ((nvar * (i + 1)) - 1)).t()).t() * Phi_draw(arma::span(0, ((nvar * nlags) - 1)), arma::span::all);
      }
    }
    if (acc) {
      for (arma::uword j = (nrow - 1); j >= 1; --j) {
        IRF_draw.row((j - 1)) += IRF_draw.row(j);
      }
    }
    IRF_chain.slice(t) = IRF_draw(arma::span(0, (nrow - 1)), arma::span::all);

  }

  //process IRF_chain
  arma::mat temp(1, nsli), temp1(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());

  arma::cube IRF(nrow, ncol, 3);
  std::fill(IRF.begin(), IRF.end(), Rcpp::NumericVector::get_na());

  const double total = double (nsli);
  const arma::uword lb = arma::uword (std::round(total * (1.0 - ci)));
  const arma::uword ub = arma::uword (std::round((total * ci) - 1.0));

  for (arma::uword j = 0; j < ncol; ++j) {

    for (arma::uword i = 0; i < nrow; ++i) {
      temp = IRF_chain.tube(i, j);
      temp1.row(i) = arma::sort(temp, "ascend", 1);
    }
    IRF.slice(0).col(j) = arma::reverse(temp1.col(lb));
    IRF.slice(1).col(j) = arma::reverse(arma::median(temp1, 1));
    IRF.slice(2).col(j) = arma::reverse(temp1.col(ub));

  }

  return IRF;
}

// Estimate Forecast Error Variance Decompositions and Process Raw Results.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
arma::cube fevd_estimates(const arma::cube& A_chain, const arma::cube& B_chain, const arma::cube& D_chain, const arma::uword nlags, const arma::uword h, const bool acc, const double ci) {
  const arma::uword nvar = A_chain.n_cols, nsli = A_chain.n_slices;
  const arma::uword nrow = (h + 1), ncol = (nvar * nvar);

  arma::cube FEVD_chain(nrow, ncol, nsli);
  std::fill(FEVD_chain.begin(), FEVD_chain.end(), Rcpp::NumericVector::get_na());

  arma::mat H_draw(nvar, nvar), FEVD_draw((h + nlags), (nvar * nvar)), IRF_draw((h + nlags), (nvar * nvar)), Phi_draw(((nvar * nlags) + 1), nvar);
  std::fill(H_draw.begin(), H_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(FEVD_draw.begin(), FEVD_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(IRF_draw.begin(), IRF_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(Phi_draw.begin(), Phi_draw.end(), Rcpp::NumericVector::get_na());

  for (arma::uword t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //impulse responses (H_chain, IRF_chain, Phi_chain)
    H_draw = arma::inv(A_chain.slice(t));
    Phi_draw = B_chain.slice(t) * H_draw;
    IRF_draw.fill(0.0);
    for (arma::uword i = 0; i < nvar; ++i) {
      IRF_draw.submat((nrow - 1), (nvar * i), (nrow - 1), ((nvar * (i + 1)) - 1)) = H_draw.row(i);
      for (arma::uword j = (nrow - 1); j >= 1; --j) {
        IRF_draw.submat((j - 1), (nvar * i), (j - 1), ((nvar * (i + 1)) - 1)) = arma::vectorise(IRF_draw.submat(j, (nvar * i), (j + nlags - 1), ((nvar * (i + 1)) - 1)).t()).t() * Phi_draw(arma::span(0, ((nvar * nlags) - 1)), arma::span::all);
      }
    }
    if (acc) {
      for (arma::uword j = (nrow - 1); j >= 1; --j) {
        IRF_draw.row((j - 1)) += IRF_draw.row(j);
      }
    }

    FEVD_draw.fill(0.0);
    for (arma::uword i = 0; i < nvar; ++i) {
      FEVD_draw.submat((nrow - 1), (nvar * i), (nrow - 1), ((nvar * (i + 1)) - 1)) = (D_chain(i, i, t)) * (IRF_draw.submat((nrow - 1), (nvar * i), (nrow - 1), ((nvar * (i + 1)) - 1)) % IRF_draw.submat((nrow - 1), (nvar * i), (nrow - 1), ((nvar * (i + 1)) - 1)));
      for (arma::uword j = (nrow - 1); j >= 1; --j) {
        FEVD_draw.submat((j - 1), (nvar * i), (j - 1), ((nvar * (i + 1)) - 1)) = FEVD_draw.submat(j, (nvar * i), j, ((nvar * (i + 1)) - 1)) + ((D_chain(i, i, t)) * (IRF_draw.submat((j - 1), (nvar * i), (j - 1), ((nvar * (i + 1)) - 1)) % IRF_draw.submat((j - 1), (nvar * i), (j - 1), ((nvar * (i + 1)) - 1))));
      }
    }

    FEVD_chain.slice(t) = FEVD_draw(arma::span(0, (nrow - 1)), arma::span::all);
  }

  //process FEVD_chain
  arma::mat temp(1, nsli), temp1(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());

  arma::cube FEVD(nrow, ncol, 3);
  std::fill(FEVD.begin(), FEVD.end(), Rcpp::NumericVector::get_na());

  const double total = double (nsli);
  const arma::uword lb = arma::uword (std::round(total * (1.0 - ci)));
  const arma::uword ub = arma::uword (std::round((total * ci) - 1.0));

  for (arma::uword j = 0; j < ncol; ++j) {

    for (arma::uword i = 0; i < nrow; ++i) {
      temp = FEVD_chain.tube(i, j);
      temp1.row(i) = arma::sort(temp, "ascend", 1);
    }
    FEVD.slice(0).col(j) = arma::reverse(temp1.col(lb));
    FEVD.slice(1).col(j) = arma::reverse(arma::median(temp1, 1));
    FEVD.slice(2).col(j) = arma::reverse(temp1.col(ub));

  }

  return FEVD;
}

// Estimate Phi and Process Raw Results.
arma::cube phi_estimates(const arma::cube& A_chain, const arma::cube& B_chain, const arma::uword nlags, const double ci) {
  const arma::uword nvar = A_chain.n_cols, nsli = A_chain.n_slices;
  const arma::uword nrow = ((nvar * nlags) + 1), ncol = nvar;

  arma::cube Phi_chain(nrow, ncol, nsli);
  std::fill(Phi_chain.begin(), Phi_chain.end(), Rcpp::NumericVector::get_na());

  for (arma::uword t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }

    Phi_chain.slice(t) = B_chain.slice(t) * arma::inv(A_chain.slice(t));

  }

  //process Phi_chain
  arma::mat temp(nrow, nsli), temp1(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());

  arma::cube Phi(nrow, ncol, 3);
  std::fill(Phi.begin(), Phi.end(), Rcpp::NumericVector::get_na());

  const double total = double (nsli);
  const arma::uword lb = arma::uword (std::round(total * (1.0 - ci)));
  const arma::uword ub = arma::uword (std::round((total * ci) - 1.0));

  for (arma::uword j = 0; j < ncol; ++j) {

    temp = Phi_chain(arma::span::all, arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    Phi.slice(0).col(j) = temp1.col(lb);
    Phi.slice(1).col(j) = arma::median(temp1, 1);
    Phi.slice(2).col(j) = temp1.col(ub);

  }

  return Phi;
}

// Estimate H and Process Raw Results.
Rcpp::List h_estimates(const arma::cube& A_chain, const arma::cube& pH, const double ci) {
  const arma::uword nvar = A_chain.n_cols, nsli = A_chain.n_slices;

  arma::cube H_chain(nvar, nvar, nsli);
  std::fill(H_chain.begin(), H_chain.end(), Rcpp::NumericVector::get_na());

  for (arma::uword t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }

    H_chain.slice(t) = arma::inv(A_chain.slice(t));

  }

  return Rcpp::List::create(Rcpp::Named("H_den") = den_function(H_chain, pH),
                            Rcpp::Named("H") = results_function(H_chain, ci),
                            Rcpp::Named("H_chain") = H_chain);
}

// Estimate det(A) and Process Raw Results.
Rcpp::List deta_estimates(const arma::cube& A_chain, const arma::cube& pdetA, const double ci) {
  const arma::uword nsli = A_chain.n_slices;

  arma::cube detA_chain(1, 1, nsli);
  std::fill(detA_chain.begin(), detA_chain.end(), Rcpp::NumericVector::get_na());

  for (arma::uword t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }

    detA_chain.slice(t) = arma::det(A_chain.slice(t));

  }

  return Rcpp::List::create(Rcpp::Named("detA_den") = den_function(detA_chain, pdetA),
                            Rcpp::Named("detA") = results_function(detA_chain, ci),
                            Rcpp::Named("detA_chain") = detA_chain);
}

// Estimate the Parameters of the BH_SBVAR Model.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List MAIN(const arma::mat& y1, const arma::mat& x1, const arma::mat& omega, const arma::mat& somega, const arma::uword nlags, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& pP, const arma::mat& pP_sig, const arma::cube& pR_sig, const arma::mat& kappa1, const arma::mat& A_start, const arma::uword itr, const arma::uword burn, const arma::uword thin, const arma::mat& scale1, const double ci) {
  const arma::uword pA_ncol = pA.n_cols, B_nrow = ((y1.n_cols * nlags) + 1);
  const double totals = double (itr - burn), totalt = double (thin), y1_nrow = double (y1.n_rows);
  const arma::uword nsli = arma::uword (std::floor((totals / totalt)));
  //start Metropolis-Hastings algorithm
  arma::field<arma::cube> list1 = MH(y1, x1, nlags, omega, somega, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, A_start, itr, burn, scale1);
  double accept_rate = list1(0)(0, 0, 0);
  
  arma::cube A_chain(pA_ncol, pA_ncol, nsli), B_chain(B_nrow, pA_ncol, nsli), zeta_chain(pA_ncol, pA_ncol, nsli), D_chain(pA_ncol, pA_ncol, nsli, arma::fill::zeros);
  std::fill(A_chain.begin(), A_chain.end(), Rcpp::NumericVector::get_na());
  std::fill(B_chain.begin(), B_chain.end(), Rcpp::NumericVector::get_na());
  std::fill(zeta_chain.begin(), zeta_chain.end(), Rcpp::NumericVector::get_na());

  if (thin > 1) {
    A_chain(arma::span::all, arma::span::all, arma::span::all) = thin_function(arma::cube (list1(1)), thin);
    B_chain(arma::span::all, arma::span::all, arma::span::all) = thin_function(arma::cube (list1(2)), thin);
    zeta_chain(arma::span::all, arma::span::all, arma::span::all) = thin_function(arma::cube (list1(3)), thin);
  } else {
    A_chain(arma::span::all, arma::span::all, arma::span::all) = list1(1);
    B_chain(arma::span::all, arma::span::all, arma::span::all) = list1(2);
    zeta_chain(arma::span::all, arma::span::all, arma::span::all) = list1(3);
  }
  
  arma::mat taustar(pA_ncol, pA_ncol), kappastar(pA_ncol, pA_ncol), M(B_nrow, pA_ncol);
  std::fill(taustar.begin(), taustar.end(), Rcpp::NumericVector::get_na());
  std::fill(kappastar.begin(), kappastar.end(), Rcpp::NumericVector::get_na());
  std::fill(M.begin(), M.end(), Rcpp::NumericVector::get_na());
  
  arma::mat Dinv_draw(pA_ncol, pA_ncol, arma::fill::zeros);
  
  arma::cube temp0(B_nrow, B_nrow, pA_ncol, arma::fill::zeros);
  for (arma::uword i = 0; i < pA_ncol; ++i) {
    temp0.slice(i) = arma::chol(arma::inv((x1.t() * x1) + pP_sig + pR_sig.slice(i))).t();
  }
  kappastar = arma::diagmat((kappa1 + ((y1_nrow / 2.0) * arma::ones(1, pA_ncol))));
  
  //estimate matrix B
  for (arma::uword t = 0; t < nsli; ++t) {
    
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    taustar = (arma::diagmat(kappa1) * arma::diagmat((A_chain.slice(t).t() * somega * A_chain.slice(t)))) + arma::diagmat((0.5 * zeta_chain.slice(t)));

    for (arma::uword i = 0; i < pA_ncol; ++i) {
      Dinv_draw(i, i) = R::rgamma((kappastar(i, i)), (1.0 / (taustar(i, i))));
    }
    D_chain.slice(t) = arma::inv(Dinv_draw);
    
    for (arma::uword i = 0; i < pA_ncol; ++i) {
      M = temp0.slice(i) * arma::randn(B_nrow, pA_ncol) * arma::inv(arma::sqrt(Dinv_draw));
      B_chain(arma::span::all, arma::span(i), arma::span(t)) += M.col(i);
    }
    
  }
  
  //construct pR
  arma::cube pR(((pA_ncol * nlags) + 1), pA_ncol, pA_ncol, arma::fill::zeros);
  for (arma::uword i = 0; i < pA_ncol; ++i) {
    for (arma::uword j = 0; j < pA_ncol; ++j) {
      if (std::isfinite(pA(j, i, 6))) {
        pR(j, i, i) = pA(j, i, 2);
      }
    }
  }
  
  Rcpp::List list2 = deta_estimates(A_chain, pdetA, ci);
  Rcpp::List list3 = h_estimates(A_chain, pH, ci);
  
  Rcpp::List list4(28);
  list4.names() = 
    Rcpp::CharacterVector(
      {"accept_rate", "y", "x", "nlags", "pA", "pdetA", "pH", "pP", "pP_sig", "pR", "pR_sig",
       "tau1", "kappa1", "A_start", "A", "detA", "H", "B", "Phi", "D",
       "A_den", "detA_den", "H_den", "A_chain", "B_chain", "D_chain", "detA_chain", "H_chain"}
    );
  
  list4["accept_rate"] = accept_rate;
  
  list4["y"] = y1;
  list4["x"] = x1;
  
  list4["nlags"] = nlags;
  
  list4["pA"] = pA;
  list4["pdetA"] = pdetA;
  list4["pH"] = pH;
  list4["pP"] = pP;
  list4["pP_sig"] = pP_sig;
  list4["pR"] = pR;
  list4["pR_sig"] = pR_sig;
  list4["tau1"] = arma::diagmat(kappa1) * arma::diagmat(pA.slice(2).t() * somega * pA.slice(2));
  list4["kappa1"] = arma::diagmat(kappa1);
  
  list4["A_start"] = A_start;
  
  list4["A"] = results_function(A_chain, ci);
  list4["detA"] = list2["detA"];
  list4["H"] = list3["H"];
  list4["B"] = results_function(B_chain, ci);
  list4["Phi"] = phi_estimates(A_chain, B_chain, nlags, ci);
  list4["D"] = results_function(D_chain, ci);
  
  list4["A_den"] = den_function(A_chain, pA);
  list4["detA_den"] = list2["detA_den"];
  list4["H_den"] = list3["H_den"];
  
  list4["A_chain"] = A_chain;
  list4["B_chain"] = B_chain;
  list4["D_chain"] = D_chain;
  list4["detA_chain"] = list2["detA_chain"];
  list4["H_chain"] = list3["H_chain"];

  return list4;
}


