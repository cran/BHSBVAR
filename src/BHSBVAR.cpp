



//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//[[Rcpp::plugins(cpp11)]]


// Non-Central T-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_nonc_t(const double a1, const double p1, const double sigma1, const double nu, const double lam1) {
  double den = R::dnt(((a1 - p1) / sigma1), nu, lam1, 0) / sigma1;
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

// Sum the Log of Prior Densities.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double sum_log_prior_densities(const arma::mat& A_test, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH) {
  int nrow = int (pA.n_rows), ncol = int (pA.n_cols);
  arma::mat H_test = arma::inv(A_test);
  double sum_log_priors = 0.0, detA_test = arma::det(A_test);
  //compute sum of log priors
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      // if pA(i,j,0) == 0, symmetric t-distribution
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
    }
  }
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
  return sum_log_priors;
}

// Likelihood Function.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double likelihood_function(const arma::mat& A_test, const arma::mat& kappa1, const arma::mat& y1, const arma::mat& omega, const arma::mat& zeta, const arma::mat& somega) {
  int kappa_ncol = int (kappa1.n_cols);
  double ynrow = double (y1.n_rows);
  double lik_A_numerator = (ynrow / 2.0) * std::log(arma::det(A_test.t() * omega * A_test));
  arma::mat tau = arma::diagmat(kappa1) * arma::diagmat(A_test.t() * somega * A_test);
  for (int i = 0; i < kappa_ncol; ++i) {
    if (kappa1(0, i) > 0.0) {
      lik_A_numerator += kappa1(0, i) * std::log(tau(i, i));
    }
  }
  double lik_A_denomenator = arma::as_scalar((kappa1 + ((ynrow / 2.0) * arma::ones(1, kappa_ncol))) * arma::log(arma::diagvec((2.0 / ynrow) * (tau + (arma::diagmat(zeta) / 2.0)))));
  double lik_A = lik_A_numerator - lik_A_denomenator;
  return lik_A;
}

// Posterior Density Function.
double post_A_function(const arma::mat& A_test, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& kappa1, const arma::mat& y1, const arma::mat& omega, const arma::mat& zeta, const arma::mat& somega) {
  double priors = sum_log_prior_densities(A_test, pA, pdetA, pH);
  double likelihood = likelihood_function(A_test, kappa1, y1, omega, zeta, somega);
  double posterior = priors + likelihood;
  return posterior;
}

// Genterate Proposal Values for A.
arma::mat proposal_function(const arma::mat& A_old, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& scale1) {
  
  int nrow = int (pA.n_rows), ncol = int (pA.n_cols);
  
  int nH = 0;
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < nrow; ++j) {
      if ((arma::is_finite(pH(j, i, 0))) && (pH(j, i, 0) == 0.0) && (arma::is_finite(pH(j, i, 1)))) {
        nH += 1;
      }
    }
  }
  
  arma::mat A_test(nrow, ncol), H_test(nrow, ncol);
  std::fill(A_test.begin(), A_test.end(), Rcpp::NumericVector::get_na());
  std::fill(H_test.begin(), H_test.end(), Rcpp::NumericVector::get_na());
  
  double detA_test = 0.0;
  
  int sign_test = 1, aa = 0, bb = 0, cc = 0;
  
  while (sign_test != 0) {
    
    sign_test = 0;
    aa = 0;
    bb += 1;
    
    if (bb == 1000) {
      std::string message = "No draws match all sign restrictions:\n  Iterations: " + std::to_string(bb);
      Rcpp::stop(message);
    }
    
    for (int i = 0; i < ncol; ++i) {
      for (int j = 0; j < nrow; ++j) {
        if (arma::is_finite(pA(j, i, 0))) {
          
          A_test(j, i) = A_old(j, i) + scale1(aa, aa) * (R::rnorm(0, 1) / std::sqrt(0.5 * (std::pow(R::rnorm(0, 1), 2) + std::pow(R::rnorm(0, 1), 2))));
          
          while ((pA(j, i, 0) == 0.0) && (arma::is_finite(pA(j, i, 1))) && (((pA(j, i, 1) > 0.0) & (A_test(j, i) < 0.0)) || ((pA(j, i, 1) < 0.0) & (A_test(j, i) > 0.0)))) {
            A_test(j, i) = A_old(j, i) + scale1(aa, aa) * (R::rnorm(0, 1) / std::sqrt(0.5 * (std::pow(R::rnorm(0, 1), 2) + std::pow(R::rnorm(0, 1), 2))));
            cc += 1;
            if (cc == 1000) {
              std::string message = "No draws match all sign restrictions:\n  Iterations: " + std::to_string(cc);
              Rcpp::stop(message);
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
      for (int i = 0; i < ncol; ++i) {
        for (int j = 0; j < nrow; ++j) {
          if ((arma::is_finite(pH(j, i, 0))) && (pH(j, i, 0) == 0.0) && (arma::is_finite(pH(j, i, 1))) && (((pH(j, i, 1) > 0.0) & (H_test(j, i) < 0.0)) || ((pH(j, i, 1) < 0.0) & (H_test(j, i) > 0.0)))) {
            sign_test += 1;
          }
        }
      }
    }
    
    if ((sign_test == 0) && (arma::is_finite(pdetA(0, 0, 0))) && (pdetA(0, 0, 0) == 0.0) && (arma::is_finite(pdetA(0, 0, 1)))) {
      detA_test = arma::det(A_test);
      if (((pdetA(0, 0, 1) > 0.0) & (detA_test < 0.0)) || ((pdetA(0, 0, 1) < 0.0) & (detA_test > 0.0))) {
        sign_test += 1;
      }
    }
    
  }
  
  return A_test;
  
}

// Start Random-Walk Metropolis-Hastigs Algorithm.
arma::field<arma::cube> MH(const arma::mat& y1, const arma::mat& x1, const int nlags, const arma::mat& omega, const arma::mat& somega, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& pP, const arma::mat& pP_sig, const arma::cube& pR_sig, const arma::mat& kappa1, arma::mat A_old, const int itr, const int burn, const arma::mat& scale1) {
  int pA_nrow = int (pA.n_rows), pA_ncol = int (pA.n_cols);
  
  A_old = proposal_function(A_old, pA, pdetA, pH, scale1);

  arma::mat B_old(((pA_ncol * nlags) + 1), pA_nrow, arma::fill::zeros), zeta_old(pA_nrow, pA_ncol, arma::fill::zeros);
  
  arma::mat A_test(pA_nrow, pA_ncol, arma::fill::zeros), B_test(((pA_ncol * nlags) + 1), pA_nrow, arma::fill::zeros), zeta_test(pA_nrow, pA_ncol, arma::fill::zeros);
  
  arma::cube A_chain(pA_nrow, pA_ncol, (itr - burn)), B_chain(((pA_ncol * nlags) + 1), pA_ncol, (itr - burn)), zeta_chain(pA_nrow, pA_ncol, (itr - burn));
  std::fill(A_chain.begin(), A_chain.end(), Rcpp::NumericVector::get_na());
  std::fill(B_chain.begin(), B_chain.end(), Rcpp::NumericVector::get_na());
  std::fill(zeta_chain.begin(), zeta_chain.end(), Rcpp::NumericVector::get_na());
  
  arma::mat Phi1(((pA_ncol * nlags) + 1), pA_ncol, arma::fill::zeros);
  
  arma::cube pR(((pA_ncol * nlags) + 1), pA_ncol, pA_ncol, arma::fill::zeros);
  
  const arma::mat temp0 = x1.t() * x1 + pP_sig;
  const arma::mat temp1 = x1.t() * y1 + pP_sig * pP;
  const arma::mat temp2 = y1.t() * y1 + pP.t() * pP_sig * pP;
  const arma::mat Phi0 = arma::solve(temp0, temp1);
  const arma::mat temp4 = temp2 - Phi0.t() * temp1;
  arma::mat temp5(((pA_ncol * nlags) + 1), pA_ncol, arma::fill::zeros);
  
  arma::vec nR(pA_ncol, arma::fill::zeros);
  for (int i = 0; i < pA_ncol; ++i) {
    if (arma::any(arma::vectorise(pR_sig.slice(i)) > 0.0)) {
      nR(i) = 1.0;
    }
  }
  
  if (arma::all(nR == 0.0)) {
    B_old = Phi0 * A_old;
    zeta_old = arma::diagmat(A_old.t() * temp4 * A_old);
  } else {
    for (int i = 0; i < pA_ncol; ++i) {
      if (nR(i) == 0.0) {
        B_old.col(i) = Phi0 * A_old.col(i);
        zeta_old(i, i) = arma::as_scalar(arma::trans(A_old.col(i)) * temp4 * A_old.col(i));
      } else {
        for (int j = 0; j < pA_nrow; ++j) {
          if (arma::is_finite(pA(j, i, 6))) {
            pR(j, i, i) = A_old(j, i);
          }
        }
        temp5 = pR_sig.slice(i) * pR.slice(i);
        Phi1 = arma::solve((temp0 + pR_sig.slice(i)), (temp1 + temp5));
        B_old.col(i) = Phi1 * A_old.col(i);
        zeta_old(i, i) = arma::as_scalar(arma::trans(A_old.col(i)) * ((temp2 + pR.slice(i).t() * temp5) - (Phi1.t() * (temp1 + temp5))) * A_old.col(i));
      }
    }
  }
  
  double post_A_old = post_A_function(A_old, pA, pdetA, pH, kappa1, y1, omega, zeta_old, somega), post_A_test = 0.0, accept = 0.0, naccept = 0.0, ru = 0.0;
  
  //start Metropolis-Hastings algorithm
  for (int c = 0; c < itr; ++c) {
    
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
      for (int i = 0; i < pA_ncol; ++i) {
        if (nR(i) == 0.0) {
          B_test.col(i) = Phi0 * A_test.col(i);
          zeta_test(i, i) = arma::as_scalar(arma::trans(A_test.col(i)) * temp4 * A_test.col(i));
        } else {
          for (int j = 0; j < pA_nrow; ++j) {
            if (arma::is_finite(pA(j, i, 6))) {
              pR(j, i, i) = A_test(j, i);
            }
          }
          temp5 = pR_sig.slice(i) * pR.slice(i);
          Phi1 = arma::solve((temp0 + pR_sig.slice(i)), (temp1 + temp5));
          B_test.col(i) = Phi1 * A_test.col(i);
          zeta_test(i, i) = arma::as_scalar(arma::trans(A_test.col(i)) * ((temp2 + pR.slice(i).t() * temp5) - (Phi1.t() * (temp1 + temp5))) * A_test.col(i));
        }
      }
    }
    
    //compute posterior density
    post_A_test = post_A_function(A_test, pA, pdetA, pH, kappa1, y1, omega, zeta_test, somega);
    //acceptance value
    accept = std::exp(post_A_test - post_A_old);
    //threshold to determine whether to accept proposals
    ru = R::runif(0, 1);
    
    //determing if the proposals will be kept
    if (ru <= accept) {
      A_old = A_test;
      B_old = B_test;
      zeta_old = zeta_test;
      post_A_old = post_A_test;
      naccept += 1.0; 
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
arma::cube thin_function(const arma::cube& chain, const int thin) {
  int nrow = int (chain.n_rows), ncol = int (chain.n_cols);
  double totals = double (chain.n_slices), totalt = double (thin);
  int nsli = int (std::floor((totals / totalt)));
  arma::cube chain1(nrow, ncol, nsli);
  std::fill(chain1.begin(), chain1.end(), Rcpp::NumericVector::get_na());
  
  //store a fraction of the estimates
  for (int i = 0; i < nsli; ++i) {
    chain1.slice(i) = chain.slice((thin * (i + 1) - 1));
  }
  return chain1;
}

// Process Raw Results.
arma::cube results_function(const arma::cube& raw, const double ci) {
  int nrow = int (raw.n_rows), ncol = int (raw.n_cols), nsli = int (raw.n_slices);
  
  arma::mat temp(nrow, nsli), temp1(nrow, 1);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());

  arma::cube results1(nrow, ncol, 3);
  std::fill(results1.begin(), results1.end(), Rcpp::NumericVector::get_na());
  
  double total = double (nsli);
  int lb = int (std::ceil(total * (1.0 - ci))) - 1;
  int ub = int (std::ceil(total * ci)) - 1;
  
  for (int j = 0; j < ncol; ++j) {
    
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
  
  //obtain environment containing function
  Rcpp::Environment stats_r("package:stats");
  //make function callable from c++
  Rcpp::Function density_r = stats_r["density"];
  
  int nrow = int (raw.n_rows), ncol = int (raw.n_cols), nsli = int (raw.n_slices);
  
  double t1 = 0.0;
  
  arma::vec t2(nsli);
  std::fill(t2.begin(), t2.end(), Rcpp::NumericVector::get_na());
  
  arma::cube hori(nrow, ncol, 512), vert(nrow, ncol, 512);
  std::fill(hori.begin(), hori.end(), Rcpp::NumericVector::get_na());
  std::fill(vert.begin(), vert.end(), Rcpp::NumericVector::get_na());
  
  Rcpp::List list1;
  
  int n = 0;
  double den = 0.0;
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (arma::is_finite(priors(i, j, 0))) {
        t1 = raw(i, j, 0);
        t2 = raw(arma::span(i), arma::span(j), arma::span(0, (nsli - 1)));
        if (arma::any(t2 != t1) == true) {
          list1 = density_r(Rcpp::_["x"] = raw(arma::span(i), arma::span(j), arma::span(0, (nsli - 1))), Rcpp::_["n"] = 512);
          for (int k = 0; k < 512; ++k) {
            hori(i, j, k) = Rcpp::as<Rcpp::NumericVector>(list1["x"])(k);
            vert(i, j, k) = Rcpp::as<Rcpp::NumericVector>(list1["y"])(k);
          }
          n = 0;
          if (priors(i, j, 1) == 1) {
            for (int k = 0; k < 512; ++k) {
              if (hori(i, j, k) < 0) {
                n += 1;
                hori(i, j, k) = 0.0;
                vert(i, j, k) = 0.0;
              }
            }
            if (n > 0) {
              den = 0.0;
              vert(i, j, (n - 1)) = vert(i, j, n) - ((vert(i, j, (n + 1)) - vert(i, j, n)) * (std::abs(hori(i, j, n) / (hori(i, j, (n + 1)) - hori(i, j, n)))));
              for (int k = n; k < 512; ++k) {
                den += (arma::sum(vert(i, j, k)) * (hori(i, j, (n + 1)) - hori(i, j, n)));
              }
              den += (vert(i, j, (n - 1)) * hori(i, j, n));
              vert(arma::span(i), arma::span(j), arma::span(0, (512 - 1))) /= den;
            }
          }
          if (priors(i, j, 1) == (-1)) {
            for (int k = 0; k < 512; ++k) {
              if (hori(i, j, k) > 0) {
                n += 1;
                hori(i, j, k) = 0.0;
                vert(i, j, k) = 0.0;
              }
            }
            if (n > 0) {
              den = 0;
              vert(i, j, (512 - n)) = vert(i, j, (512 - n - 1)) - ((vert(i, j, (512 - n - 2)) - vert(i, j, (512 - n - 1))) * (std::abs(hori(i, j, (512 - n-1)) / (hori(i, j, (512 - n - 1)) - hori(i, j, (512 - n - 2))))));
              for (int k = n; k < 512; ++k) {
                den += (arma::sum(vert(i, j, k)) * (hori(i, j, (512 - n - 1)) - hori(i, j, (512 - n - 2))));
              }
              den += (vert(i, j, (512 - n)) * (hori(i, j, (512 - n - 1)) * (-1.0)));
              vert(arma::span(i), arma::span(j), arma::span(0, (512 - 1))) /= den;
            }
          }
        } else {
          hori(arma::span(i), arma::span(j), arma::span(0, (512 - 1))).fill(t1);
        }
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("hori") = hori,
                            Rcpp::Named("vert") = vert);
}

// Line Plots.
void line_plots(const arma::cube& raw, const arma::cube& priors, const Rcpp::StringVector& prior_name, const Rcpp::Function& line_plot) {
  int nrow = int (raw.n_rows), ncol = int (raw.n_cols);
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (arma::is_finite(priors(i, j, 0))) {
        line_plot(Rcpp::_["data1"] = raw(arma::span(i), arma::span(j), arma::span::all), Rcpp::_["prior_name"] = prior_name, Rcpp::_["i"] = (i + 1), Rcpp::_["j"] = (j + 1));
      }
    }
  }
}

// ACF Plots.
void acf_plots(const arma::cube& raw, const arma::cube& priors, const Rcpp::StringVector& prior_name, const Rcpp::Function& acf_plot) {
  int nrow = int (raw.n_rows), ncol = int (raw.n_cols);
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (arma::is_finite(priors(i, j, 0))) {
        acf_plot(Rcpp::_["data1"] = raw(arma::span(i), arma::span(j), arma::span::all), Rcpp::_["prior_name"] = prior_name, Rcpp::_["i"] = (i + 1), Rcpp::_["j"] = (j + 1));
      }
    }
  }
}

// Estimate Historical Decompositions and Process Raw Results.
arma::cube hd_estimates(const arma::cube& A_chain, const arma::cube& B_chain, const arma::mat& y1, const arma::mat& x1, const int pA_ncol, const int nlags, const int nsli, const double ci) {
  int nrow = int (y1.n_rows), ncol = (pA_ncol * pA_ncol);
  
  arma::cube HD_chain(nrow, ncol, nsli);
  std::fill(HD_chain.begin(), HD_chain.end(), Rcpp::NumericVector::get_na());
  
  arma::mat HD_draw((nrow + nlags), (pA_ncol * pA_ncol), arma::fill::zeros);
  arma::mat u1(nrow, pA_ncol), H_draw(pA_ncol, pA_ncol), Phi_draw(((pA_ncol * nlags) + 1), pA_ncol);
  std::fill(u1.begin(), u1.end(), Rcpp::NumericVector::get_na());
  std::fill(H_draw.begin(), H_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(Phi_draw.begin(), Phi_draw.end(), Rcpp::NumericVector::get_na());
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    //historical decomposition (HD_chain, U_chain, H_chain, Phi_chain)
    u1 = (y1 * A_chain.slice(t)) - (x1 * B_chain.slice(t));
    H_draw = arma::inv(A_chain.slice(t));
    Phi_draw = B_chain.slice(t) * arma::inv(A_chain.slice(t));
    for (int i = 0; i < pA_ncol; ++i) {
      HD_draw(arma::span(nlags, (nrow + nlags - 1)), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) = u1.col(i) * H_draw.row(i);
      for (int j = (nlags + 1); j < (nrow + nlags); ++j) {
        for (int k = 0; k < nlags; ++k) {
          HD_draw(arma::span(j), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) += HD_draw(arma::span(j - k - 1), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) * Phi_draw(arma::span((pA_ncol * k), ((pA_ncol * (k + 1)) - 1)), arma::span(0, (pA_ncol - 1)));
        }
      }
    }
    HD_chain.slice(t) = HD_draw(arma::span(nlags, (nrow + nlags - 1)), arma::span(0, ((pA_ncol * pA_ncol) - 1)));
    
  }
  
  //process HD_chain
  arma::mat temp(nrow, nsli), temp1(nrow, 1);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());
  
  arma::cube HD(nrow, ncol, 3);
  std::fill(HD.begin(), HD.end(), Rcpp::NumericVector::get_na());
  
  double total = double (nsli);
  int lb = int (std::ceil(total * (1.0 - ci))) - 1;
  int ub = int (std::ceil(total * ci)) - 1;

  for (int j = 0; j < ncol; ++j) {
    
    temp = HD_chain(arma::span::all,arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    HD.slice(0).col(j) = temp1.col(lb);
    HD.slice(1).col(j) = arma::median(temp1, 1);
    HD.slice(2).col(j) = temp1.col(ub);
  
  }
  
  return HD;
}

// Estimate Impulse Responses and Process Raw Results.
arma::cube irf_estimates(const arma::cube& A_chain, const arma::cube& B_chain, const int pA_ncol, const int nlags, const int nsli, const int h1_irf, const bool acc_irf, const double ci) {
  int nrow = (h1_irf + 1), ncol = (pA_ncol * pA_ncol);
  
  arma::cube IRF_chain(nrow, ncol, nsli);
  std::fill(IRF_chain.begin(), IRF_chain.end(), Rcpp::NumericVector::get_na());
  
  arma::mat H_draw(pA_ncol, pA_ncol), IRF_draw((h1_irf + nlags), (pA_ncol * pA_ncol)), Phi_draw(((pA_ncol * nlags) + 1), pA_ncol);
  std::fill(H_draw.begin(), H_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(IRF_draw.begin(), IRF_draw.end(), Rcpp::NumericVector::get_na());
  std::fill(Phi_draw.begin(), Phi_draw.end(), Rcpp::NumericVector::get_na());
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    //impulse responses (H_chain, IRF_chain, Phi_chain)
    Phi_draw = B_chain.slice(t) * arma::inv(A_chain.slice(t));
    H_draw = arma::inv(A_chain.slice(t));
    IRF_draw.fill(0.0);
    for (int i = 0; i < pA_ncol; ++i) {
      IRF_draw(arma::span(nlags - 1), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) = H_draw(arma::span(i), arma::span(0, (pA_ncol - 1)));
      for (int j = 0; j < h1_irf; ++j) {
        for (int k = 0; k < nlags; ++k) {
          IRF_draw(arma::span(j + nlags), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) = IRF_draw(arma::span(j + nlags), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) + IRF_draw(arma::span(j + nlags - k - 1), arma::span((pA_ncol * (i)), ((pA_ncol * (i + 1)) - 1))) * Phi_draw(arma::span((pA_ncol * k), ((pA_ncol * (k + 1)) - 1)), arma::span(0, (pA_ncol - 1)));
        }
      }
      if (acc_irf == true) {
        for (int j = 0; j < h1_irf; ++j) {
          IRF_draw(arma::span(j + nlags), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) = IRF_draw(arma::span(j + nlags), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1))) + IRF_draw(arma::span(j + nlags - 1), arma::span((pA_ncol * i), ((pA_ncol * (i + 1)) - 1)));
        }
      }
    }
    IRF_chain.slice(t) = IRF_draw(arma::span((nlags - 1), (h1_irf + nlags - 1)), arma::span(0, ((pA_ncol * pA_ncol) - 1)));
    
  }
  
  //process IRF_chain
  arma::mat temp(nrow, nsli), temp1(nrow, 1);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());
  
  arma::cube IRF(nrow, ncol, 3);
  std::fill(IRF.begin(), IRF.end(), Rcpp::NumericVector::get_na());
  
  double total = double (nsli);
  int lb = int (std::ceil(total * (1.0 - ci))) - 1;
  int ub = int (std::ceil(total * ci)) - 1;

  for (int j = 0; j < ncol; ++j) {
    
    temp = IRF_chain(arma::span::all, arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    IRF.slice(0).col(j) = temp1.col(lb);
    IRF.slice(1).col(j) = arma::median(temp1, 1);
    IRF.slice(2).col(j) = temp1.col(ub);
    
  }
  
  return IRF;
}

// Estimate Phi and Process Raw Results.
arma::cube phi_estimates(const arma::cube& A_chain, const arma::cube& B_chain, const int pA_ncol, const int nlags, const int nsli, const double ci) {
  int nrow = ((pA_ncol * nlags) + 1), ncol = pA_ncol;
  
  arma::cube Phi_chain(nrow, ncol, nsli);
  std::fill(Phi_chain.begin(), Phi_chain.end(), Rcpp::NumericVector::get_na());
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    Phi_chain.slice(t) = B_chain.slice(t) * arma::inv(A_chain.slice(t));
    
  }
  
  //process Phi_chain
  arma::mat temp(nrow, nsli), temp1(nrow, 1);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  std::fill(temp1.begin(), temp1.end(), Rcpp::NumericVector::get_na());
  
  arma::cube Phi(nrow, ncol, 3);
  std::fill(Phi.begin(), Phi.end(), Rcpp::NumericVector::get_na());
  
  double total = double (nsli);
  int lb = int (std::ceil(total * (1.0 - ci))) - 1;
  int ub = int (std::ceil(total * ci)) - 1;

  for (int j = 0; j < ncol; ++j) {
    
    temp = Phi_chain(arma::span::all, arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    Phi.slice(0).col(j) = temp1.col(lb);
    Phi.slice(1).col(j) = arma::median(temp1, 1);
    Phi.slice(2).col(j) = temp1.col(ub);
    
  }
  
  return Phi;
}

// Estimate H and Process Raw Results.
Rcpp::List h_estimates(const arma::cube& A_chain, const int pA_ncol, const int nsli, const arma::cube& pH, const Rcpp::Function& line_plot, const Rcpp::Function& acf_plot, const double ci) {
  arma::cube H_chain(pA_ncol, pA_ncol, nsli);
  std::fill(H_chain.begin(), H_chain.end(), Rcpp::NumericVector::get_na());
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    H_chain.slice(t) = arma::inv(A_chain.slice(t));
    
  }
  
  acf_plots(H_chain, pH, "pH", acf_plot);
  line_plots(H_chain, pH, "pH", line_plot);
  
  return Rcpp::List::create(Rcpp::Named("H_den") = den_function(H_chain, pH),
                            Rcpp::Named("H") = results_function(H_chain, ci));
}

// Estimate det(A) and Process Raw Results.
Rcpp::List deta_estimates(const arma::cube& A_chain, const int nsli, const arma::cube& pdetA, const Rcpp::Function& line_plot, const Rcpp::Function& acf_plot, const double ci) {
  arma::cube detA_chain(1, 1, nsli);
  std::fill(detA_chain.begin(), detA_chain.end(), Rcpp::NumericVector::get_na());
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    detA_chain.slice(t) = arma::det(A_chain.slice(t));
    
  }
  
  acf_plots(detA_chain, pdetA, "pdetA", acf_plot);
  line_plots(detA_chain, pdetA, "pdetA", line_plot);
  
  return Rcpp::List::create(Rcpp::Named("detA_den") = den_function(detA_chain, pdetA),
                            Rcpp::Named("detA") = results_function(detA_chain, ci));
}

// Run the BH_SBVAR Model.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List MAIN(const arma::mat& y1, const arma::mat& x1, const arma::mat& omega, const arma::mat& somega, const int nlags, const arma::cube& pA, const arma::cube& pdetA, const arma::cube& pH, const arma::mat& pP, const arma::mat& pP_sig, const arma::cube& pR_sig, const arma::mat& kappa1, const arma::mat& A_start, const int itr, const int burn, const int thin, const arma::mat& scale1, const int h1_irf, const bool acc_irf, const double ci, const Rcpp::StringVector& varnames, const Rcpp::Function& line_plot, const Rcpp::Function& acf_plot) {
  int pA_ncol = int (pA.n_cols), B_nrow = ((int (y1.n_cols) * nlags) + 1);
  double totals = double (itr - burn), totalt = double (thin), y1_nrow = double (y1.n_rows);
  int nsli = int (std::floor((totals / totalt)));
  //start Metropolis-Hastings algorithm
  arma::field<arma::cube> list1 = MH(y1, x1, nlags, omega, somega, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, A_start, itr, burn, scale1);
  double accept_rate = list1(0)(0, 0, 0);
  arma::cube A_chain = list1(1);
  arma::cube B_chain = list1(2);
  arma::cube zeta_chain = list1(3);
  
  if (thin > 1) {
    A_chain = thin_function(A_chain, thin);
    B_chain = thin_function(B_chain, thin);
    zeta_chain = thin_function(zeta_chain, thin);
  }
  
  arma::mat taustar(pA_ncol, pA_ncol), kappastar(pA_ncol, pA_ncol), M(B_nrow, pA_ncol);
  std::fill(taustar.begin(), taustar.end(), Rcpp::NumericVector::get_na());
  std::fill(kappastar.begin(), kappastar.end(), Rcpp::NumericVector::get_na());
  std::fill(M.begin(), M.end(), Rcpp::NumericVector::get_na());
  
  arma::mat Dinv_draw(pA_ncol, pA_ncol, arma::fill::zeros);
  
  //estimate matrix B
  for (int t = 0; t < nsli; ++t) {
    
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    taustar = (arma::diagmat(kappa1) * arma::diagmat(A_chain.slice(t).t() * somega * A_chain.slice(t))) + arma::diagmat((1.0 / 2.0) * zeta_chain.slice(t));
    kappastar = arma::diagmat(kappa1 + ((y1_nrow / 2.0) * arma::ones(1, pA_ncol)));
    
    for (int i = 0; i < pA_ncol; ++i) {
      Dinv_draw(i, i) = R::rgamma((double (kappastar(i, i))), (double (1.0 / (double (taustar(i, i))))));
    }
    
    for (int i = 0; i < pA_ncol; ++i) {
      M = arma::chol(arma::inv((x1.t() * x1) + pP_sig + pR_sig.slice(i))).t() * arma::randn(B_nrow,pA_ncol) * arma::inv(arma::sqrt(Dinv_draw));
      B_chain(arma::span::all, arma::span(i), arma::span(t)) += M.col(i);
    }
    
  }
  
  //construct pR
  arma::mat pR(((pA_ncol * nlags) + 1), pA_ncol, arma::fill::zeros);
  for (int i = 0; i < pA_ncol; ++i) {
    for (int j = 0; j < pA_ncol; ++j) {
      if (arma::is_finite(pA(j, i, 6))) {
        pR(j, i) = pA(j, i, 6);
      }
    }
  }
  
  acf_plots(A_chain, pA, "pA", acf_plot);
  line_plots(A_chain, pA, "pA", line_plot);
  
  Rcpp::List list2 = deta_estimates(A_chain, nsli, pdetA, line_plot, acf_plot, ci);
  Rcpp::List list3 = h_estimates(A_chain, pA_ncol, nsli, pH, line_plot, acf_plot, ci);
  
  Rcpp::List list4;
  
  list4["accept_rate"] = accept_rate;
  
  list4["y"] = y1;
  list4["x"] = x1;
  
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
  list4["Phi"] = phi_estimates(A_chain, B_chain, pA_ncol, nlags, nsli, ci);
  
  list4["HD"] = hd_estimates(A_chain, B_chain, y1, x1, pA_ncol, nlags, nsli, ci);
  list4["IRF"] = irf_estimates(A_chain, B_chain, pA_ncol, nlags, nsli, h1_irf, acc_irf, ci);
  
  list4["A_den"] = den_function(A_chain, pA);
  list4["detA_den"] = list2["detA_den"];
  list4["H_den"] = list3["H_den"];
  
  return list4;
}


