

#include "RcppArmadillo.h"
//[[Rcpp::depends(RcppArmadillo)]]





// Non-Central T-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_nonc_t(double x1, double c1, double sigma1, double nu, double lam1) {
  double den = R::dnt(((x1 - c1) / sigma1), nu, lam1, false) / sigma1;
  return den;
}

// Negative T-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_t_n(double x1, double c1, double sigma1, double nu) {
  double den = R::dt(((x1 - c1) / sigma1), nu, false) / (sigma1 * (R::pt((((-1) * c1) / sigma1), nu, true, false)));
  return den;
}

// Positive T-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_t_p(double x1, double c1, double sigma1, double nu) {
  double den = R::dt(((x1 - c1) / sigma1), nu, false) / (sigma1 * (1 - R::pt((((-1) * c1) / sigma1), nu, true, false)));
  return den;
}

// Indifferent T-Distribution.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double prior_t(double x1, double c1, double sigma1, double nu) {
  double den = R::dt(((x1 - c1) / sigma1), nu, false) / sigma1;
  return den;
}

// Sum the Log of Prior Densities.
double sum_log_prior_densities(arma::mat A_test, arma::cube pA, arma::cube pdetA, arma::cube pH) {
  //dimension parameters
  int nrow = pA.n_rows, ncol = pA.n_cols;
  //construct H_test
  arma::mat H_test = arma::inv(A_test);
  //construct sum_log_priors and detA_test
  double sum_log_priors = 0, detA_test = arma::det(A_test);
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

// Likelihood Function for the Structural Parameters in Matrix A
double likelihood_function(arma::mat A_test, arma::mat kappa1, arma::mat y1, arma::mat omega, arma::mat zeta, arma::mat somega) {
  int kappa_ncol = kappa1.n_cols;
  double ynrow = y1.n_rows;
  double lik_A_numerator = (ynrow / 2.0) * std::log(arma::det(A_test.t() * omega * A_test));
  arma::mat tau = arma::diagmat(kappa1) * arma::diagmat(A_test.t() * somega * A_test);
  for (int i = 0; i < kappa_ncol; ++i) {
    if (kappa1(0, i) > 0) {
      lik_A_numerator += kappa1(0, i) * std::log(tau(i, i));
    }
  }
  arma::mat lik_A_denomenator = (kappa1 + ((ynrow / 2.0) * arma::ones(1, kappa_ncol))) * arma::log(arma::diagvec((2.0 / ynrow) * (tau + (arma::diagmat(zeta) / 2.0))));
  double lik_A = lik_A_numerator - lik_A_denomenator(0, 0);
  return lik_A;
}

// Construct zeta and B
Rcpp::List zB_function(arma::mat y1, arma::mat x1, int nlags, arma::mat A_test, arma::cube pA, arma::mat pP, arma::mat pP_sig, arma::cube pR_sig) {
  //dimension parameters
  int pA_nrow = pA.n_rows, pA_ncol = pA.n_cols;
  //create storage matrices and cubes
  arma::mat B(((pA_ncol * nlags) + 1), pA_ncol);
  std::fill(B.begin(), B.end(), Rcpp::NumericVector::get_na());
  arma::mat zeta(pA_nrow, pA_ncol, arma::fill::zeros);
  arma::cube pR(((pA_ncol * nlags) + 1), pA_ncol, pA_ncol, arma::fill::zeros);
  arma::mat B_temp = B;
  arma::mat zeta_temp(pA_nrow, pA_ncol, arma::fill::zeros);
  //construce pR
  for (int i = 0; i < pA_ncol; ++i) {
    for (int j = 0; j < pA_nrow; ++j) {
      if (!R_IsNA(pA(j, i, 6))) {
        pR(j, i, i) = A_test(j, i);
      }
    }
    //compute B and zeta estimates
    B_temp = arma::solve((x1.t() * x1 + pP_sig + pR_sig.slice(i)), ((x1.t() * y1 + pP_sig * pP + pR_sig.slice(i) * pR.slice(i)) * A_test));
    zeta_temp = A_test.t() * ((y1.t() * y1 + pP.t() * pP_sig * pP + pR.slice(i).t() * pR_sig.slice(i) * pR.slice(i)) - ((y1.t() * x1 + pP.t() * pP_sig + pR.slice(i).t() * pR_sig.slice(i)) * arma::inv(x1.t() * x1 + pP_sig + pR_sig.slice(i)) * (y1.t() * x1 + pP.t() * pP_sig + pR.slice(i).t() * pR_sig.slice(i)).t())) * A_test;
    //store estimates
    B.col(i) = B_temp.col(i);
    zeta(i, i) = zeta_temp(i, i);
    //reset temp matrices
    B_temp.reset();
    zeta_temp.reset();
  }
  //load a list
  Rcpp::List list1;
  list1["B"] = B;
  list1["zeta"] = zeta;
  return list1;
}

// Posterior Function
double post_A_function(arma::mat A_test, arma::cube pA, arma::cube pdetA, arma::cube pH, arma::mat kappa1, arma::mat y1, arma::mat omega, arma::mat zeta, arma::mat somega) {
  double priors = sum_log_prior_densities(A_test, pA, pdetA, pH);
  double likelihood = likelihood_function(A_test, kappa1, y1, omega, zeta, somega);
  double posterior = priors + likelihood;
  return posterior;
}

// Posterior Function used to Determine Starting Values for the Matrix A
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
double post_A_optim(Rcpp::NumericVector par, arma::cube pA, arma::cube pdetA, arma::cube pH, arma::mat pP, arma::mat pP_sig, arma::cube pR_sig, arma::mat kappa1, arma::mat y1, arma::mat x1, arma::mat omega, arma::mat somega, int nlags) {
  
  //check for interruptions
  Rcpp::checkUserInterrupt();
  
  //dimension parameters
  int pA_nrow = pA.n_rows, pA_ncol = pA.n_cols;
  //construct Anew
  arma::mat Anew(pA_nrow, pA_ncol);
  std::fill(Anew.begin(), Anew.end(), Rcpp::NumericVector::get_na());
  //construct n
  int n = 0;
  //fill in Anew
  for (int i = 0; i < pA_ncol; ++i) {
    for (int j = 0; j < pA_nrow; ++j) {
      if (!R_IsNA(pA(j, i, 0))) {
        Anew(j, i) = par[n];
        n += 1;
      } else {
        Anew(j, i) = pA(j, i, 2);
      }
    }
  }
  //construct zeta
  Rcpp::List list1 = zB_function(y1, x1, nlags, Anew, pA, pP, pP_sig, pR_sig);
  arma::mat zeta = list1["zeta"];
  //construct priors likelihood, and posterior
  double priors, likelihood, posterior;
  //compute posterior given sign restrictions
  priors = sum_log_prior_densities(Anew, pA, pdetA, pH);
  likelihood = likelihood_function(Anew, kappa1, y1, omega, zeta, somega);
  posterior = -(priors + likelihood);

  return posterior;
}

// Genterate Proposal Values for the Elements of Matrix A
Rcpp::List proposal(arma::mat A_old, arma::cube pA, arma::cube pdetA, arma::cube pH, arma::mat scale1) {
  //dimension parameters
  int nrow = pA.n_rows, ncol = pA.n_cols;
  //identify the number of priors for A
  int nA = arma::size(pA.slice(0).elem( arma::find_finite(pA.slice(0)) ), 0);
  //construct Adraw, Aold, Adist, and Asign
  arma::mat Adraw(nA, 1);
  std::fill(Adraw.begin(), Adraw.end(), Rcpp::NumericVector::get_na());
  arma::mat Aold(nA, 1);
  std::fill(Aold.begin(), Aold.end(), Rcpp::NumericVector::get_na());
  arma::mat Adist(nA, 1);
  std::fill(Adist.begin(), Adist.end(), Rcpp::NumericVector::get_na());
  arma::mat Asign(nA, 1);
  std::fill(Asign.begin(), Asign.end(), Rcpp::NumericVector::get_na());
  //fill in Adraw, Aold, Adist, and Asign
  int a = 0;
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < nrow; ++j) {
      if (arma::is_finite(pA(j, i, 0))) {
        Adraw(a, 0) = R::rnorm(0, 1) / std::sqrt(0.5 * (std::pow(R::rnorm(0, 1), 2) + std::pow(R::rnorm(0, 1), 2)));
        Aold(a, 0) = A_old(j, i);
        Adist(a, 0) = pA(j, i, 0);
        Asign(a, 0) = pA(j, i, 1);
        a += 1;
      }
    }
  }
  //construct Anew
  arma::mat Anew = Aold + (scale1 * Adraw);
  //test sign restrictions for Anew if there are any
  int sign_test = 0;
  for (int i = 0; i < nA; ++i) {
    if (Adist(i, 0) == 0) { // if Adist(i, 0) == 0, symmetric t-distribution
      if ((Asign(i, 0) == 1) & (Anew(i, 0) < 0)) { // if Asign(i, 0) == 1, positive sign restriction
        sign_test += 1;
      }
      if ((Asign(i, 0) == (-1)) & (Anew(i, 0) > 0)) { // if Asign(i, 0) == 1, positive sign restriction
        sign_test += 1;
      }
    }
  }
  //construct A_test
  arma::mat A_test(nrow, ncol);
  //fill in A_test
  a = 0;
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < nrow; ++j) {
      if (arma::is_finite(pA(j, i, 0))) {
        A_test(j, i) = Anew(a, 0);
        a += 1;
      } else {
        A_test(j, i) = pA(j, i, 2);
      }
    }
  }
  //construct H_test
  arma::mat H_test = arma::inv(A_test);
  //identify the number of priors for H
  int nH = 0;
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < nrow; ++j) {
      if (arma::is_finite(pH(j, i, 0))) {
        nH += 1;
      }
    }
  }
  if (nH > 0) {
    //construct Hnew, Hdist, and Hsign
    arma::mat Hnew(nH, 1);
    std::fill(Hnew.begin(), Hnew.end(), Rcpp::NumericVector::get_na());
    arma::mat Hdist(nH, 1);
    std::fill(Hdist.begin(), Hdist.end(), Rcpp::NumericVector::get_na());
    arma::mat Hsign(nH, 1);
    std::fill(Hsign.begin(), Hsign.end(), Rcpp::NumericVector::get_na());
    //fill in Hnew Hdist, and Hsign
    int h = 0;
    for (int i = 0; i < ncol; ++i) {
      for (int j = 0; j < nrow; ++j) {
        if (arma::is_finite(pH(j, i, 0))) {
          Hnew(h, 0) = H_test(j, i);
          Hdist(h, 0) = pH(j, i, 0);
          Hsign(h, 0) = pH(j, i, 1);
          h += 1;
        }
      }
    }
    //test sign restrictions for Hnew if there are any
    for (int i = 0; i < nH; ++i) {
      if (Hdist(i, 0) == 0) { // if Hdist(i, 0) == 0, symmetric t-distribution
        if ((Hsign(i, 0) == 1) & (Hnew(i, 0) < 0)) { // if Hsign(i, 0) == 1, positive sign restriction
          sign_test += 1;
        }
        if ((Hsign(i, 0) == (-1)) & (Hnew(i, 0) > 0)) { // if Hsign(i, 0) == 1, positive sign restriction
          sign_test += 1;
        }
      }
    }
  }
  //construct detA_test
  double detA_test = arma::det(A_test);
  //test sign restrictions for detA_test if there are any
  if (pdetA(0, 0, 0) == 0) { // if pdetA(0, 0, 0) == 0, symmetric t-distribution
    if ((pdetA(0, 0, 1) == 1) & (detA_test < 0)) { // if pdetA(0, 0, 1) == 1, positive sign restriction
      sign_test += 1;
    }
    if ((pdetA(0, 0, 1) == (-1)) & (detA_test > 0)) { // if pdetA(0, 0, 1) == 1, positive sign restriction
      sign_test += 1;
    }
  }
  //load a list
  Rcpp::List list1;
  list1["A_test"] = A_test;
  list1["sign_test"] = sign_test;
  return list1;
}

// Test Proposal Values for the Elements of Matrix A
arma::mat proposal_function(arma::mat A_old, arma::cube pA, arma::cube pdetA, arma::cube pH, arma::mat scale1) {
  //int nrow = pA.n_rows, ncol = pA.n_cols;
  Rcpp::List test = proposal(A_old, pA, pdetA, pH, scale1);
  arma::mat A_test = test["A_test"];
  int sign_test = test["sign_test"];
  while (sign_test != 0) {
    test = proposal(A_old, pA, pdetA, pH, scale1);
    A_test = Rcpp::as<arma::mat>(test["A_test"]);
    sign_test = test["sign_test"];
  }
  return A_test;
}

// Start Random-Walk Metropolis-Hastigs Algorithm
Rcpp::List MH(arma::mat y1, arma::mat x1, int nlags, arma::mat omega, arma::mat somega, arma::cube pA, arma::cube pdetA, arma::cube pH, arma::mat pP, arma::mat pP_sig, arma::cube pR_sig, arma::mat kappa1, arma::mat A_old, int itr, int burn, arma::mat scale1) {
  //dimension parameters
  int pA_nrow = pA.n_rows, pA_ncol = pA.n_cols;
  //construct cubes to store estimates
  arma::cube A_chain(pA_nrow, pA_ncol, (itr - burn));
  std::fill(A_chain.begin(), A_chain.end(), Rcpp::NumericVector::get_na());
  arma::cube B_chain(((pA_ncol * nlags) + 1), pA_ncol, (itr - burn));
  std::fill(B_chain.begin(), B_chain.end(), Rcpp::NumericVector::get_na());
  arma::cube zeta_chain(pA_nrow, pA_ncol, (itr - burn));
  std::fill(zeta_chain.begin(), zeta_chain.end(), Rcpp::NumericVector::get_na());
  //construct B and zeta from starting values for A
  Rcpp::List list1 = zB_function(y1, x1, nlags, A_old, pA, pP, pP_sig, pR_sig);
  arma::mat B_old = list1["B"];
  arma::mat zeta_old = list1["zeta"];
  //construct matrices to test proposals for A
  arma::mat A_test, B_test, zeta_test;
  //construct parameters used in Metropolis-Hastings algorithm
  double post_A_old = post_A_function(A_old, pA, pdetA, pH, kappa1, y1, omega, zeta_old, somega), post_A_test = 0.0, accept = 0.0, naccept = 0.0, accept_rate = 0.0, ru = 0.0;
  //start Metropolis-Hastings algorithm
  for (int c = 0; c < itr; ++c) {
    //check for interruptions
    if (c % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    //random draw proposal for A, B, and zeta
    A_test = proposal_function(A_old, pA, pdetA, pH, scale1);
    list1 = zB_function(y1, x1, nlags, A_test, pA, pP, pP_sig, pR_sig);
    B_test = Rcpp::as<arma::mat>(list1["B"]);
    zeta_test = Rcpp::as<arma::mat>(list1["zeta"]);
    //posterior values for A_test
    post_A_test = post_A_function(A_test, pA, pdetA, pH, kappa1, y1, omega, zeta_test, somega);
    //acceptance value to be tested agains a threshold to determing if the proposals for A will be kept
    accept = std::exp(post_A_test - post_A_old);
    //threshold to determine whether to accept proposal for A
    ru = R::runif(0, 1);
    //determing if the proposals for A will be kept
    if (ru <= accept) {
      //proposals for A, B, and zeta become old proposals as the chain moves forward
      A_old = A_test;
      B_old = B_test;
      zeta_old = zeta_test;
      post_A_old = post_A_test;
      if (c >= burn) {
        //store new proposals for A, B, and zeta
        A_chain.slice((c - burn)) = A_test;
        B_chain.slice((c - burn)) = B_test;
        zeta_chain.slice((c - burn)) = zeta_test;
        naccept += 1;
      }
    } else {
      if (c >= burn) {
        //store old proposals for A, B, and zeta
        A_chain.slice((c - burn)) = A_old;
        B_chain.slice((c - burn)) = B_old;
        zeta_chain.slice((c - burn)) = zeta_old;
      }
    }
  }
  //compute acceptance rate
  double total = (itr - burn);
  accept_rate = (naccept / total);
  //load a list
  Rcpp::List list2;
  list2["accept_rate"] = accept_rate;
  list2["A_chain"] = A_chain;
  list2["B_chain"] = B_chain;
  list2["zeta_chain"] = zeta_chain;
  return list2;
}

// Thin the Markov Chains
arma::cube thin_function(arma::cube chain, int thin) {
  //dimension parameters
  int nrow = chain.n_rows, ncol = chain.n_cols;
  double totals = chain.n_slices, totalt = thin;
  int nsli = std::floor((totals / totalt));
  //construct cube to stor estimates
  arma::cube chain1(nrow, ncol, nsli);
  std::fill(chain1.begin(), chain1.end(), Rcpp::NumericVector::get_na());
  //store a fraction of the estimates
  for (int i = 0; i < nsli; ++i) {
    chain1.slice(i) = chain.slice((thin * (i + 1) - 1));
  }
  return chain1;
}

// Process raw results into lower bound, median, and upper bound.
arma::cube results_function(arma::cube raw, double ci) {
  int nrow = raw.n_rows, ncol = raw.n_cols, nsli = raw.n_slices;
  
  arma::mat temp(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  
  arma::mat temp1(nrow, 1);
  
  arma::cube results1(nrow, ncol, 3);
  std::fill(results1.begin(), results1.end(), Rcpp::NumericVector::get_na());
  
  double total = nsli;
  int lb = std::ceil(total * (1.0 - ci)) - 1;
  int ub = std::ceil(total * ci) - 1;

  for (int j = 0; j < ncol; ++j) {
    
    temp = raw(arma::span::all, arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    results1(arma::span::all, arma::span(j), arma::span(0)) = temp1(arma::span::all, arma::span(lb));
    results1(arma::span::all, arma::span(j), arma::span(1)) = median(temp1, 1);
    results1(arma::span::all, arma::span(j), arma::span(2)) = temp1(arma::span::all, arma::span(ub));
  
  }
  
  return results1;
  
}

// Estimate density coordinates
Rcpp::List den_function(arma::cube raw, arma::cube priors) {
  
  // obtain environment containing function
  Rcpp::Environment stats_r("package:stats");
  // make function callable from c++
  Rcpp::Function density_r = stats_r["density"];
  
  int nrow = raw.n_rows, ncol = raw.n_cols, nsli = raw.n_slices;
  
  double t1;
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
      if (!R_IsNA(priors(i, j, 0))) {
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
  
  Rcpp::List list2;
  list2["hori"] = hori;
  list2["vert"] = vert;
  
  return list2;
}

// Line plots
void line_plots(arma::cube raw, arma::cube priors, Rcpp::StringVector prior_name, Rcpp::Function line_plot) {
  int nrow = raw.n_rows, ncol = raw.n_cols;
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (!R_IsNA(priors(i, j, 0))) {
        line_plot(Rcpp::_["data1"] = raw(arma::span(i), arma::span(j), arma::span::all), Rcpp::_["prior_name"] = prior_name, Rcpp::_["i"] = (i + 1), Rcpp::_["j"] = (j + 1));
      }
    }
  }
}

// ACF plots
void acf_plots(arma::cube raw, arma::cube priors, Rcpp::StringVector prior_name, Rcpp::Function acf_plot) {
  int nrow = raw.n_rows, ncol = raw.n_cols;
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (!R_IsNA(priors(i, j, 0))) {
        acf_plot(Rcpp::_["data1"] = raw(arma::span(i), arma::span(j), arma::span::all), Rcpp::_["prior_name"] = prior_name, Rcpp::_["i"] = (i + 1), Rcpp::_["j"] = (j + 1));
      }
    }
  }
}

// Estimate lower bounds, medians, and upper bounds for the elements in HD
arma::cube hd_estimates(arma::cube A_chain, arma::cube B_chain, arma::mat y1, arma::mat x1, int pA_ncol, int nlags, int nsli, double ci) {
  int nrow = y1.n_rows, ncol = (pA_ncol * pA_ncol);
  //create cubes to store parameters
  arma::cube HD_chain(nrow, ncol, nsli);
  std::fill(HD_chain.begin(), HD_chain.end(), Rcpp::NumericVector::get_na());
  
  arma::mat HD_draw((nrow + nlags), (pA_ncol * pA_ncol), arma::fill::zeros), u1, H_draw, Phi_draw;
  
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
  // Process HD_chain
  arma::mat temp(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  
  arma::mat temp1(nrow, 1);
  
  arma::cube HD(nrow, ncol, 3);
  std::fill(HD.begin(), HD.end(), Rcpp::NumericVector::get_na());
  
  double total = nsli;
  int lb = std::ceil(total * (1.0 - ci)) - 1;
  int ub = std::ceil(total * ci) - 1;

  for (int j = 0; j < ncol; ++j) {
    
    temp = HD_chain(arma::span::all,arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    HD(arma::span::all, arma::span(j), arma::span(0)) = temp1(arma::span::all, arma::span(lb));
    HD(arma::span::all, arma::span(j), arma::span(1)) = median(temp1, 1);
    HD(arma::span::all, arma::span(j), arma::span(2)) = temp1(arma::span::all, arma::span(ub));
  
  }
  
  return HD;
}

// Estimate lower bounds, medians, and upper bounds for the elements in IRF
arma::cube irf_estimates(arma::cube A_chain, arma::cube B_chain, int pA_ncol, int nlags, int nsli, int h1_irf, bool acc_irf, double ci) {
  int nrow = (h1_irf + 1), ncol = (pA_ncol * pA_ncol);
  //create cubes to store parameters
  arma::cube IRF_chain(nrow, ncol, nsli);
  std::fill(IRF_chain.begin(), IRF_chain.end(), Rcpp::NumericVector::get_na());
  
  arma::mat H_draw, IRF_draw, Phi_draw;
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    //impulse responses (H_chain, IRF_chain, Phi_chain)
    Phi_draw = B_chain.slice(t) * arma::inv(A_chain.slice(t));
    H_draw = arma::inv(A_chain.slice(t));
    IRF_draw = arma::zeros((h1_irf + nlags), (pA_ncol * pA_ncol));
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
  // Process IRF_chain
  arma::mat temp(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  
  arma::mat temp1(nrow, 1);
  
  arma::cube IRF(nrow, ncol, 3);
  std::fill(IRF.begin(), IRF.end(), Rcpp::NumericVector::get_na());
  
  double total = nsli;
  int lb = std::ceil(total * (1.0 - ci)) - 1;
  int ub = std::ceil(total * ci) - 1;

  for (int j = 0; j < ncol; ++j) {
    
    temp = IRF_chain(arma::span::all, arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    IRF(arma::span::all, arma::span(j), arma::span(0)) = temp1(arma::span::all, arma::span(lb));
    IRF(arma::span::all, arma::span(j), arma::span(1)) = median(temp1, 1);
    IRF(arma::span::all, arma::span(j), arma::span(2)) = temp1(arma::span::all, arma::span(ub));
  
  }
  
  return IRF;
}

// Estimate lower bounds, medians, and upper bounds for the elements in Phi
arma::cube phi_estimates(arma::cube A_chain, arma::cube B_chain, int pA_ncol, int nlags, int nsli, double ci) {
  int nrow = ((pA_ncol * nlags) + 1), ncol = pA_ncol;
  //create cubes to store parameters
  arma::cube Phi_chain(nrow, ncol, nsli);
  std::fill(Phi_chain.begin(), Phi_chain.end(), Rcpp::NumericVector::get_na());
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    Phi_chain.slice(t) = B_chain.slice(t) * arma::inv(A_chain.slice(t));
    
  }
  // Process Phi_chain
  arma::mat temp(nrow, nsli);
  std::fill(temp.begin(), temp.end(), Rcpp::NumericVector::get_na());
  
  arma::mat temp1(nrow, 1);
  
  arma::cube Phi(nrow, ncol, 3);
  std::fill(Phi.begin(), Phi.end(), Rcpp::NumericVector::get_na());
  
  double total = nsli;
  int lb = std::ceil(total * (1.0 - ci)) - 1;
  int ub = std::ceil(total * ci) - 1;

  for (int j = 0; j < ncol; ++j) {
    
    temp = Phi_chain(arma::span::all, arma::span(j), arma::span::all);
    temp1 = arma::sort(temp, "ascend", 1);
    Phi(arma::span::all, arma::span(j), arma::span(0)) = temp1(arma::span::all, arma::span(lb));
    Phi(arma::span::all, arma::span(j), arma::span(1)) = median(temp1, 1);
    Phi(arma::span::all, arma::span(j), arma::span(2)) = temp1(arma::span::all, arma::span(ub));
  
  }
  
  return Phi;
}

// Estimate lower bounds, medians, and upper bounds for the elements in H
Rcpp::List h_estimates(arma::cube A_chain, int pA_ncol, int nsli, arma::cube pH, Rcpp::Function line_plot, Rcpp::Function acf_plot, double ci) {
  //create cubes to store parameters
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
  
  Rcpp::List list1;
  list1["H_den"] = den_function(H_chain, pH);
  list1["H"] = results_function(H_chain, ci);
  
  return list1;
}

// Estimate lower bound, median, and upper bound for det(A)
Rcpp::List deta_estimates(arma::cube A_chain, int nsli, arma::cube pdetA, Rcpp::Function line_plot, Rcpp::Function acf_plot, double ci) {
  //create cubes to store parameters
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
  
  Rcpp::List list1;
  list1["detA_den"] = den_function(detA_chain, pdetA);
  list1["detA"] = results_function(detA_chain, ci);
  
  return list1;
}

// Main Function that Runs the BH_SBVAR model.
//' @useDynLib BHSBVAR, .registration = TRUE
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List MAIN(arma::mat y1, arma::mat x1, arma::mat omega, arma::mat somega, int nlags, arma::cube pA, arma::cube pdetA, arma::cube pH, arma::mat pP, arma::mat pP_sig, arma::cube pR_sig, arma::mat kappa1, arma::mat A_start, int itr, int burn, int thin, arma::mat scale1, int h1_irf, bool acc_irf, double ci, Rcpp::StringVector varnames, Rcpp::Function line_plot, Rcpp::Function acf_plot) {
  //dimension parameters
  int pA_ncol = pA.n_cols, B_nrow = ((y1.n_cols * nlags) + 1);
  double totals = itr - burn, totalt = thin, y1_nrow = y1.n_rows;
  int nsli = std::floor((totals / totalt));
  //metropolis hastings algorithm
  Rcpp::List list1 = MH(y1, x1, nlags, omega, somega, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, A_start, itr, burn, scale1);
  double accept_rate = list1["accept_rate"];
  arma::cube A_chain = thin_function(list1["A_chain"], thin);
  arma::cube B_chain = thin_function(list1["B_chain"], thin);
  arma::cube zeta_chain = thin_function(list1["zeta_chain"], thin);
  
  //create matrices
  arma::mat taustar, kappastar, Dinv_draw(pA_ncol, pA_ncol, arma::fill::zeros), M;
  
  //shape1 and rate1 are required because doubles are needed for R::gamma function
  double shape1 = 0, rate1 = 0;
  //estimate matrix B, Phi, Impulse responses, and historical decomposition
  
  for (int t = 0; t < nsli; ++t) {
    //check for interruptions
    if (t % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    taustar = (arma::diagmat(kappa1) * arma::diagmat(A_chain.slice(t).t() * somega * A_chain.slice(t))) + arma::diagmat((1.0 / 2.0) * zeta_chain.slice(t));
    kappastar = arma::diagmat(kappa1 + ((y1_nrow / 2.0) * arma::ones(1, pA_ncol)));
    for (int i = 0; i < pA_ncol; ++i) {
      //shape1 and rate1 are required because doubles are needed for R::gamma function
      shape1 = kappastar(i, i);
      rate1 = taustar(i, i);
      Dinv_draw(i, i) = R::rgamma(shape1, (1.0 / rate1));
    }
    for (int i = 0; i < pA_ncol; ++i) {
      M = arma::chol(arma::inv((x1.t() * x1) + pP_sig + pR_sig.slice(i))).t() * arma::randn(B_nrow,pA_ncol) * arma::inv(arma::sqrt(Dinv_draw));
      B_chain(arma::span::all, arma::span(i), arma::span(t)) += M.col(i);
    }
    
  }
  
  //construce pR with pA
  arma::mat pR(((pA_ncol * nlags) + 1), pA_ncol, arma::fill::zeros);
  for (int i = 0; i < pA_ncol; ++i) {
    for (int j = 0; j < pA_ncol; ++j) {
      if (!R_IsNA(pA(j, i, 6))) {
        pR(j, i) = pA(j, i, 6);
      }
    }
  }
  
  acf_plots(A_chain, pA, "pA", acf_plot);
  line_plots(A_chain, pA, "pA", line_plot);
  
  Rcpp::List list2;
  list2 = deta_estimates(A_chain, nsli, pdetA, line_plot, acf_plot, ci);
  Rcpp::List list3;
  list3 = h_estimates(A_chain, pA_ncol, nsli, pH, line_plot, acf_plot, ci);
  
  //load a list
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


