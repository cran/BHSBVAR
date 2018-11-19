
# Create matrices containing dependent and independet variables.
#' @keywords internal
getXY <- function(data,nlags) {
  data <- data - (matrix(data = 1, nrow = nrow(data), ncol = ncol(data)) %*% diag(x = colMeans(x = data, na.rm = FALSE, dims = 1), names = TRUE))
  X <- matrix(data = NA_real_, nrow = (nrow(data) - nlags), ncol = (ncol(data) * nlags))
  for (k in 1:nlags) {
    X[, (ncol(data) * (k - 1) + 1):(ncol(data) * k)] <- data[(nlags - k + 1):(nrow(data) - k), ]
  }
  X <- cbind(X, 1)
  colnames(X) <- c(paste(rep(colnames(data), nlags), ".L", sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(data)), decreasing = FALSE), sep = ""), "cons")
  Y <- data[(nlags + 1):nrow(data), ]
  list1 <- list(X, Y)
  return(list1)
}

# Check arguments from the BH_SBVAR function that should be integers.
#' @keywords internal
check_integers <- function(list1) {
  #testing inputs that should be integers
  for (i in 1:length(list1)) {
    if (((class(list1[[i]]) != "numeric") & (class(list1[[i]]) != "integer")) || (!is.finite(list1[[i]]))) {
      return(paste(names(list1[i]), ": Must be finite 'numeric' or 'integer'.", sep = ""))
    }
    if ((list1[[i]] %% 1) != 0) {
      return(paste(names(list1[i]), ": Must be a whole number.", sep = ""))
    }
    if ((names(list1[i]) == "nlags") && (list1[[i]] <= 0)) {
      return(paste(names(list1[i]), ": Must be greater than 0.", sep = ""))
    }
    if ((names(list1[i]) == "itr") && (list1[[i]] < 100)) {
      return(paste(names(list1[i]), ": Must be greater than 100.", sep = ""))
    }
    if ((names(list1[i]) == "burn") && (list1[[i]] < 0)) {
      return(paste(names(list1[i]), ": Must be greater than or equal to 0.", sep = ""))
    }
    if ((names(list1[i]) == "thin") && (list1[[i]] <= 0)) {
      return(paste(names(list1[i]), ": Must be greater than 0.", sep = ""))
    }
    if ((names(list1[i]) == "h1_irf") && (list1[[i]] < 4)) {
      return(paste(names(list1[i]), ": Must be greater than or equal to 3.", sep = ""))
    }
  }
  return("pass")
}

# Check arguments from the BH_SBVAR function that should be doubles.
#' @keywords internal
check_doubles <- function(list1) {
  #testing inputs that could be doubles
  for (i in 1:length(list1)) {
    if (((class(list1[[i]]) != "numeric") & (class(list1[[i]]) != "integer")) || (!is.finite(list1[[i]]))) {
      return(paste(names(list1[i]), ": Must be finite 'numeric' or 'integer'.", sep = ""))
    }
    if ((names(list1[i]) == "ci") && ((list1[[i]] < 0.7) | (list1[[i]] > 1))) {
      return(paste(names(list1[i]), ": Must be greater than or equal to 0.7 and less than or equal to 1.", sep = ""))
    }
  }
  return("pass")
}

# Check arguments from the BH_SBVAR function that should be matrices.
#' @keywords internal
check_matrices <- function(list1, nlags) {
  #testing inputs that should be matrices
  for (i in 1:length(list1)) {
    if (!is.matrix(list1[[i]])) {
      return(paste(names(list1[i]), ": Must be a matrix.", sep = ""))
    }
    if (any(!is.finite(list1[[i]]))) {
      return(paste(names(list1[i]), ": Elements must be finite numeric values", sep = ""))
    }
    if ((names(list1[i]) == "y") && (nrow(list1[[i]]) <= ncol(list1[[i]]))) {
      return("y: The number of rows must be greater than the number of columns.")
    }
    if ((names(list1[i]) == "y") && (ncol(list1[[i]]) < 2)) {
      return(paste("y: The number of columns or endogenous variables must be greater than 1.", sep = ""))
    }
    if ((names(list1[i]) == "y") && (((ncol(list1[[i]]) * nlags) + 1) >= (nrow(list1[[i]])))) {
      return(paste("y: The number observations must be greater than ", ((ncol(list1[[i]]) * nlags) + 1),". Reduce the number of lags or increase the number of observations.", sep = ""))
    }
    if ((names(list1[i]) == "pP") && (nrow(list1[[i]]) != ((nlags * ncol(list1$y)) + 1))) {
      return(paste("pP: The number of rows must equal ", ((nlags * ncol(list1$y)) + 1), ".", sep = ""))
    }
    if ((names(list1[i]) == "pP") && (ncol(list1[[i]]) != ncol(list1$y))) {
      return(paste("pP: The number of columns must equal ", (ncol(list1$y)), ".", sep = ""))
    }
    if ((names(list1[i]) == "pP_sig") && (nrow(list1[[i]]) != ((nlags * ncol(list1$y)) + 1))) {
      return(paste("pP_sig: The number of rows must equal ", ((nlags * ncol(list1$y)) + 1), ".", sep = ""))
    }
    if ((names(list1[i]) == "pP_sig") && (ncol(list1[[i]]) != ((nlags * ncol(list1$y)) + 1))) {
      return(paste("pP_sig: The number of columns must equal ",((nlags * ncol(list1$y)) + 1), ".", sep = ""))
    }
    if ((names(list1[i]) == "pP_sig") && (any(list1[[i]] < 0))) {
      return(paste("pP_sig: Elements must be greater than or equal to 0.", sep = ""))
    }
    if ((names(list1[i]) == "pP_sig") && (!isSymmetric(list1[[i]]))) {
      return(paste("pP_sig: Must be symmetric.", sep = ""))
    }
    if ((names(list1[i]) == "kappa1") && (nrow(list1[[i]]) != 1)) {
      return(paste("kappa1: The number of rows must equal 1.", sep = ""))
    }
    if ((names(list1[i]) == "kappa1") && (ncol(list1[[i]]) != ncol(list1$y))) {
      return(paste("kappa1: The number of columns must equal ", ncol(list1$y), ".", sep = ""))
    }
    if ((names(list1[i]) == "kappa1") && (any(list1[[i]] < 0))) {
      return(paste("kappa1: Elements must be greater than or equal to 0.", sep = ""))
    }
  }
  return("pass")
}

# Check arguments from the BH_SBVAR function that should be arrays.
#' @keywords internal
check_arrays <- function(list1, y) {
  for (i in 1:length(list1)) {
    if (!is.array(list1[[i]])) {
      return(paste(names(list1[i]), ": Must be an array.", sep = ""))
    }
    if (!is.numeric(list1[[i]])) {
      return(paste(names(list1[i]), ": Should contain 'numeric' elements for arrays specifying prior distributions. Use 'NA_real_' for elements in arrays that contain all NAs.", sep = ""))
    }
    if ((names(list1[i]) == "pA") && (all(is.na(list1[[i]][,,1])))) {
      return(paste(names(list1[i]), "[,,1]: Should indicate at least one parameter to be estimated.", sep = ""))
    }
    if ((names(list1[i]) == "pA") && ((dim(list1[[i]])[1] != ncol(y)) | (dim(list1[[i]])[2] != ncol(y)) | (dim(list1[[i]])[3] != 8))) {
      return(paste(names(list1[i]), ": Should be an (", ncol(y), ", ", ncol(y), ", 8) array.", sep = ""))
    }
    if ((names(list1[i]) == "pdetA") && ((dim(list1[[i]])[1] != 1) | (dim(list1[[i]])[2] != 1) | (dim(list1[[i]])[3] != 6))) {
      return(paste(names(list1[i]), ": Should be an (1, 1, 6) array.", sep = ""))
    }
    if ((names(list1[i]) == "pH") && ((dim(list1[[i]])[1] != ncol(y)) | (dim(list1[[i]])[2] != ncol(y)) | (dim(list1[[i]])[3] != 6))) {
      return(paste(names(list1[i]), ": Should be an (", ncol(y), ", ", ncol(y), ", 6) array.", sep = ""))
    }
    for (j in 1:(dim(list1[[i]])[1])) {
      for (k in 1:(dim(list1[[i]])[2])) {
        if (is.na(list1[[i]][j, k, 1])) { #if distribution is not specified, no other parameters should be specified
          if ((names(list1[i]) == "pA") && (any(!is.na(list1[[i]][j, k, c(1:2, 4:8)])))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 1]: Indicates no prior distribution so sign, scale, degrees of freedom, skew, long-run restriction, and proposal scaling parameter (", names(list1[i]),"[", j, ", ", k, ", c(2,4:7)]) should all be NA.", sep = ""))
          }
          if ((names(list1[i]) == "pA") && (!is.finite(list1[[i]][j, k, 3]))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 1]: Indicates no prior distribution so position (", names(list1[i]), "[", j, ", ", k, ", 3]) should be some constant value.", sep = ""))
          }
          if ((names(list1[i]) != "pA") && (any(!is.na(list1[[i]][j, k, 1:6])))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 1]: Indicates no prior distribution so sign, position, scale, degrees of freedom, skew (", names(list1[i]), "[", j, ", ", k, ", 1:6]) should all be NA.", sep = ""))
          }
        } else if (list1[[i]][j,k,1] == 0) { #if distribution is 0 (symmetric t-distribution), parameters in slices 2:5 must be specified
          if ((!is.na(list1[[i]][j, k, 2])) && ((list1[[i]][j, k, 2] != 1) & (list1[[i]][j, k, 2] != (-1)))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 2]: Sign should be indicated with a NA, 1, or -1.", sep = ""))
          }
          if (!is.finite(list1[[i]][j, k, 3])) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 3]: Position should be indicated with a finite number.", sep = ""))
          }
          if ((!is.na(list1[[i]][j, k, 2])) && ((list1[[i]][j, k, 3]) != 0) && ((list1[[i]][j, k, 2]) != ((list1[[i]][j, k, 3]) / abs(list1[[i]][j, k, 3])))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 3]: Position should have the same sign as sign (", names(list1[i]), "[", j, ", ", k, ", 2]).", sep = ""))
          }
          if ((!is.finite(list1[[i]][j, k, 4])) || (list1[[i]][j, k, 4] <= 0)) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 4]: Scale should be indicated with a finite number greater than 0.", sep = ""))
          }
          if ((!is.finite(list1[[i]][j, k, 5])) || (list1[[i]][j, k, 5] <= 2)) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 5]: Degrees of freedom should be indicated with a finite number greater than 2.", sep = ""))
          }
          if (any(!is.na(list1[[i]][j, k, 6]))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 6]: Skew should be NA.", sep = ""))
          }
          if ((names(list1[i]) == "pA") && ((!is.na(list1[[i]][j, k, 7])) && ((!is.finite(list1[[i]][j, k, 7])) || (list1[[i]][j, k, 7] != 1)))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 7]: Long-run restriction should be indicated with an NA (no long-run restriction) or a 1 (long-run restriction).", sep = ""))
          }
          if ((names(list1[i]) == "pA") && ((is.na(list1[[i]][j, k, 8])) || (!is.finite(list1[[i]][j, k, 8])) || (list1[[i]][j, k, 8] < 0.1))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 8]: Proposal scaling parameter should be greater than or equal to 0.1.", sep = ""))
          }
        } else if (list1[[i]][j, k, 1] == 1) { #if distribution is 1 (non-central t-distribution), parameters in slices 2:6 must be specified
          if (!is.na(list1[[i]][j, k, 2])) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 2]: Sign should be NA.", sep = ""))
          }
          if (!is.finite(list1[[i]][j, k, 3])) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 3]: Position should be indicated with a finite number.", sep = ""))
          }
          if ((!is.finite(list1[[i]][j, k, 4])) || (list1[[i]][j, k, 4] <= 0)) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 4]: Scale should be indicated with a finite number greater than 0.", sep = ""))
          }
          if ((!is.finite(list1[[i]][j, k, 5])) || (list1[[i]][j, k, 5] <= 2)) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 5]: Degrees of freedom should be indicated with a finite number greater than 2.", sep = ""))
          }
          if (!is.finite(list1[[i]][j, k, 6])) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 6]: Skew should be indicated with a finite number.", sep = ""))
          }
          if (((list1[[i]][j, k, 6] == 0) & (list1[[i]][j, k, 3] != 0)) | ((list1[[i]][j, k, 6] != 0) & (list1[[i]][j, k, 3] == 0))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 6]: Skew should be zero if position (", names(list1[i]), "[", j, ", ", k, ", 3]) is zero.", sep = ""))
          }
          if ((list1[[i]][j, k, 6] != 0) && (list1[[i]][j, k, 3] != 0) && (((list1[[i]][j, k, 6]) / abs(list1[[i]][j, k, 6])) != ((list1[[i]][j, k, 3]) / abs(list1[[i]][j, k, 3])))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 6]: Skew should have the same sign as position (", names(list1[i]), "[", j, ", ", k, ", 3]).", sep = ""))
          }
          if ((names(list1[i]) == "pA") && ((!is.na(list1[[i]][j, k, 7])) && ((!is.finite(list1[[i]][j, k, 7])) || (list1[[i]][j, k, 7] != 1)))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 7]: Long-run restriction should be indicated with an NA (no long-run restriction) or a 1 (long-run restriction).", sep = ""))
          }
          if ((names(list1[i]) == "pA") && ((is.na(list1[[i]][j, k, 8])) || (!is.finite(list1[[i]][j, k, 8])) || (list1[[i]][j, k, 8] < 0.1))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 8]: Proposal scaling parameter should be greater than or equal to 0.1.", sep = ""))
          }
        } else {
          return(paste(names(list1[i]), "[", j, ", ", k, ", 1]: Distribution should be indicated with a NA (no prior), 0 (symetric t-distribution), or 1 (non-central t-distribution).", sep = ""))
        }
      }
    }
  }
  return("pass")
}

# Check arguments from the BH_SBVAR function
#' @keywords internal
arguments_check <- function(y, nlags, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, itr, burn, thin, acc_irf, h1_irf, ci) {
  test <- check_integers(list(nlags = nlags, itr = itr, burn = burn, thin = thin, h1_irf = h1_irf))
  if (test != "pass") {
    return(test)
  }
  test <- check_doubles(list(ci = ci))
  if (test != "pass") {
    return(test)
  }
  if ((!is.logical(acc_irf)) || (is.na(acc_irf))) {
    return(paste("acc_irf: Must be logical 'TRUE' or 'FALSE'.", sep = ""))
  }
  if (floor((itr - burn) / thin) < 5000) {
    return(paste("'floor((itr-burn)/thin)' must be greater than or equal to 5000.", sep = ""))
  }
  test <- check_matrices(list(y = y, pP = pP, pP_sig = pP_sig, kappa1 = kappa1), nlags)
  if (test != "pass") {
    return(test)
  }
  test <- check_arrays(list(pA = pA, pdetA = pdetA, pH = pH), y)
  if (test != "pass") {
    return(test)
  }
  # check pR_sig
  if (!is.array(pR_sig)) {
    return("pR_sig: Must be an array.")
  }
  if (any(!is.finite(pR_sig)) || (any(pR_sig < 0))) {
    return("pR_sig: Mulst contain finite values greater than or equal to 0.")
  }
  if ((dim(pR_sig)[1] != ((nlags * ncol(y)) + 1)) | (dim(pR_sig)[2] != ((nlags * ncol(y)) + 1)) | (dim(pR_sig)[3] != ncol(y))) {
    return(paste("pR_sig: Dimensions should be (", ((nlags * ncol(y)) + 1), ", ", ((nlags * ncol(y)) + 1), ", ", (ncol(y)), ").", sep = ""))
  }
  for (i in 1:ncol(y)) {
    if (any(is.finite(pA[, i, 7]))) {
      n <- which(is.finite(pA[, i, 7]))
      for (j in n) {
        if (pR_sig[j, j, i] <= 0) {
          return(paste("pR_sig: The value at pR_sig[", j, ", ", j, ", ", i, "] should be a finite value greater than 0 since pA[",j,", ", i,", "," 7] indicates a long-run restriction.", sep = ""))
        }
      }
      if (any(!is.finite(pA[, i, 7]))) {
        n <- which(!is.finite(pA[, i, 7]))
        for (j in n) {
          for (k in seq(from = j, to = (dim(pR_sig)[1] - 1), by = ncol(y))) {
            for (l in seq(from = j, to = (dim(pR_sig)[2] - 1), by = ncol(y))) {
              if (pR_sig[k, l, i] != 0) {
                return(paste("pR_sig: The vlue at pR_sig[", k, ", ", l, ", ", i, "] should be 0.", sep = ""))
              }
            }
          }
        }
      }
      if ((any(pR_sig[(dim(pR_sig)[1]),, i] != 0)) | (any(pR_sig[, (dim(pR_sig)[2]), i] != 0))) {
        return(paste("pR_sig: The values at pR_sig[", (dim(pR_sig)[1]), ", ,", i, "] and pR_sig[,",(dim(pR_sig)[2]), ", ", i, "] should be 0.", sep = ""))
      }
    } else {
      if (any(pR_sig[,, i] != 0)) {
        return(paste("pR_sig: Each element in pR_sig[,,", i,"] should be 0 since there were no long-run restrictions indicated for equation ", i, ".", sep = ""))
      }
    }
    if (!isSymmetric(pR_sig[,,i])) {
      return(paste("pR_sig[,,",i,"]: Must be symmetric.", sep = ""))
    }
  }
  return("pass")
}

#' Structural Bayesian Vector Autoregression
#' 
#' Runs a Structural Bayesian Vector Autoregression model with the method developed by Baumeister and Hamilton (2015,2017,and 2018).
#' @author Paul Richardson
#' @export
#' @import Rcpp
#' @name BH_SBVAR
#' @param y Matrix containing the endogenous variables.
#' @param nlags Integer specifying the lag order.
#' @param pA Array where each slice of the third dimension contains the prior distributions, sign restrictions, distribution positions, distribution scales, distribution degrees of freedom, distribution skew, long-run restriction scale parameters, and random-walk proposal scale parameters for the coefficient matrix 'A', respectively.
#' @param pdetA Array where each slice of the third dimension contains the prior distributions, sign restrictions, distribution positions, distribution scales, distribution degrees of freedom, and distribution skew parameters for the determinant of 'A', respectively (default=NULL). 'NULL' indicates no priors for the determinant of 'A'.
#' @param pH Array where each slice of the third dimension contains the prior distributions, sign restrictions, distribution positions, distribution scales, distribution degrees of freedom, distribution skew parameters for the inverse of 'A', respectively (default=NULL). 'NULL' indicates no priors for the inverse of 'A'.
#' @param pP Matrix containing the prior position parameters for the reduced form lagged coefficient matrix '\eqn{\Phi}' (default=NULL). 'NULL' indicates no priors for '\eqn{\Phi}'.
#' @param pP_sig Matrix containing values indicating confidence in 'pP' (default=NULL). 'NULL' indicates no priors for '\eqn{\Phi}'.
#' @param pR_sig Array containing values indicating confidence in long-run restrictions on the reduced form lagged coefficient matrix '\eqn{\Phi}' (default=NULL). 'Null' indicates no long-run restrictions.
#' @param kappa1 Matrix containing values indicating confidence in priors for the structural variances (default=NULL). 'NULL' indicates no priors for structural variances. 'kappa1' should have a number of columns equal to the number of endogenous variables in, or columns of, 'y'. 'kappa1' should have 1 row.
#' @param itr Integer specifying the total number of iterations for the algorithm (default=5000).
#' @param burn Integer specifying the number of draws to through out at the beginning of the algorithm (default=0).
#' @param thin Integer specifying the thinning parameter (default=1). All draws beyond burn are kept when thin = 1. Draw 1, draw 3, etc. beyond burn are kept when thin = 2.
#' @param acc_irf Boolean indicating whether or not accumulated impulse responses are to be returned (default=TRUE).
#' @param h1_irf Integer specifying the time horizon for computing impulse responses (default=12).
#' @param ci Numeric value indicating credibility intervals for the estimates to be returned (default=0.975).
#' @details Runs a Structural Bayesian Vector Autoregression model with the method developed in Baumeister and Hamilton (2015,2017,and 2018). The function returns a list containing the results.
#' @return A list containing the following:
#' @return accept_rate: Acceptance rate of the algorithm.
#' @return y and x: Matrices containing the endogenous variables and their lags.
#' @return pA, pdetA, pH, pP, pP_sig, pR, pR_sig: Matrices and arrays containing prior information.
#' @return A_max: Matrix containing estimates of the parameters in 'A' from the optimization routine.
#' @return A, detA, H, B, Phi: Arrays containing estimates of the model parameters. The first, second, and third slices of the third dimension are lower, median, and upper bounds of the estimates.
#' @return HD and IRF: Arrays containing historical decomposition of structural shocks and impulse response functions. The first, second, and thrird slices of the third dimension are lower, median, and upper bounds of the estimates.
#' @return A_den, detA_den, and H_den: Lists containing the horizontal and vertical axis coordinates of posterior densities of A, detA, and H.
#' @return Line and ACF plots of the estimates for A, det(A), and H.
#' @references Baumeister, C., and Hamilton, J.D. (2015). Sign restrictions, structural vector autoregressions, and useful prior information. \emph{Econometrica}, 83(5), 1963-1999.
#' @references Baumeister, C., and Hamilton, J.D. (2017). Structural interpretation of vector autoregressions with incomplete identification: Revisiting the role of oil supply and demand shocks (No. w24167). National Bureau of Economic Research.
#' @references Baumeister, C., and Hamilton, J.D. (2018). Inference in structural vector autoregressions when the identifying assumptions are not fully believed: Re-evaluating the role of monetary policy in economic fluctuations. \emph{Journal of Monetary Economics},
#' @examples
#' # Import data
#' library(BHSBVAR)
#' set.seed(123)
#' data(USLMData)
#' y <- matrix(data = c(USLMData$Wage, USLMData$Employment), ncol = 2)
#' colnames(y) <- c("Wage", "Employment")
#' 
#' # Set function arguments
#' nlags <- 4
#' itr <- 5000
#' burn <- 0
#' thin <- 1
#' acc_irf <- TRUE
#' h1_irf <- 20
#' ci <- 0.975
#' 
#' # Priors for A
#' pA <- array(data = NA, dim = c(2, 2, 8))
#' pA[, , 1] <- c(0, NA, 0, NA)
#' pA[, , 2] <- c(1, NA, -1, NA)
#' pA[, , 3] <- c(0.6, 1, -0.6, 1)
#' pA[, , 4] <- c(0.6, NA, 0.6, NA)
#' pA[, , 5] <- c(3, NA, 3, NA)
#' pA[, , 6] <- c(NA, NA, NA, NA)
#' pA[, , 7] <- c(NA, NA, 1, NA)
#' pA[, , 8] <- c(2.4, NA, 2.4, NA)
#' 
#' # Position priors for P
#' pP <- matrix(data = 0, nrow = ((nlags * ncol(pA)) + 1), ncol = ncol(pA))
#' pP[1:nrow(pA), 1:ncol(pA)] <-
#'   diag(x = 1, nrow = nrow(pA), ncol = ncol(pA))
#' 
#' # Confidence in the priors for P
#' x1 <- 
#'   matrix(data = NA, nrow = (nrow(y) - nlags), 
#'          ncol = (ncol(y) * nlags))
#' for (k in 1:nlags) {
#'   x1[, (ncol(y) * (k - 1) + 1):(ncol(y) * k)] <-
#'     y[(nlags - k + 1):(nrow(y) - k),]
#' }
#' x1 <- cbind(x1, 1)
#' colnames(x1) <- 
#'   c(
#'     paste(
#'       rep(colnames(y), nlags), ".L",
#'       sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(y)),
#'            decreasing = FALSE),
#'       sep = ""
#'       ),
#'     "cons"
#'     )
#' y1 <- y[(nlags + 1):nrow(y),]
#' ee <- matrix(data = NA, nrow = nrow(y1), ncol = ncol(y1))
#' for (i in 1:ncol(y1)) {
#'   xx <- cbind(x1[, seq(from = i, to = (ncol(x1) - 1), by = ncol(y1))], 1)
#'   yy <- matrix(data = y1[, i], ncol = 1)
#'   phi <- solve(t(xx) %*% xx, t(xx) %*% yy)
#'   ee[, i] <- yy - (xx %*% phi)
#' }
#' somega <- (t(ee) %*% ee) / nrow(ee)
#' lambda0 <- 0.2
#' lambda1 <- 1
#' lambda3 <- 100
#' v1 <- matrix(data = (1:nlags), nrow = nlags, ncol = 1)
#' v1 <- v1^((-2) * lambda1)
#' v2 <- matrix(data = diag(solve(diag(diag(somega)))), ncol = 1)
#' v3 <- kronecker(v1, v2)
#' v3 <- (lambda0^2) * rbind(v3, (lambda3^2))
#' v3 <- 1 / v3
#' pP_sig <- diag(x = 1, nrow = nrow(v3), ncol = nrow(v3))
#' diag(pP_sig) <- v3
#' 
#' # Confidence in priors for R (long-run restrictions)
#' pR_sig <-
#'   array(data = 0,
#'         dim = c(((nlags * ncol(y)) + 1),
#'                 ((nlags * ncol(y)) + 1),
#'                 ncol(y)))
#' Ri <-
#'   cbind(
#'     kronecker(matrix(data = 1, nrow = 1, ncol = nlags),
#'               matrix(data = c(1, 0), nrow = 1)),
#'     0)
#' pR_sig[,,2] <- (t(Ri) %*% Ri) / 0.1
#' 
#' # Confidence in priors for structural variances
#' kappa1 <- matrix(data = 2, nrow = 1, ncol = ncol(y))
#' 
#' # Set graphical parameters
#' par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
#'     mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
#' 
#' # Run the model and estimate the model parameters
#' results1 <- 
#'   BH_SBVAR(y = y, nlags = nlags, pA = pA, pP = pP, pP_sig = pP_sig,
#'            pR_sig = pR_sig, kappa1 = kappa1, itr = itr, burn = burn,
#'            thin = thin, acc_irf = acc_irf,
#'            h1_irf = h1_irf, ci = ci)
BH_SBVAR <- function(y, nlags, pA, pdetA = NULL, pH = NULL, pP = NULL, pP_sig = NULL, pR_sig = NULL, kappa1 = NULL, itr = 5000, burn = 0, thin = 1, acc_irf = TRUE, h1_irf = 12, ci = 0.975) {
  
  #construct objects from NULL inputs
  if (is.null(pdetA)) {
    pdetA <- array(data = NA_real_, dim = c(1, 1, 6))
  }
  if (is.null(pH)) {
    pH <- array(data = NA_real_, dim = c(ncol(y), ncol(y), 6))
  }
  if (is.null(pP) | is.null(pP_sig)) {
    pP <- matrix(data = 0, nrow = ((nlags * ncol(y)) + 1), ncol = ncol(y))
    pP_sig <- matrix(data = 0, nrow = ((nlags * ncol(y)) + 1), ncol = ((nlags * ncol(y)) + 1))
  }
  if (is.null(pR_sig)) {
    pR_sig <- array(data = 0, dim = c(((nlags * ncol(y)) + 1), ((nlags * ncol(y)) + 1), ncol(y)))
  }
  if (is.null(kappa1)) {
    kappa1 <- matrix(data = 0, nrow = 1, ncol = ncol(y))
  }
  
  #check BH_SBVAR function arguments
  test <- arguments_check(y, nlags, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, itr, burn, thin, acc_irf, h1_irf, ci)
  if (test != "pass") {
    stop(test)
  }
  
  #create proposal scale matrix
  scale_ar <-  diag(x = c(pA[, , 8])[which(!is.na(pA[, , 1]))], nrow = length(which(!is.na(pA[, , 1]))), ncol = length(which(!is.na(pA[, , 1]))))
  
  #trim pA
  pA <- pA[,, 1:7]
  
  #check for variable names
  if (is.null(colnames(y))) {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
  } else {
    colnames(y) <- make.names(colnames(y), unique = TRUE)
  }
  rownames(y) <- NULL
  
  #get variable names
  VarNames <- colnames(y)
  
  #get x and y data matrices
  list1 <- getXY(y, nlags)
  x1 <- list1[[1]]
  y1 <- list1[[2]]
  
  #omega
  omega <- ((t(y1) %*% y1) - (t(y1) %*% x1) %*% solve(t(x1) %*% x1) %*% t(t(y1) %*% x1)) / nrow(y1)
  
  #somega
  ee <- matrix(data = NA_real_, nrow = nrow(y1), ncol = ncol(y1), dimnames = list(rownames(y1), colnames(y1)))
  for (i in 1:ncol(y1)) {
    xx <- cbind(x1[, seq(from = i, to = (ncol(x1) - 1), by = ncol(y1))], 1)
    yy <- matrix(data = y1[, i], ncol = 1)
    phi <- solve((t(xx) %*% xx), (t(xx) %*% yy))
    ee[, i] <- yy - (xx %*% phi)
  }
  somega <- (t(ee) %*% ee) / nrow(ee)
  
  #optimization
  startvalues <- matrix(data = c(pA[, , 3])[c(which(!is.na(pA[, , 1])))], ncol = 1)
  A_optim <- stats::optim(par = c(startvalues), fn = post_A_max, pA = pA, pdetA = pdetA, pH = pH, pP = pP, pP_sig = pP_sig, pR_sig = pR_sig, kappa1 = kappa1, y1 = y1, x1 = x1, omega = omega, somega = somega, nlags = nlags, method = "BFGS", hessian = TRUE)
  #optimum values in A
  A_max <- matrix(data = NA_real_, nrow = (nrow(pA) * ncol(pA)), ncol = 1)
  A_max[c(which(!is.na(pA[, , 1]))), 1] <- A_optim[[1]]
  A_max[c(which(is.na(pA[, , 1]))), 1] <- c(pA[, , 3])[c(which(is.na(pA[, , 1])))]
  A_max <- matrix(data = A_max, nrow = nrow(pA), ncol = ncol(pA), dimnames = list(colnames(y1), colnames(y1)))
  
  #test possible sign restrictions for optimized starting values
  H_max <- solve(A_max)
  for (i in 1:nrow(pA)) {
    for (j in 1:ncol(pA)) {
      if ((pA[i, j, 1] == 0) & (!is.na(pA[i, j, 2])) & (pA[i, j, 2] != (A_max[i, j] / abs(A_max[i, j])))) {
        stop("Optimization routine produces values for the elements in A that are not cosistent with sign restrictions. Reconsider your priors.")
      }
      if ((pH[i, j, 1] == 0) & (!is.na(pH[i, j, 2])) & (pH[i, j, 2] != (H_max[i, j] / abs(H_max[i, j])))) {
        stop("Optimization routine produces values for the elements in H that are not cosistent with sign restrictions. Reconsider your priors.")
      }
    }
  }
  if ((pdetA[1, 1, 1] == 0) & (!is.na(pdetA[1, 1, 2])) & (pdetA[1, 1, 2] != (det(A_max) / abs(det(A_max))))) {
    stop("Optimization routine produces values for the determinant of A that are not consistent with sign restrictions. Reconsider your priors.")
  }
  
  #scale
  H0 <- A_optim[[6]]
  if (min(eigen(solve(H0))[[1]]) > 0) {
    PH <- t(chol(solve(H0)))
  } else {
    PH <- diag(x = 1, nrow = nrow(H0))
  }
  scale1 <- PH * scale_ar
  
  #Metropolis-Hastings Algorithm
  results <- MAIN(y1, x1, omega, somega, nlags, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, A_max, itr, burn, thin, scale1, h1_irf, acc_irf, ci, VarNames, line_plot, acf_plot)

  return(results)
}

# Check arguments from the IRF_Plots, HD_Plots, Dist_Plots functions.
#' @keywords internal
check_results <- function(results, xlab, ylab) {
  if ((!is.list(results)) || (length(results) == 0)) {
    return(paste("results: Must be a list of arrays obtained from running BH_SBVAR() function.", sep = ""))
  }
  if ((is.null(results$y)) || (!is.matrix(results$y)) || (any(!is.finite(results$y)))) {
    return(paste("results: y from BH_SBVAR() is not present", sep = ""))
  }
  if ((is.null(results$A)) || (!is.array(results$A)) || (any(!is.finite(results$A))) || (dim(results$A)[1] != dim(results$A)[2]) || (dim(results$A)[3] < 3) || (dim(results$A)[2] != ncol(results$y)) || (dim(results$A)[2] < 2)) {
    return(paste("results: A from BH_SBVAR() is not present", sep = ""))
  }
  if ((is.null(results$IRF)) || (!is.array(results$IRF)) || (any(!is.finite(results$IRF))) || (dim(results$IRF)[1] < 4) || (dim(results$IRF)[3] < 3) || (dim(results$IRF)[2] != ((dim(results$A)[2])^2))) {
    return(paste("results: IRF from BH_SBVAR() is not present", sep = ""))
  }
  if ((is.null(results$HD)) || (!is.array(results$HD)) || (dim(results$HD)[2] < 4) || (dim(results$HD)[3] < 3) || (dim(results$HD)[2] != ((dim(results$A)[2])^2))) {
    return(paste("results: HD from BH_SBVAR() is not present", sep = ""))
  }
  if ((!is.null(xlab)) && ((class(xlab) != "character") || (length(xlab) != 1))) {
    return(paste("xlab: Must be a character vector containing the label for the horizontal axis", sep = ""))
  }
  if ((!is.null(ylab)) && ((class(ylab) != "character") || (length(ylab) != 1))) {
    return(paste("ylab: Must be a character vector containing the label for the vertical axis", sep = ""))
  }
  return("pass")
}

#' Plot Impulse Responses
#' 
#' Plot Impulse Responses.
#' @author Paul Richardson
#' @export
#' @name IRF_Plots
#' @param results List containing the results from running BH_SBVAR().
#' @param varnames Character vector containing the names of the endogenous variables.
#' @param shocknames Character vector containing the names of the shocks.
#' @param xlab Character label for the horizontal axis of impulse response plots (default = NULL). Default produces plots without a label for the horizontal axis.
#' @param ylab Character label for the vertical axis of impulse response plots (default = NULL). Default produces plots without a label for the vertical axis.
#' @details Plots impulse responses and returns a list containing the actual processed data used to create the plots.
#' @return A list containing impulse responses:
#' @examples
#' # Import data
#' library(BHSBVAR)
#' set.seed(123)
#' data(USLMData)
#' y <- matrix(data = c(USLMData$Wage, USLMData$Employment), ncol = 2)
#' colnames(y) <- c("Wage", "Employment")
#' 
#' # Set function arguments
#' nlags <- 4
#' itr <- 5000
#' burn <- 0
#' thin <- 1
#' acc_irf <- TRUE
#' h1_irf <- 20
#' ci <- 0.975
#' 
#' # Priors for A
#' pA <- array(data = NA, dim = c(2, 2, 8))
#' pA[, , 1] <- c(0, NA, 0, NA)
#' pA[, , 2] <- c(1, NA, -1, NA)
#' pA[, , 3] <- c(0.6, 1, -0.6, 1)
#' pA[, , 4] <- c(0.6, NA, 0.6, NA)
#' pA[, , 5] <- c(3, NA, 3, NA)
#' pA[, , 6] <- c(NA, NA, NA, NA)
#' pA[, , 7] <- c(NA, NA, 1, NA)
#' pA[, , 8] <- c(2.4, NA, 2.4, NA)
#' 
#' # Position priors for P
#' pP <- matrix(data = 0, nrow = ((nlags * ncol(pA)) + 1), ncol = ncol(pA))
#' pP[1:nrow(pA), 1:ncol(pA)] <-
#'   diag(x = 1, nrow = nrow(pA), ncol = ncol(pA))
#' 
#' # Confidence in the priors for P
#' x1 <- 
#'   matrix(data = NA, nrow = (nrow(y) - nlags), 
#'          ncol = (ncol(y) * nlags))
#' for (k in 1:nlags) {
#'   x1[, (ncol(y) * (k - 1) + 1):(ncol(y) * k)] <-
#'     y[(nlags - k + 1):(nrow(y) - k),]
#' }
#' x1 <- cbind(x1, 1)
#' colnames(x1) <- 
#'   c(
#'     paste(
#'       rep(colnames(y), nlags), ".L",
#'       sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(y)),
#'            decreasing = FALSE),
#'       sep = ""
#'       ),
#'     "cons"
#'     )
#' y1 <- y[(nlags + 1):nrow(y),]
#' ee <- matrix(data = NA, nrow = nrow(y1), ncol = ncol(y1))
#' for (i in 1:ncol(y1)) {
#'   xx <- cbind(x1[, seq(from = i, to = (ncol(x1) - 1), by = ncol(y1))], 1)
#'   yy <- matrix(data = y1[, i], ncol = 1)
#'   phi <- solve(t(xx) %*% xx, t(xx) %*% yy)
#'   ee[, i] <- yy - (xx %*% phi)
#' }
#' somega <- (t(ee) %*% ee) / nrow(ee)
#' lambda0 <- 0.2
#' lambda1 <- 1
#' lambda3 <- 100
#' v1 <- matrix(data = (1:nlags), nrow = nlags, ncol = 1)
#' v1 <- v1^((-2) * lambda1)
#' v2 <- matrix(data = diag(solve(diag(diag(somega)))), ncol = 1)
#' v3 <- kronecker(v1, v2)
#' v3 <- (lambda0^2) * rbind(v3, (lambda3^2))
#' v3 <- 1 / v3
#' pP_sig <- diag(x = 1, nrow = nrow(v3), ncol = nrow(v3))
#' diag(pP_sig) <- v3
#' 
#' # Confidence in priors for R (long-run restrictions)
#' pR_sig <-
#'   array(data = 0,
#'         dim = c(((nlags * ncol(y)) + 1),
#'                 ((nlags * ncol(y)) + 1),
#'                 ncol(y)))
#' Ri <-
#'   cbind(
#'     kronecker(matrix(data = 1, nrow = 1, ncol = nlags),
#'               matrix(data = c(1, 0), nrow = 1)),
#'     0)
#' pR_sig[,,2] <- (t(Ri) %*% Ri) / 0.1
#' 
#' # Confidence in priors for structural variances
#' kappa1 <- matrix(data = 2, nrow = 1, ncol = ncol(y))
#' 
#' # Set graphical parameters
#' par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
#'     mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
#' 
#' # Run the model and estimate the model parameters
#' results1 <- 
#'   BH_SBVAR(y = y, nlags = nlags, pA = pA, pP = pP, pP_sig = pP_sig,
#'            pR_sig = pR_sig, kappa1 = kappa1, itr = itr, burn = burn,
#'            thin = thin, acc_irf = acc_irf,
#'            h1_irf = h1_irf, ci = ci)
#' 
#' # Plot impulse responses
#' VarNames <- colnames(USLMData)[2:3]
#' ShockNames <- c("Labor Demand","Labor Supply")
#' irf_results <- 
#'   IRF_Plots(results = results1, varnames = VarNames,
#'             shocknames = ShockNames)
IRF_Plots <- function(results, varnames, shocknames = NULL, xlab = NULL, ylab = NULL) {
  
  #test arguments
  test <- check_results(results, xlab, ylab)
  if (test != "pass") {
    stop(test)
  }
  if ((class(varnames) != "character") || (length(varnames) != dim(results$A)[2])) {
    return(paste("varnames: Must be a character vector containing the names of the endogenous variables", sep = ""))
  }
  if (is.null(shocknames)) {
    shocknames <- varnames
  }
  if ((class(shocknames) != "character") || (length(shocknames) != dim(results$A)[2])) {
    stop(paste("shocknames: Must be a character vector containing the names of the shocks", sep = ""))
  }
  
  if (is.null(xlab)) {
    xlab <- ""
  }
  if (is.null(ylab)) {
    ylab <- ""
  }
  
  IRF <- results$IRF
  nvar <- dim(results$A)[1]
  xticks <- floor(dim(IRF)[1] / 4)
  
  #store results from impulse responses
  irf_results <- list()
  for (j in 1:nvar) {
    for (i in 1:nvar) {
      #impulse responses
      irf_name <- paste("Impulse_", c(varnames[j]), "_Response_", c(varnames[i]), sep = "")
      irf_results[[(length(irf_results) + 1)]] <- IRF[, ((nvar * (j - 1)) + i), ]
      names(irf_results)[length(irf_results)] <- irf_name[1]
      #impulse response plots
      mat_ts <- stats::ts(cbind(0, irf_results[[length(irf_results)]]))
      colnames(mat_ts) <- c("Series1", "Series2", "Series3", "Series4")
      stats::ts.plot(mat_ts, col = c("black", "red", "black", "red"), gpars = list(xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "r", xaxt = "n", lty = c(1, 2, 1, 2))) + 
        graphics::title(main = paste("Response of ", varnames[i], " to ", shocknames[j], sep = ""), col.main = "black") +
        graphics::axis(side = 1, at = seq(from = 1, to = nrow(mat_ts), by = xticks), labels = seq(from = 0, to = (nrow(mat_ts) - 1),by = xticks))
    }
  }
  return(irf_results)
}

#' Plot Historical Decompositions
#' 
#' Plot Historical Decompositions.
#' @author Paul Richardson
#' @export
#' @name HD_Plots
#' @param results List containing the results from running BH_SBVAR().
#' @param varnames Character vector containing the names of the endogenous variables.
#' @param shocknames Character vector containing the names of the shocks.
#' @param xlab Character label for the horizontal axis of historical decomposition plots (default = NULL). Default produces plots without a label for the horizontal axis.
#' @param ylab Character label for the vertical axis of historical decomposition plots (default = NULL). Default produces plots without a label for the vertical axis.
#' @param freq Numeric value indicating the frequency of the data.
#' @param start_date Numeric vector indicating the starting date.
#' @details Plots historical decompositions and returns a list containing the actual processed data used to create the plots.
#' @return A list containing historical decompositions:
#' @examples
#' # Import data
#' library(BHSBVAR)
#' set.seed(123)
#' data(USLMData)
#' y <- matrix(data = c(USLMData$Wage, USLMData$Employment), ncol = 2)
#' colnames(y) <- c("Wage", "Employment")
#' 
#' # Set function arguments
#' nlags <- 4
#' itr <- 5000
#' burn <- 0
#' thin <- 1
#' acc_irf <- TRUE
#' h1_irf <- 20
#' ci <- 0.975
#' 
#' # Priors for A
#' pA <- array(data = NA, dim = c(2, 2, 8))
#' pA[, , 1] <- c(0, NA, 0, NA)
#' pA[, , 2] <- c(1, NA, -1, NA)
#' pA[, , 3] <- c(0.6, 1, -0.6, 1)
#' pA[, , 4] <- c(0.6, NA, 0.6, NA)
#' pA[, , 5] <- c(3, NA, 3, NA)
#' pA[, , 6] <- c(NA, NA, NA, NA)
#' pA[, , 7] <- c(NA, NA, 1, NA)
#' pA[, , 8] <- c(2.4, NA, 2.4, NA)
#' 
#' # Position priors for P
#' pP <- matrix(data = 0, nrow = ((nlags * ncol(pA)) + 1), ncol = ncol(pA))
#' pP[1:nrow(pA), 1:ncol(pA)] <-
#'   diag(x = 1, nrow = nrow(pA), ncol = ncol(pA))
#' 
#' # Confidence in the priors for P
#' x1 <- 
#'   matrix(data = NA, nrow = (nrow(y) - nlags), 
#'          ncol = (ncol(y) * nlags))
#' for (k in 1:nlags) {
#'   x1[, (ncol(y) * (k - 1) + 1):(ncol(y) * k)] <-
#'     y[(nlags - k + 1):(nrow(y) - k),]
#' }
#' x1 <- cbind(x1, 1)
#' colnames(x1) <- 
#'   c(
#'     paste(
#'       rep(colnames(y), nlags), ".L",
#'       sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(y)),
#'            decreasing = FALSE),
#'       sep = ""
#'       ),
#'     "cons"
#'     )
#' y1 <- y[(nlags + 1):nrow(y),]
#' ee <- matrix(data = NA, nrow = nrow(y1), ncol = ncol(y1))
#' for (i in 1:ncol(y1)) {
#'   xx <- cbind(x1[, seq(from = i, to = (ncol(x1) - 1), by = ncol(y1))], 1)
#'   yy <- matrix(data = y1[, i], ncol = 1)
#'   phi <- solve(t(xx) %*% xx, t(xx) %*% yy)
#'   ee[, i] <- yy - (xx %*% phi)
#' }
#' somega <- (t(ee) %*% ee) / nrow(ee)
#' lambda0 <- 0.2
#' lambda1 <- 1
#' lambda3 <- 100
#' v1 <- matrix(data = (1:nlags), nrow = nlags, ncol = 1)
#' v1 <- v1^((-2) * lambda1)
#' v2 <- matrix(data = diag(solve(diag(diag(somega)))), ncol = 1)
#' v3 <- kronecker(v1, v2)
#' v3 <- (lambda0^2) * rbind(v3, (lambda3^2))
#' v3 <- 1 / v3
#' pP_sig <- diag(x = 1, nrow = nrow(v3), ncol = nrow(v3))
#' diag(pP_sig) <- v3
#' 
#' # Confidence in priors for R (long-run restrictions)
#' pR_sig <-
#'   array(data = 0,
#'         dim = c(((nlags * ncol(y)) + 1),
#'                 ((nlags * ncol(y)) + 1),
#'                 ncol(y)))
#' Ri <-
#'   cbind(
#'     kronecker(matrix(data = 1, nrow = 1, ncol = nlags),
#'               matrix(data = c(1, 0), nrow = 1)),
#'     0)
#' pR_sig[,,2] <- (t(Ri) %*% Ri) / 0.1
#' 
#' # Confidence in priors for structural variances
#' kappa1 <- matrix(data = 2, nrow = 1, ncol = ncol(y))
#' 
#' # Set graphical parameters
#' par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
#'     mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
#' 
#' # Run the model and estimate the model parameters
#' results1 <- 
#'   BH_SBVAR(y = y, nlags = nlags, pA = pA, pP = pP, pP_sig = pP_sig,
#'            pR_sig = pR_sig, kappa1 = kappa1, itr = itr, burn = burn,
#'            thin = thin, acc_irf = acc_irf,
#'            h1_irf = h1_irf, ci = ci)
#' 
#' # Plot historical decompositions
#' VarNames <- colnames(USLMData)[2:3]
#' ShockNames <- c("Labor Demand","Labor Supply")
#' hd_results <- 
#'   HD_Plots(results  = results1, varnames = VarNames,
#'            shocknames = ShockNames,
#'            freq = 4, start_date = c(1971, 2))
HD_Plots <- function(results, varnames, shocknames = NULL, xlab = NULL, ylab = NULL, freq, start_date) {
  
  #test arguments
  test <- check_results(results, xlab, ylab)
  if (test != "pass") {
    stop(test)
  }
  if ((class(varnames) != "character") || (length(varnames) != dim(results$A)[2])) {
    return(paste("varnames: Must be a character vector containing the names of the endogenous variables", sep = ""))
  }
  if (is.null(shocknames)) {
    shocknames <- varnames
  }
  if ((class(shocknames) != "character") || (length(shocknames) != dim(results$A)[2])) {
    stop(paste("shocknames: Must be a character vector containing the names of the shocks", sep = ""))
  }
  if ((class(freq) != "numeric") || (!is.finite(freq)) || (length(freq) != 1) || ((freq %% 1) != 0) || (freq < 1)) {
    stop("freq: Must be a finite whole number grater than 0.")
  }
  if ((class(start_date) != "numeric") || (!is.finite(start_date)) || (length(start_date) != 2) || ((start_date %% 1) != 0) || (any(start_date < 0))) {
    stop("start_date: Must be a numeric vector containing finite whole numbers greater than or equal to 0.")
  }
  
  if (is.null(xlab)) {
    xlab <- ""
  }
  if (is.null(ylab)) {
    ylab <- ""
  }
  
  y <- results$y
  HD <- results$HD
  nvar <- dim(results$A)[1]
  
  #store results from histroical decompositions
  hd_results <- list()
  for (j in 1:nvar) {
    for (i in 1:nvar) {
      #historical decomposition
      hd_name <- paste("Contribution_", c(varnames[j]), "_On_", c(varnames[i]), sep = "")
      hd_results[[(length(hd_results) + 1)]] <- HD[, ((nvar * (j - 1)) + i), ]
      names(hd_results)[length(hd_results)] <- hd_name[1]
      #historical decomposition plots
      mat_ts <- stats::ts(cbind(0, y[,i], hd_results[[length(hd_results)]]),frequency = freq, start = start_date)
      colnames(mat_ts) <- c("Series1", "Series2", "Series3", "Series4", "Series5")
      stats::ts.plot(mat_ts, col = c("black", "black", "red", "red", "red"), gpars = list(xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "r", lty = c(1, 1, 2, 1, 2))) + 
        graphics::title(main = paste("Contribution of ", shocknames[j], " Shocks on ", varnames[i], sep = ""), col.main = "black")
    }
  }
  return(hd_results)
}

# Line Plots
#' @keywords internal
line_plot <- function(data, prior_name, i, j) {
  if (prior_name == "pA") {
    elast = -1
  } else {
    elast = 1
  }
  graphics::plot(x = (elast * data), type = "l", col = "black", yaxs = "r", xaxs = "i", xlab = "Iteration", ylab = "Estimate")
  if (prior_name == "pA") {
    graphics::title(main = paste("-A(", i, ",", j, ")", sep = ""), col.main = "black")
  } else if (prior_name == "pH") {
    graphics::title(main = paste("H(", i, ",", j, ")", sep = ""), col.main = "black")
  } else if (prior_name == "pdetA") {
    graphics::title(main = paste("Determinant of A"), col.main = "black")
  }
  Sys.sleep(0.25)
}

# ACF Plots
#' @keywords internal
acf_plot <- function(data, prior_name, i, j) {
  stats::acf(stats::ts(data), lag.max = NULL, plot = TRUE, type = c("correlation"), demean = TRUE, main = "", xlab = "Lag Length", ylab = "Correlation")
  if (prior_name == "pA") {
    graphics::title(main = paste("-A(", i, ",", j, ")", sep = ""), col.main = "black")
  } else if (prior_name == "pH") {
    graphics::title(main = paste("H(", i, ",", j, ")", sep = ""), col.main = "black")
  } else if (prior_name == "pdetA") {
    graphics::title(main = paste("Determinant of A"), col.main = "black")
  }
  Sys.sleep(0.25)
}

# Density Plots
#' @keywords internal
den_plot <- function(list2, den1, elast, lb, ub, nticks0, A_titles, H_titles, xlab, ylab, k, j, i) {
  yticks <- signif(((max(den1[, 2]) - min(den1[, 2])) / nticks0), 2)
  graphics::plot(x = (elast * den1[, 1]), y = den1[, 2], type = "l", col = "black", yaxs = "i", xaxs = "r", yaxt = "n", xlab = xlab, ylab = ylab, xlim = c(lb, ub), ylim = c(0, (yticks * (nticks0 + 1))))
  if (names(list2[k]) == "pA") {
    graphics::title(main = A_titles[i, j], col.main = "black")
  } else if (names(list2[k]) == "pH") {
    graphics::title(main = H_titles[i, j], col.main = "black")
  } else if (names(list2[k]) == "pdetA") {
    graphics::title(main = paste("Determinant of A"), col.main = "black")
  }
  graphics::axis(side = 2, at = seq(from = -yticks, to = (nticks0 * yticks), by = yticks), labels = seq(from = -yticks, to = (nticks0 * yticks), by = yticks))
  graphics::polygon(x = (elast * den1[, 1]), y = den1[, 2], col = "blue")
}

#' Plot Posterior Distributions Against Priors
#' 
#' Plot Posterior Distributions Against Priors.
#' @author Paul Richardson
#' @export
#' @import Rcpp
#' @name Dist_Plots
#' @param results List containing the results from running BH_SBVAR().
#' @param A_titles Matrix containing the titles for the plots of the elements in the coefficient matrix A.
#' @param H_titles Matrix containing the titles for the plots of the elements in the coefficient matrix H (default = NULL).
#' @param xlab Character label for the horizontal axis of historical decomposition plots (default = NULL). Default produces plots without a label for the horizontal axis.
#' @param ylab Character label for the vertical axis of historical decomposition plots (default = NULL). Default produces plots without a label for the vertical axis.
#' @details Plots posterior distributions against prior distributions.
#' @examples
#' # Import data
#' library(BHSBVAR)
#' set.seed(123)
#' data(USLMData)
#' y <- matrix(data = c(USLMData$Wage, USLMData$Employment), ncol = 2)
#' colnames(y) <- c("Wage", "Employment")
#' 
#' # Set function arguments
#' nlags <- 4
#' itr <- 5000
#' burn <- 0
#' thin <- 1
#' acc_irf <- TRUE
#' h1_irf <- 20
#' ci <- 0.975
#' 
#' # Priors for A
#' pA <- array(data = NA, dim = c(2, 2, 8))
#' pA[, , 1] <- c(0, NA, 0, NA)
#' pA[, , 2] <- c(1, NA, -1, NA)
#' pA[, , 3] <- c(0.6, 1, -0.6, 1)
#' pA[, , 4] <- c(0.6, NA, 0.6, NA)
#' pA[, , 5] <- c(3, NA, 3, NA)
#' pA[, , 6] <- c(NA, NA, NA, NA)
#' pA[, , 7] <- c(NA, NA, 1, NA)
#' pA[, , 8] <- c(2.4, NA, 2.4, NA)
#' 
#' # Position priors for P
#' pP <- matrix(data = 0, nrow = ((nlags * ncol(pA)) + 1), ncol = ncol(pA))
#' pP[1:nrow(pA), 1:ncol(pA)] <-
#'   diag(x = 1, nrow = nrow(pA), ncol = ncol(pA))
#' 
#' # Confidence in the priors for P
#' x1 <- 
#'   matrix(data = NA, nrow = (nrow(y) - nlags), 
#'          ncol = (ncol(y) * nlags))
#' for (k in 1:nlags) {
#'   x1[, (ncol(y) * (k - 1) + 1):(ncol(y) * k)] <-
#'     y[(nlags - k + 1):(nrow(y) - k),]
#' }
#' x1 <- cbind(x1, 1)
#' colnames(x1) <- 
#'   c(
#'     paste(
#'       rep(colnames(y), nlags), ".L",
#'       sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(y)),
#'            decreasing = FALSE),
#'       sep = ""
#'       ),
#'     "cons"
#'     )
#' y1 <- y[(nlags + 1):nrow(y),]
#' ee <- matrix(data = NA, nrow = nrow(y1), ncol = ncol(y1))
#' for (i in 1:ncol(y1)) {
#'   xx <- cbind(x1[, seq(from = i, to = (ncol(x1) - 1), by = ncol(y1))], 1)
#'   yy <- matrix(data = y1[, i], ncol = 1)
#'   phi <- solve(t(xx) %*% xx, t(xx) %*% yy)
#'   ee[, i] <- yy - (xx %*% phi)
#' }
#' somega <- (t(ee) %*% ee) / nrow(ee)
#' lambda0 <- 0.2
#' lambda1 <- 1
#' lambda3 <- 100
#' v1 <- matrix(data = (1:nlags), nrow = nlags, ncol = 1)
#' v1 <- v1^((-2) * lambda1)
#' v2 <- matrix(data = diag(solve(diag(diag(somega)))), ncol = 1)
#' v3 <- kronecker(v1, v2)
#' v3 <- (lambda0^2) * rbind(v3, (lambda3^2))
#' v3 <- 1 / v3
#' pP_sig <- diag(x = 1, nrow = nrow(v3), ncol = nrow(v3))
#' diag(pP_sig) <- v3
#' 
#' # Confidence in priors for R (long-run restrictions)
#' pR_sig <-
#'   array(data = 0,
#'         dim = c(((nlags * ncol(y)) + 1),
#'                 ((nlags * ncol(y)) + 1),
#'                 ncol(y)))
#' Ri <-
#'   cbind(
#'     kronecker(matrix(data = 1, nrow = 1, ncol = nlags),
#'               matrix(data = c(1, 0), nrow = 1)),
#'     0)
#' pR_sig[,,2] <- (t(Ri) %*% Ri) / 0.1
#' 
#' # Confidence in priors for structural variances
#' kappa1 <- matrix(data = 2, nrow = 1, ncol = ncol(y))
#' 
#' # Set graphical parameters
#' par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
#'     mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
#' 
#' # Run the model and estimate the model parameters
#' results1 <- 
#'   BH_SBVAR(y = y, nlags = nlags, pA = pA, pP = pP, pP_sig = pP_sig,
#'            pR_sig = pR_sig, kappa1 = kappa1, itr = itr, burn = burn,
#'            thin = thin, acc_irf = acc_irf,
#'            h1_irf = h1_irf, ci = ci)
#' 
#' # Plot Posterior and Prior Densities
#' A_titles <- 
#'   matrix(data = NA_character_, nrow = dim(pA)[1], ncol = dim(pA)[2])
#' A_titles[1, 1] <- "Wage Elasticity of Labor Demand"
#' A_titles[1, 2] <- "Wage Elasticity of Labor Supply"
#' par(mfcol = c(1, 2))
#' dist_results <- 
#'   Dist_Plots(results = results1, A_titles = A_titles)
Dist_Plots <- function(results, A_titles, H_titles = NULL, xlab = NULL, ylab = NULL) {
  
  #test arguments
  test <- check_results(results, xlab, ylab)
  if (test != "pass") {
    stop(test)
  }
  
  if (is.null(xlab)) {
    xlab <- ""
  }
  if (is.null(ylab)) {
    ylab <- ""
  }
  
  pA <- results$pA
  pdetA <- results$pdetA
  pH <- results$pH
  
  A_den <- results$A_den
  detA_den <- results$detA_den
  H_den <- results$H_den
  
  if (!is.matrix(A_titles) || ((nrow(A_titles) != dim(pA)[1]) | (ncol(A_titles) != dim(pA)[2]))) {
    stop(paste("A_titles: Must be a matrix with row and column length each equal to the number of endogenous variables.", sep = ""))
  }
  if (is.null(H_titles)) {
    H_titles <- matrix(data = NA_character_, nrow = dim(pA)[1], ncol = dim(pA)[2])
  }
  if (!is.matrix(H_titles) || ((nrow(H_titles) != dim(pH)[1]) | (ncol(H_titles) != dim(pH)[2]))) {
    stop(paste("H_titles: Must be a matrix with row and column length each equal to the number of endogenous variables.", sep = ""))
  }
  for (i in 1:dim(pA)[1]) {
    for (j in 1:dim(pA)[2]) {
      if ((is.na(pA[i,j,1])) && (!is.na(A_titles[i,j]))) {
        stop(paste("A_titles: A_titles[", i, ", ", j, "] should be empty since pA[", i, ", ", j, ", ", 1, "] is empty.", sep = ""))
      }
      if ((!is.na(pA[i,j,1])) && (is.na(A_titles[i,j]))) {
        stop(paste("A_titles: A_titles[", i, ", ", j, "] is missing.", sep = ""))
      }
      if ((is.na(pH[i,j,1])) && (!is.na(H_titles[i,j]))) {
        stop(paste("H_titles: H_titles[", i, ", ", j, "] should be empty since pH[", i, ", ", j, ", ", 1, "] is empty.", sep = ""))
      }
      if ((!is.na(pH[i,j,1])) && (is.na(H_titles[i,j]))) {
        stop(paste("H_titles: H_titles[", i, ", ", j, "] is missing.", sep = ""))
      }
    }
  }
  
  nticks0 <- 3
  
  list1 <- list(A_den = A_den, H_den = H_den)
  list2 <- list(pA = pA, pH = pH)
  for (k in 1:length(list1)) {
    if (names(list2[k]) == "pA") {
      elast <- -1
    } else {
      elast <- 1
    }
    max_distance <- 0
    distance <- 0
    for (j in 1:(dim(list2[[k]])[2])) { #equations are by column
      for (i in 1:(dim(list2[[k]])[1])) {
        if (any(!is.na(list1[[k]]$hori[i, j,]))) {
          distance <- ceiling(max(list1[[k]]$hori[i, j,], na.rm = TRUE) - min(list1[[k]]$hori[i, j,], na.rm = TRUE))
        }
        if (distance > max_distance) {
          max_distance <- distance
        }
      }
    }
    for (j in 1:(dim(list2[[k]])[2])) { #equations are by column
      for (i in 1:(dim(list2[[k]])[1])) {
        if (!is.na(list2[[k]][i, j, 1])) {
          if (is.na(list2[[k]][i,j,2])) {
            ub <- (elast * stats::median(list1[[k]]$hori[i,j,])) + (max_distance * 0.5)
            lb <- (elast * stats::median(list1[[k]]$hori[i,j,])) - (max_distance * 0.5)
          } else if (list2[[k]][i,j,2] == 1) {
            if (names(list2[k]) == "pA") {
              ub <- 0
              lb <- (-1) * max_distance
            } else {
              ub <- max_distance
              lb <- 0
            }
          } else if (list2[[k]][i,j,2] == (-1)) {
            if (names(list2[k]) == "pA") {
              ub <- max_distance
              lb <- 0
            } else {
              ub <- 0
              lb <- (-1) * max_distance
            }
          }
          den1 <- cbind(list1[[k]]$hori[i,j,],list1[[k]]$vert[i,j,])
          den_plot(list2, den1, elast, lb, ub, nticks0, A_titles, H_titles, xlab, ylab, k, j, i)
          if (list2[[k]][i, j, 1] == 0) {
            if (is.na(list2[[k]][i, j, 2])) {
              
              prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
              for (h in 1:nrow(prior_den)) {
                prior_den[h, 2] <- prior_t(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              }
              graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
              
            } else if (list2[[k]][i, j, 2] == 1) {
              
              prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
              for (h in 1:nrow(prior_den)) {
                prior_den[h, 2] <- prior_t_p(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              }
              graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
              
            } else if (list2[[k]][i, j, 2] == (-1)) {
              
              prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
              for (h in 1:nrow(prior_den)) {
                prior_den[h, 2] <- prior_t_n(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              }
              graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
              
            }
          } else if (list2[[k]][i, j, 1] == 1) {
            
            prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
            for (h in 1:nrow(prior_den)) {
              prior_den[h, 2] <- prior_nonc_t(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5], list2[[k]][i, j, 6])
            }
            graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
            
          }
        }
      }
    }
  }
  
  list2 <- list(pdetA = pdetA)
  elast <- 1
  if (!is.na(pdetA[1, 1, 1])) {
    max_distance <- ceiling(max(detA_den$hori[1, 1,], na.rm = TRUE) - min(detA_den$hori[1, 1,], na.rm = TRUE))
    
    if (list2[[1]][1,1,2] == 1) {
      ub <- max_distance
      lb <- 0
    } else if (list2[[1]][1,1,2] == (-1)) {
      ub <- 0
      lb <- (-1) * max_distance
    } else {
      ub <- (elast * stats::median(detA_den$hori[1,1,])) + (max_distance * 0.5)
      lb <- (elast * stats::median(detA_den$hori[1,1,])) - (max_distance * 0.5)
    }
    
    den1 <- cbind(detA_den$hori[1,1,],detA_den$vert[1,1,])
    den_plot(list2, den1, elast, lb, ub, nticks0, A_titles, H_titles, xlab, ylab, 1, 1, 1)
    if (pdetA[1, 1, 1] == 0) {
      if (is.na(pdetA[1, 1, 2])) {
        
        prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
        for (h in 1:nrow(prior_den)) {
          prior_den[h, 2] <- prior_t(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        }
        graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
        
      } else if (pdetA[1, 1, 2] == 1) {
        
        prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
        for (h in 1:nrow(prior_den)) {
          prior_den[h, 2] <- prior_t_p(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        }
        graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
        
      } else if (pdetA[1, 1, 2] == (-1)) {
        
        prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
        for (h in 1:nrow(prior_den)) {
          prior_den[h, 2] <- prior_t_n(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        }
        graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
        
      }
    } else if (pdetA[1, 1, 1] == 1) {
      
      prior_den <- matrix(data = seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500)), nrow = 501, ncol = 2)
      for (h in 1:nrow(prior_den)) {
        prior_den[h, 2] <- prior_nonc_t(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5], list2[[1]][1, 1, 6])
      }
      graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
      
    }
  }
}

