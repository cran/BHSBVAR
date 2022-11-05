 
# Create matrices containing dependent and independent variables.
#' @keywords internal
getXY <- function(data1, nlags) {
  varnames <- colnames(data1)
  data1 <- data1 - (matrix(data = 1, nrow = nrow(data1), ncol = ncol(data1)) %*% diag(x = colMeans(x = data1, na.rm = FALSE, dims = 1)))
  colnames(data1) <- varnames
  X <- matrix(data = NA_real_, nrow = (nrow(data1) - nlags), ncol = (ncol(data1) * nlags))
  for (k in 1:nlags) {
    X[, (ncol(data1) * (k - 1) + 1):(ncol(data1) * k)] <- data1[(nlags - k + 1):(nrow(data1) - k), ]
  }
  X <- cbind(X, 1)
  colnames(X) <- c(paste(rep(colnames(data1), nlags), "_L", sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(data1)), decreasing = FALSE), sep = ""), "cons")
  Y <- data1[(nlags + 1):nrow(data1), ]
  list1 <- list("X" = X, "Y" = Y)
  return(list1)
}

# Posterior Density Function used to Determine Starting Values for A.
#' @keywords  internal
#' @import Rcpp
post_A_optim <- function(par, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, y1, x1, omega, somega, nlags) {
  
  nrow1 <- dim(pA)[1]
  ncol1 <- dim(pA)[2]
  
  A_temp <- c(pA[, , 3])
  A_temp[which(c(!is.na(pA[, , 1])))] <- par
  A_test <- matrix(data = A_temp, nrow = nrow1, ncol = ncol1)
  
  B_test <- matrix(data = NA_real_, nrow = ((ncol1 * nlags) + 1), ncol = ncol1)
  zeta_test <- matrix(data = 0, nrow = nrow1, ncol = ncol1)
  
  pR <- array(data = 0, dim = c(((ncol1 * nlags) + 1), ncol1, ncol1))
  
  temp0 <- t(x1) %*% x1 + pP_sig
  temp1 <- t(x1) %*% y1 + pP_sig %*% pP
  temp2 <- t(y1) %*% y1 + t(pP) %*% pP_sig %*% pP
  Phi0 <- solve(temp0, temp1)
  temp4 <- temp2 - t(Phi0) %*% temp1
  temp5 <- matrix(data = NA_real_, nrow = ((ncol1 * nlags) + 1), ncol = 1)
  
  nR <- matrix(data = 0, nrow = ncol1, ncol = 1)
  for (i in 1:ncol1) {
    if (any(pR_sig[, , i] > 0)) {
      nR[i, 1] <- 1
    }
  }
  
  if (all(nR == 0)) {
    B_test <- Phi0 %*% A_test
    diag(zeta_test) <- diag(t(A_test) %*% temp4 %*% A_test)
  } else {
    for (i in 1:ncol1) {
      if (nR[i, 1] == 0) {
        B_test[, i] <- Phi0 %*% A_test[, i]
        zeta_test[i, i] <- c(t(A_test[, i]) %*% temp4 %*% A_test[, i])
      } else {
        for (j in 1:nrow1) {
          if (is.finite(pA[j, i, 7])) {
            pR[j, i, i] <- A_test[j, i]
          }
        }
        temp5[,] <- pR_sig[, , i] %*% pR[, i, i]
        B_test[, i] <- solve((temp0 + pR_sig[, , i])) %*% ((temp1 %*% A_test[, i]) + temp5)
        zeta_test[i, i] <- ((t(A_test[, i]) %*% temp2 %*% A_test[, i]) + (t(pR[,i,i]) %*% temp5)) - (t((temp1 %*% A_test[, i]) + temp5) %*% solve(temp0 + pR_sig[,,i]) %*% ((temp1 %*% A_test[, i]) + temp5))
      }
    }
  }
  
  #compute posterior density
  priors <- sum_log_prior_densities(A_test = A_test, pA = pA, pdetA = pdetA, pH = pH)
  likelihood <- log_likelihood_function(A_test = A_test, kappa1 = kappa1, y1 = y1, omega = omega, zeta_test = zeta_test, somega = somega)
  posterior <- -(priors + likelihood)
  
  return(posterior)
  
}

# Check arguments from the BH_SBVAR function that should be integers.
#' @keywords internal
check_integers <- function(list1) {
  #testing inputs that should be integers
  for (i in 1:length(list1)) {
    if (!is.null(list1[[i]])) {
      if (((!is.numeric(list1[[i]])) & (!is.integer(list1[[i]]))) || (!is.finite(list1[[i]]))) {
        return(paste(names(list1[i]), ": Must be finite 'numeric' or 'integer'.", sep = ""))
      }
      if (length(list1[[i]]) != 1) {
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
    } else {
      return(paste(names(list1[i]), ": Must be finite 'numeric' or 'integer'.", sep = ""))
    }
  }
  return("pass")
}

# Check arguments from the BH_SBVAR function that should be doubles.
#' @keywords internal
check_doubles <- function(list1) {
  #testing inputs that could be doubles
  for (i in 1:length(list1)) {
    if (!is.null(list1[[i]])) {
      if (((!is.numeric(list1[[i]])) & (!is.integer(list1[[i]]))) || (!is.finite(list1[[i]]))) {
        return(paste(names(list1[i]), ": Must be finite 'numeric' or 'integer'.", sep = ""))
      }
      if (length(list1[[i]]) != 1) {
        return(paste(names(list1[i]), ": Must be finite 'numeric' or 'integer'.", sep = ""))
      }
      if ((names(list1[i]) == "ci") && ((list1[[i]] < 0.7) | (list1[[i]] > 1))) {
        return(paste(names(list1[i]), ": Must be greater than or equal to 0.7 and less than or equal to 1.", sep = ""))
      }
      if ((names(list1[i]) == "cri") && ((list1[[i]] < 0.4) | (list1[[i]] > 1))) {
        return(paste(names(list1[i]), ": Must be greater than or equal to 0.4 and less than or equal to 1.", sep = ""))
      }
    } else {
      if (names(list1[i]) == "cri") {
        return(paste(names(list1[i]), ": Must be finite 'numeric' or 'integer'.", sep = ""))
      }
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
    if ((names(list1[i]) == "pA") && (all(is.na(list1[[i]][, , 1])))) {
      return(paste(names(list1[i]), "[, , 1]: Should indicate at least one parameter to be estimated.", sep = ""))
    }
    if ((names(list1[i]) == "pA") && ((dim(list1[[i]])[1] != ncol(y)) | (dim(list1[[i]])[2] != ncol(y)) | (dim(list1[[i]])[3] != 8))) {
      return(paste(names(list1[i]), ": Should be a (", ncol(y), ", ", ncol(y), ", 8) array.", sep = ""))
    }
    if ((names(list1[i]) == "pdetA") && ((dim(list1[[i]])[1] != 1) | (dim(list1[[i]])[2] != 1) | (dim(list1[[i]])[3] != 6))) {
      return(paste(names(list1[i]), ": Should be a (1, 1, 6) array.", sep = ""))
    }
    if ((names(list1[i]) == "pH") && ((dim(list1[[i]])[1] != ncol(y)) | (dim(list1[[i]])[2] != ncol(y)) | (dim(list1[[i]])[3] != 6))) {
      return(paste(names(list1[i]), ": Should be a (", ncol(y), ", ", ncol(y), ", 6) array.", sep = ""))
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
        } else if (list1[[i]][j, k, 1] == 0) { #if distribution is 0 (symmetric t-distribution), parameters in slices 2:5 must be specified
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
          if ((!is.na(list1[[i]][j, k, 6]))) {
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
          if ((names(list1[i]) == "pA") && ((!is.na(list1[[i]][j, k, 7])) && ((!is.finite(list1[[i]][j, k, 7])) || (list1[[i]][j, k, 7] != 1)))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 7]: Long-run restriction should be indicated with an NA (no long-run restriction) or a 1 (long-run restriction).", sep = ""))
          }
          if ((names(list1[i]) == "pA") && ((is.na(list1[[i]][j, k, 8])) || (!is.finite(list1[[i]][j, k, 8])) || (list1[[i]][j, k, 8] < 0.1))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 8]: Proposal scaling parameter should be greater than or equal to 0.1.", sep = ""))
          }
        } else if ((list1[[i]][j, k, 1] == 2) | (list1[[i]][j, k, 1] == 3)) { #if distribution is 2 or 3 (inverted beta-distribution or beta-distribution), parameters in slices 2:5 must be specified
          if ((is.na(list1[[i]][j, k, 2])) || ((list1[[i]][j, k, 2] != 1) & (list1[[i]][j, k, 2] != (-1)))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 2]: Sign should be indicated with a 1, or -1.", sep = ""))
          }
          if ((!is.finite(list1[[i]][j, k, 4])) || ((list1[[i]][j, k, 4] <= 1))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 4]: Shape1 should be indicated with a finite number greater than 1.", sep = ""))
          }
          if ((!is.finite(list1[[i]][j, k, 5])) || ((list1[[i]][j, k, 5] <= 1))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 5]: Shape2 should be indicated with a finite number greater than 1.", sep = ""))
          }
          if (!is.na(list1[[i]][j, k, 6])) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 6]: Should be NA.", sep = ""))
          }
          if (list1[[i]][j, k, 1] == 2) {
            if (round(x = list1[[i]][j, k, 3], digits = 2) != round(x = (list1[[i]][j, k, 2] * ((list1[[i]][j, k, 4] + list1[[i]][j, k, 5]) / (list1[[i]][j, k, 5] + 1))), digits = 2)) {
              return(paste(names(list1[i]), "[", j, ", ", k, ", 3]: Position should be (", list1[[i]][j, k, 2], ") * (", names(list1[i]), "[", j, ", ", k, ", 4] + ", names(list1[i]), "[", j, ", ", k, ", 5]) / (", names(list1[i]), "[", j, ", ", k, ", 5] + 1).", sep = ""))
            }
          } else if (list1[[i]][j, k, 1] == 3) {
            if (round(x = list1[[i]][j, k, 3], digits = 2) != round(x = (list1[[i]][j, k, 2] * ((list1[[i]][j, k, 4] - 1) / (list1[[i]][j, k, 4] + list1[[i]][j, k, 5] - 2))), digits = 2)) {
              return(paste(names(list1[i]), "[", j, ", ", k, ", 3]: Position should be (", list1[[i]][j, k, 2], ") * (", names(list1[i]), "[", j, ", ", k, ", 4] - 1) / (", names(list1[i]), "[", j, ", ", k, ", 4] + ", names(list1[i]), "[", j, ", ", k, ", 5] - 2).", sep = ""))
            }
          }
          if ((names(list1[i]) == "pA") && ((!is.na(list1[[i]][j, k, 7])) && ((!is.finite(list1[[i]][j, k, 7])) || (list1[[i]][j, k, 7] != 1)))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 7]: Long-run restriction should be indicated with an NA (no long-run restriction) or a 1 (long-run restriction).", sep = ""))
          }
          if ((names(list1[i]) == "pA") && ((is.na(list1[[i]][j, k, 8])) || (!is.finite(list1[[i]][j, k, 8])) || (list1[[i]][j, k, 8] < 0.1))) {
            return(paste(names(list1[i]), "[", j, ", ", k, ", 8]: Proposal scaling parameter should be greater than or equal to 0.1.", sep = ""))
          }
        } else {
          return(paste(names(list1[i]), "[", j, ", ", k, ", 1]: Distribution should be indicated with a NA (no prior), 0 (symetric t-distribution), 1 (non-central t-distribution), or 2 (inverted beta-distribution).", sep = ""))
        }
      }
    }
  }
  
  if ((is.finite(list1[["pdetA"]][1, 1, 3])) && (list1[["pdetA"]][1, 1, 3] == 0)) {
    return("pdetA[1, 1, 3]: Should be a finite value other than 0.")
  }
  list2 <- list(pH = list1[["pH"]], pdetA = list1[["pdetA"]])
  list3 <- list(H = solve(list1[["pA"]][, , 3]), detA = matrix(data = det(list1[["pA"]][, , 3])))
  
  den1 <- 0
  for (i in 1:length(list2)) {
    if (any(is.finite(list2[[i]][, , 1]))) {
      for (j in 1:(dim(list2[[i]])[1])) {
        for (k in 1:(dim(list2[[i]])[2])) {
          
          if (is.finite(list2[[i]][j, k, 1])) {
            if (list2[[i]][j, k, 1] == 0) {
              if (!is.finite(list2[[i]][j, k, 2])) {
                den1 <- prior_t(list3[[i]][j, k], list2[[i]][j, k, 3], list2[[i]][j, k, 4], list2[[i]][j, k, 5])
                if (den1 == 0) {
                  if (names(list2[i]) == "pH") {
                    return(paste0("pH[", j, ", ", k, ", ]: The (", j, ", ", k, ") element of the inverse of pA[, , 3] produces a prior density of 0."))
                  }
                  if (names(list2[i]) == "pdetA") {
                    return(paste0("pdetA[", j, ", ", k, ", ]: The determinant of pA[, , 3] produces a prior density of 0."))
                  }
                }
              } else {
                if ((list3[[i]][j, k] != 0) && (list2[[i]][j, k, 2] != (list3[[i]][j, k] / abs(list3[[i]][j, k])))) {
                  if (names(list2[i]) == "pH") {
                    return(paste0("pH[", j, ", ", k, ", 2:3]: Must have the same sign as the (", j, ", ", k, ") element of the inverse of pA[, , 3]."))
                  }
                  if (names(list2[i]) == "pdetA") {
                    return(paste0("pdetA[", j, ", ", k, ", 2:3]: Must have the same sign as the determinant of pA[, , 3]."))
                  }
                }
                if (list2[[i]][j, k, 2] == 1) {
                  den1 <- prior_t_p(list3[[i]][j, k], list2[[i]][j, k, 3], list2[[i]][j, k, 4], list2[[i]][j, k, 5])
                }
                if (list2[[i]][j, k, 2] == -1) {
                  den1 <- prior_t_n(list3[[i]][j, k], list2[[i]][j, k, 3], list2[[i]][j, k, 4], list2[[i]][j, k, 5])
                }
                if (den1 == 0) {
                  if (names(list2[i]) == "pH") {
                    return(paste0("pH[", j, ", ", k, ", ]: The (", j, ", ", k, ") element of the inverse of pA[, , 3] produces a prior density of 0."))
                  }
                  if (names(list2[i]) == "pdetA") {
                    return(paste0("pdetA[", j, ", ", k, ", ]: The determinant of pA[, , 3] produces a prior density of 0."))
                  }
                }
              }
            }
            
            if (list2[[i]][j, k, 1] == 1) {
              den1 <- prior_nonc_t(list3[[i]][j, k], list2[[i]][j, k, 3], list2[[i]][j, k, 4], list2[[i]][j, k, 5], list2[[i]][j, k, 6])
              if (den1 == 0) {
                if (names(list2[i]) == "pH") {
                  return(paste0("pH[", j, ", ", k, ", ]: The (", j, ", ", k, ") element of the inverse of pA[, , 3] produces a prior density of 0."))
                }
                if (names(list2[i]) == "pdetA") {
                  return(paste0("pdetA[", j, ", ", k, ", ]: The determinant of pA[, , 3] produces a prior density of 0."))
                }
              }
            }
            
            if ((list2[[i]][j, k, 1] == 2) | (list2[[i]][j, k, 1] == 3)) {
              if ((list3[[i]][j, k] != 0) && (list2[[i]][j, k, 2] != (list3[[i]][j, k] / abs(list3[[i]][j, k])))) {
                if (names(list2[i]) == "pH") {
                  return(paste0("pH[", j, ", ", k, ", 2:3]: Must have the same sign as the (", j, ", ", k, ") element of the inverse of pA[, , 3]."))
                }
                if (names(list2[i]) == "pdetA") {
                  return(paste0("pdetA[", j, ", ", k, ", 2:3]: Must have the same sign as the determinant of pA[, , 3]."))
                }
              }
              if (list2[[i]][j, k, 1] == 2) {
                den1 <- prior_ibeta((list2[[i]][j, k, 2] * list3[[i]][j, k]), list2[[i]][j, k, 4], list2[[i]][j, k, 5])
              }
              if (list2[[i]][j, k, 1] == 3) {
                den1 <- prior_beta((list2[[i]][j, k, 2] * list3[[i]][j, k]), list2[[i]][j, k, 4], list2[[i]][j, k, 5])
              }
              if (den1 == 0) {
                if (names(list2[i]) == "pH") {
                  return(paste0("pH[", j, ", ", k, ", ]: The (", j, ", ", k, ") element of the inverse of pA[, , 3] produces a prior density of 0."))
                }
                if (names(list2[i]) == "pdetA") {
                  return(paste0("pdetA[", j, ", ", k, ", ]: The determinant of pA[, , 3] produces a prior density of 0."))
                }
              }
            }
            
            den1 <- 0
          }
          
        }
      }
    }
  }
  
  return("pass")
}

# Check arguments from the BH_SBVAR function
#' @keywords internal
arguments_check <- function(y, nlags, pA, pdetA, pH, pP, pP_sig, pR_sig, kappa1, itr, burn, thin, cri) {
  test <- check_integers(list1 = list(nlags = nlags, itr = itr, burn = burn, thin = thin))
  if (test != "pass") {
    return(test)
  }
  test <- check_doubles(list1 = list(cri = cri))
  if (test != "pass") {
    return(test)
  }
  if (floor((itr - burn) / thin) < 5000) {
    return(paste("'floor((itr-burn)/thin)' must be greater than or equal to 5000.", sep = ""))
  }
  test <- check_matrices(list1 = list(y = y, pP = pP, pP_sig = pP_sig, kappa1 = kappa1), nlags)
  if (test != "pass") {
    return(test)
  }
  test <- check_arrays(list1 = list(pA = pA, pdetA = pdetA, pH = pH), y)
  if (test != "pass") {
    return(test)
  }
  
  # check pR_sig
  if (!is.array(pR_sig)) {
    return("pR_sig: Must be an array.")
  }
  if (any(!is.finite(pR_sig)) || (any(pR_sig < 0))) {
    return("pR_sig: Must contain finite values greater than or equal to 0.")
  }
  if ((dim(pR_sig)[1] != ((nlags * ncol(y)) + 1)) | (dim(pR_sig)[2] != ((nlags * ncol(y)) + 1)) | (dim(pR_sig)[3] != ncol(y))) {
    return(paste("pR_sig: Dimensions should be (", ((nlags * ncol(y)) + 1), ", ", ((nlags * ncol(y)) + 1), ", ", (ncol(y)), ").", sep = ""))
  }
  for (i in 1:ncol(y)) {
    if (any(is.finite(pA[, i, 7]))) {
      n <- which(is.finite(pA[, i, 7]))
      for (j in 1:ncol(y)) {
        if (any(n == j)) {
          if (any(pR_sig[-seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), i] != 0)) {
            return(paste0("pR_sig[-c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), ", i, "]: Must be 0."))
          }
          if (any(pR_sig[seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), -seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), i] != 0)) {
            return(paste0("pR_sig[c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), -c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), ", i, "]: Must be 0."))
          }
          if (any(pR_sig[seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), i] == 0)) {
            return(paste0("pR_sig[-c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), ", i, "]: Must be greater than 0 since pA[", j, ", ", i, ", 7] = ", pA[j, i, 7], "."))
          }
        } else {
          if (any(pR_sig[, seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), i] != 0)) {
            return(paste0("pR_sig[, c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), ", i, "]: Must be 0 since pA[", j, ", ", i, ", 7] = ", pA[j, i, 7], "."))
          }
          if (any(pR_sig[seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), , i] != 0)) {
            return(paste0("pR_sig[c(", gsub(pattern = " ", replacement = ", ", x = trimws(paste0(seq(from = j, to = (nlags * ncol(y)), by = ncol(y)), " ", collapse = ""), which = "right")), "), ,", i, "]: Must be 0 since pA[", j, ", ", i, ", 7] = ", pA[j, i, 7], "."))
          }
        }
      }
    } else {
      if (any(pR_sig[, , i] != 0)) {
        return(paste0("pR_sig[, , ", i, "]: Must be 0 since pA[, ", i, ", 7] = ", pA[, i, 7], "."))
      }
    }
    if (!isSymmetric(pR_sig[, , i])) {
      return(paste0("pR_sig[, , ", i, "]: Must be symmetric."))
    }
  }
  
  return("pass")
}

# Line Plot
#' @keywords internal
line_plot <- function(data1, prior_name, i, j) {
  if (any(data1 != data1[1])) {
    if (prior_name == "pA") {
      elast = -1
    } else {
      elast = 1
    }
    graphics::plot(x = (elast * data1), type = "l", col = "black", yaxs = "r", xaxs = "i", xlab = "Iteration", ylab = "Estimate")
    if (prior_name == "pA") {
      graphics::title(main = paste("-A(", i, "," , j, ")", sep = ""), col.main = "black")
    } else if (prior_name == "pH") {
      graphics::title(main = paste("H(", i, "," , j, ")", sep = ""), col.main = "black")
    } else if (prior_name == "pdetA") {
      graphics::title(main = paste("Determinant of A"), col.main = "black")
    }
    Sys.sleep(0.25)
  }
}

# Line Plots
#' @keywords internal
Line_Plots <- function(raw_array, priors_array, prior_name) {
  n_rows <- nrow(raw_array)
  n_cols <- ncol(raw_array)
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      if (all(is.finite(raw_array[i, j, ]))) {
        line_plot(data1 = raw_array[i, j, ], prior_name = prior_name, i = i, j = j)
      }
    }
  }
}

# Autocorrelation Plot
#' @keywords internal
acf_plot <- function(data1, prior_name, i, j) {
  if (any(data1 != data1[1])) {
    stats::acf(x = stats::ts(data1), lag.max = NULL, plot = TRUE, type = c("correlation"), demean = TRUE, main = "", xlab = "Lag Length", ylab = "Correlation", ci = 0)
    if (prior_name == "pA") {
      graphics::title(main = paste("-A(", i, "," , j, ")", sep = ""), col.main = "black")
    } else if (prior_name == "pH") {
      graphics::title(main = paste("H(", i, "," , j, ")", sep = ""), col.main = "black")
    } else if (prior_name == "pdetA") {
      graphics::title(main = paste("Determinant of A"), col.main = "black")
    }
    Sys.sleep(0.25)
  }
}

# Autocorrelation Plots
#' @keywords internal
ACF_Plots <- function(raw_array, priors_array, prior_name) {
  n_rows <- nrow(raw_array)
  n_cols <- ncol(raw_array)
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      if (all(is.finite(raw_array[i, j, ]))) {
        acf_plot(data1 = raw_array[i, j, ], prior_name = prior_name, i = i, j = j)
      }
    }
  }
}

#' Structural Bayesian Vector Autoregression
#' 
#' Estimates the parameters of a Structural Bayesian Vector Autoregression model with the method developed by Baumeister and Hamilton (2015/2017/2018).
#' @author Paul Richardson
#' @export
#' @import Rcpp
#' @name BH_SBVAR
#' @param y \emph{(T x n)} matrix containing the endogenous variables. \emph{T} is the number of observations and \emph{n} is the number of endogenous variables.
#' @param nlags Integer specifying the lag order.
#' @param pA \emph{(n x n x 8)} array where \emph{n} is the number of endogenous variables and each slice of the third dimension contains the prior distributions (NA - no prior, 0 - symmetric t-distribution, 1 - non-central t-distribution, 2 - inverted beta distribution, 3 - beta distribution), sign restrictions (NA - no restriction, 1 - positive restriction, -1 - negative restriction), distribution position parameters, distribution scale or shape1 parameters for t-distributions or inverted beta and beta distributions, distribution degrees of freedom or shape2 parameters for t-distributions or inverted beta and beta distributions, distribution skew parameters for t-distributions, indication for long-run restrictions (NA - no long-run restriction, 1 - long-run restriction), and random-walk proposal scale parameters for \emph{A}, respectively.
#' @param pdetA \emph{(1 x 1 x 6)} array where each slice of the third dimension contains the prior distributions (NA - no prior, 0 - symmetric t-distribution, 1 - non-central t-distribution, 2 - inverted beta distribution, 3 - beta distribution), sign restrictions (NA - no restriction, 1 - positive restriction, -1 - negative restriction), distribution position parameters, distribution scale or shape1 parameters for t-distributions or inverted beta and beta distributions, distribution degrees of freedom or shape2 parameters for t-distributions or inverted beta and beta distributions, and distribution skew parameters for t-distributions for the determinant of \emph{A}, respectively (default = NULL). NULL indicates no priors for the determinant of \emph{A}.
#' @param pH \emph{(n x n x 6)} array where \emph{n} is the number of endogenous variables and each slice of the third dimension contains the prior distributions (NA - no prior, 0 - symmetric t-distribution, 1 - non-central t-distribution, 2 - inverted beta distribution, 3 - beta distribution), sign restrictions (NA - no restriction, 1 - positive restriction, -1 - negative restriction), distribution position parameters, distribution scale or shape1 parameters for t-distributions or inverted beta and beta distributions, distribution degrees of freedom or shape2 parameters for t-distributions or inverted beta and beta distributions, and distribution skew parameters for t-distributions for \emph{H}, the inverse of \emph{A}, respectively (default = NULL). NULL indicates no priors for the inverse of \emph{A}.
#' @param pP \emph{(k x n)} matrix containing the prior position parameters for the reduced form lagged coefficient matrix \emph{\eqn{\Phi}} (default = NULL). \emph{\eqn{k = n L + 1}}, \emph{n} is the number of endogenous variables, and \emph{L} is the lag length. NULL indicates no priors for \emph{\eqn{\Phi}}.
#' @param pP_sig \emph{(k x k)} matrix containing values indicating confidence in the priors for \emph{\eqn{\Phi}} (default = NULL). \emph{\eqn{k = n L + 1}}, \emph{n} is the number of endogenous variables, and \emph{L} is the lag length. NULL indicates no priors for \emph{\eqn{\Phi}}.
#' @param pR_sig \emph{(k x k x n)} array containing values indicating confidence in long-run restrictions on the lagged structural coefficient matrix \emph{B} (default = NULL). \emph{\eqn{k = n L + 1}}, \emph{n} is the number of endogenous variables, and \emph{L} is the lag length. NULL indicates no long-run restrictions.
#' @param kappa1 \emph{(1 x n)} matrix containing values indicating confidence in priors for the structural variances (default = NULL). \emph{n} is the number of endogenous variables. NULL indicates no priors for structural variances.
#' @param itr Integer specifying the total number of iterations for the algorithm (default = 5000).
#' @param burn Integer specifying the number of draws to throw out at the beginning of the algorithm (default = 0).
#' @param thin Integer specifying the thinning parameter (default = 1). All draws beyond burn are kept when thin = 1. Draw 1, draw 3, etc. beyond burn are kept when thin = 2.
#' @param cri credibility intervals for the estimates to be returned (default = 0.95). A value of 0.95 will return 95\% credibility intervals. A value of 0.90 will return 90\% credibility intervals.
#' @details Estimates the parameters of a Structural Bayesian Vector Autoregression model with the method developed in Baumeister and Hamilton (2015/2017/2018). The function returns a list containing the results.
#' @return A list containing the following:
#' @return accept_rate: Acceptance rate of the algorithm.
#' @return y and x: Matrices containing the endogenous variables and their lags.
#' @return nlags: Numeric value indicating the number of lags included in the model.
#' @return pA, pdetA, pH, pP, pP_sig, pR, pR_sig, tau1, and kappa1: Matrices and arrays containing prior information.
#' @return A_start: Matrix containing estimates of the parameters in \emph{A} from the optimization routine.
#' @return A, detA, H, B, Phi, and D: Arrays containing estimates of the model parameters. The first, second, and third slices of the third dimension are lower, median, and upper bounds of the estimates.
#' @return A_den, detA_den, and H_den: Lists containing the horizontal and vertical axis coordinates of posterior densities of \emph{A}, \emph{det(A)}, and \emph{H}.
#' @return A_chain, B_chain, D_chain, detA_chain, H_chain: Arrays containing the raw results for \emph{A}, \emph{B}, \emph{D}, \emph{detA}, \emph{H}.
#' @return Line and ACF plots of the estimates for \emph{A}, \emph{det(A)}, and \emph{H}.
#' @references Baumeister, C., & Hamilton, J.D. (2015). Sign restrictions, structural vector autoregressions, and useful prior information. \emph{Econometrica}, 83(5), 1963-1999.
#' @references Baumeister, C., & Hamilton, J.D. (2017). Structural interpretation of vector autoregressions with incomplete identification: Revisiting the role of oil supply and demand shocks (No. w24167). National Bureau of Economic Research.
#' @references Baumeister, C., & Hamilton, J.D. (2018). Inference in structural vector autoregressions when the identifying assumptions are not fully believed: Re-evaluating the role of monetary policy in economic fluctuations. \emph{Journal of Monetary Economics}, 100, 48-65.
#' @seealso Dr. Christiane Baumeister's website \href{https://sites.google.com/site/cjsbaumeister/}{https://sites.google.com/site/cjsbaumeister/}.
#' @seealso Dr. James D. Hamilton's website \href{https://econweb.ucsd.edu/~jhamilton/}{https://econweb.ucsd.edu/~jhamilton/}.
#' @examples
#' # Import data
#' library(BHSBVAR)
#' set.seed(123)
#' data(USLMData)
#' y0 <- matrix(data = c(USLMData$Wage, USLMData$Employment), ncol = 2)
#' y <- y0 - (matrix(data = 1, nrow = nrow(y0), ncol = ncol(y0)) %*% 
#'              diag(x = colMeans(x = y0, na.rm = FALSE, dims = 1)))
#' colnames(y) <- c("Wage", "Employment")
#' 
#' # Set function arguments
#' nlags <- 8
#' itr <- 5000
#' burn <- 0
#' thin <- 1
#' cri <- 0.95
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
#' pA[, , 8] <- c(2, NA, 2, NA)
#' 
#' # Position priors for Phi
#' pP <- matrix(data = 0, nrow = ((nlags * ncol(pA)) + 1), ncol = ncol(pA))
#' pP[1:nrow(pA), 1:ncol(pA)] <-
#'   diag(x = 1, nrow = nrow(pA), ncol = ncol(pA))
#' 
#' # Confidence in the priors for Phi
#' x1 <- 
#'   matrix(data = NA, nrow = (nrow(y) - nlags), 
#'          ncol = (ncol(y) * nlags))
#' for (k in 1:nlags) {
#'   x1[, (ncol(y) * (k - 1) + 1):(ncol(y) * k)] <-
#'     y[(nlags - k + 1):(nrow(y) - k),]
#' }
#' x1 <- cbind(x1, 1)
#' colnames(x1) <- 
#'   c(paste(rep(colnames(y), nlags),
#'           "_L",
#'           sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(y)),
#'                decreasing = FALSE),
#'           sep = ""),
#'     "cons")
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
#' # Confidence in long-run restriction priors
#' pR_sig <-
#'   array(data = 0,
#'         dim = c(((nlags * ncol(y)) + 1),
#'                 ((nlags * ncol(y)) + 1),
#'                 ncol(y)))
#' Ri <-
#'   cbind(kronecker(matrix(data = 1, nrow = 1, ncol = nlags),
#'                   matrix(data = c(1, 0), nrow = 1)),
#'         0)
#' pR_sig[, , 2] <- (t(Ri) %*% Ri) / 0.1
#' 
#' # Confidence in priors for D
#' kappa1 <- matrix(data = 2, nrow = 1, ncol = ncol(y))
#' 
#' # Set graphical parameters
#' par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
#'     mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
#' 
#' # Estimate the parameters of the model
#' results1 <- 
#'   BH_SBVAR(y = y, nlags = nlags, pA = pA, pP = pP, pP_sig = pP_sig,
#'            pR_sig = pR_sig, kappa1 = kappa1, itr = itr, burn = burn,
#'            thin = thin, cri = cri)
BH_SBVAR <- function(y, nlags, pA, pdetA = NULL, pH = NULL, pP = NULL, pP_sig = NULL, pR_sig = NULL, kappa1 = NULL, itr = 5000, burn = 0, thin = 1, cri = 0.95) {
  
  results <- tryCatch(
    expr = {
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
      test <- arguments_check(y = y, nlags = nlags, pA = pA, pdetA = pdetA, pH = pH, pP = pP, pP_sig = pP_sig, pR_sig = pR_sig, kappa1 = kappa1, itr = itr, burn = burn, thin = thin, cri = cri)
      if (test != "pass") {
        stop(test)
      }
      
      ci <- (1.0 - ((1.0 - cri) / 2.0))
      
      #create proposal scale matrix
      scale_ar <-  diag(x = c(pA[, , 8])[which(!is.na(c(pA[, , 1])))], nrow = length(which(!is.na(c(pA[, , 1])))), ncol = length(which(!is.na(c(pA[, , 1])))))
      
      #trim pA
      pA <- pA[, , 1:7]
      
      #check for variable names
      if (is.null(colnames(y))) {
        colnames(y) <- paste("y", 1:ncol(y), sep = "")
      } else {
        colnames(y) <- make.names(names = colnames(y), unique = TRUE)
      }
      rownames(y) <- NULL
      
      #get variable names
      varnames <- colnames(y)
      
      #get x and y data matrices
      list1 <- getXY(data1 = y, nlags = nlags)
      x1 <- list1$X
      y1 <- list1$Y
      
      #omega
      omega <- ((t(y1) %*% y1) - (t(y1) %*% x1) %*% solve(t(x1) %*% x1) %*% t(t(y1) %*% x1)) / nrow(y1)
      
      #somega
      ee <- matrix(data = NA_real_, nrow = nrow(y1), ncol = ncol(y1), dimnames = list(rownames(y1), varnames))
      for (i in 1:ncol(y1)) {
        xx <- cbind(x1[, seq(from = i, to = (ncol(x1) - 1), by = ncol(y1))], 1)
        yy <- matrix(data = y1[, i], ncol = 1)
        phi <- solve((t(xx) %*% xx), (t(xx) %*% yy))
        ee[, i] <- yy - (xx %*% phi)
      }
      somega <- (t(ee) %*% ee) / nrow(ee)
      
      #optimization
      startvalues <- c(pA[, , 3])[which(!is.na(c(pA[, , 1])))]
      
      A_optim <- list(par = NULL, value = NULL, counts = NULL, convergence = 1, message = NULL, hessian = diag(x = 1, nrow = length(startvalues), ncol = length(startvalues)))
      A_optim[c("par", "value", "counts", "convergence", "message")] <- stats::optim(par = startvalues, fn = post_A_optim, pA = pA, pdetA = pdetA, pH = pH, pP = pP, pP_sig = pP_sig, pR_sig = pR_sig, kappa1 = kappa1, y1 = y1, x1 = x1, omega = omega, somega = somega, nlags = nlags, method = "Nelder-Mead", control = list(maxit = 2500))[c("par", "value", "counts", "convergence", "message")]
      if (A_optim$convergence != 0) {
        if ((all(!is.null(A_optim$par))) && (all(is.finite(A_optim$par)))) {
          A_optim[c("par", "value", "counts", "convergence", "message")] <- stats::optim(par = A_optim$par, fn = post_A_optim, pA = pA, pdetA = pdetA, pH = pH, pP = pP, pP_sig = pP_sig, pR_sig = pR_sig, kappa1 = kappa1, y1 = y1, x1 = x1, omega = omega, somega = somega, nlags = nlags, method = "Nelder-Mead", control = list(maxit = 500))[c("par", "value", "counts", "convergence", "message")]
        } else {
          stop("Optimization routine convergence was not successful.")
        }
      }
      if (A_optim$convergence != 0) {
        stop("Optimization routine convergence was not successful.")
      }
      A_optim$hessian <- stats::optimHess(par = A_optim$par, fn = post_A_optim, pA = pA, pdetA = pdetA, pH = pH, pP = pP, pP_sig = pP_sig, pR_sig = pR_sig, kappa1 = kappa1, y1 = y1, x1 = x1, omega = omega, somega = somega, nlags = nlags)
      
      #optimum values in A
      A_temp <- c(pA[, , 3])
      A_temp[which(!is.na(c(pA[, , 1])))] <- A_optim$par
      A_start <- matrix(data = A_temp, nrow = dim(pA)[1], ncol = dim(pA)[2], dimnames = list(varnames, varnames))
      
      #test that optimized starting values are consistent with sign restrictions
      H_max <- solve(A_start)
      for (i in 1:nrow(pA)) {
        for (j in 1:ncol(pA)) {
          if ((!is.na(pA[i, j, 1])) && ((pA[i, j, 1] == 0) | (pA[i, j, 1] == 2) | (pA[i, j, 1] == 3)) && (!is.na(pA[i, j, 2])) && (A_start[i, j] != 0) && (pA[i, j, 2] != (A_start[i, j] / abs(A_start[i, j])))) {
            warning("Optimization routine produced values for the elements in A that are not consistent with sign restrictions.", immediate. = TRUE)
          }
          if ((!is.na(pH[i, j, 1])) && ((pH[i, j, 1] == 0) | (pH[i, j, 1] == 2) | (pH[i, j, 1] == 3)) && (!is.na(pH[i, j, 2])) && (H_max[i, j] != 0) && (pH[i, j, 2] != (H_max[i, j] / abs(H_max[i, j])))) {
            warning("Optimization routine produced values for the elements in H that are not consistent with sign restrictions.", immediate. = TRUE)
          }
        }
      }
      if ((!is.na(pdetA[1, 1, 1])) && ((pdetA[1, 1, 1] == 0) | (pdetA[1, 1, 1] == 2) | (pdetA[1, 1, 1] == 3)) && (!is.na(pdetA[1, 1, 2])) && (pdetA[1, 1, 2] != (det(A_start) / abs(det(A_start))))) {
        warning("Optimization routine produced values for the determinant of A that are not consistent with sign restrictions.", immediate. = TRUE)
      }
      
      #scale
      H0 <- A_optim$hessian
      if (min(eigen(solve(H0))[[1]]) > 0) {
        PH <- t(chol(solve(H0)))
      } else {
        PH <- diag(x = 1, nrow = nrow(H0))
      }
      scale1 <- PH * scale_ar
      
      #Metropolis-Hastings Algorithm
      results <- MAIN(y1 = y1, x1 = x1, omega = omega, somega = somega, nlags = nlags, pA = pA, pdetA = pdetA, pH = pH, pP = pP, pP_sig = pP_sig, pR_sig = pR_sig, kappa1 = kappa1, A_start = A_start, itr = itr, burn = burn, thin = thin, scale1 = scale1, ci = ci)
      
      dimnames(results$y) <- dimnames(y1)
      dimnames(results$x) <- dimnames(x1)
      
      dimnames(results$pA) <- list(varnames, paste0(varnames, "_Eq"), c("Dist", "Sign", "Posn", "Dist_Arg1", "Dist_Arg2", "Dist_Arg3", "LR"))
      dimnames(results$pdetA) <- list("detA", "detA", c("Dist", "Sign", "Posn", "Dist_Arg1", "Dist_Arg2", "Dist_Arg3"))
      dimnames(results$pH) <- list(varnames, paste0(varnames, "_Eq"), c("Dist", "Sign", "Posn", "Dist_Arg1", "Dist_Arg2", "Dist_Arg3"))
      
      dimnames(results$pP) <- list(colnames(x1), paste0(varnames, "_Eq"))
      dimnames(results$pP_sig) <- list(colnames(x1), colnames(x1))
      
      dimnames(results$pR) <- list(colnames(x1), paste0(varnames, "_Eq"), paste0(varnames, "_Eq"))
      dimnames(results$pR_sig) <- list(colnames(x1), colnames(x1), paste0(varnames, "_Eq"))
      
      dimnames(results$tau1) <- list(varnames, varnames)
      dimnames(results$kappa1) <- list(varnames, varnames)
      
      dimnames(results$A_start) <- list(varnames, paste0(varnames, "_Eq"))
      
      dimnames(results$A) <- list(varnames, paste0(varnames, "_Eq"), paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
      dimnames(results$detA) <- list("detA", "detA", paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
      dimnames(results$H) <- list(varnames, paste0(varnames, "_Eq"), paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
      
      dimnames(results$B) <- list(colnames(x1), paste0(colnames(y1), "_Eq"), paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
      dimnames(results$Phi) <- list(colnames(x1), paste0(colnames(y1), "_Eq"), paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
      dimnames(results$D) <- list(varnames, paste0(varnames, "_Eq"), paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
      
      dimnames(results$A_den$hori) <- list(varnames, paste0(varnames, "_Eq"), NULL)
      dimnames(results$A_den$vert) <- list(varnames, paste0(varnames, "_Eq"), NULL)
      dimnames(results$detA_den$hori) <- list("detA", "detA", NULL)
      dimnames(results$detA_den$vert) <- list("detA", "detA", NULL)
      dimnames(results$H_den$hori) <- list(varnames, paste0(varnames, "_Eq"), NULL)
      dimnames(results$H_den$vert) <- list(varnames, paste0(varnames, "_Eq"), NULL)
      
      Line_Plots(raw_array = results$A_chain, priors_array = results$pA, prior_name = "pA")
      ACF_Plots(raw_array = results$A_chain, priors_array = results$pA, prior_name = "pA")
      
      return(results)
    },
    error = function(e) {e}
  )
  
  return(results)
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
#' @param A_titles \emph{(n x n)} matrix containing the titles for the plots of the estimated parameters in the coefficient matrix \emph{A}. \emph{n} is the number of endogenous variables.
#' @param H_titles \emph{(n x n)} matrix containing the titles for the plots of the estimated parameters in the coefficient matrix \emph{H} (default = NULL). \emph{n} is the number of endogenous variables.
#' @param xlab Character label for the horizontal axis of historical decomposition plots (default = NULL). Default produces plots without a label for the horizontal axis.
#' @param ylab Character label for the vertical axis of historical decomposition plots (default = NULL). Default produces plots without a label for the vertical axis.
#' @details Plots posterior distributions against prior distributions.
#' @examples
#' # Import data
#' library(BHSBVAR)
#' set.seed(123)
#' data(USLMData)
#' y0 <- matrix(data = c(USLMData$Wage, USLMData$Employment), ncol = 2)
#' y <- y0 - (matrix(data = 1, nrow = nrow(y0), ncol = ncol(y0)) %*% 
#'              diag(x = colMeans(x = y0, na.rm = FALSE, dims = 1)))
#' colnames(y) <- c("Wage", "Employment")
#' 
#' # Set function arguments
#' nlags <- 8
#' itr <- 5000
#' burn <- 0
#' thin <- 1
#' acc <- TRUE
#' h <- 20
#' cri <- 0.95
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
#' pA[, , 8] <- c(2, NA, 2, NA)
#' 
#' # Position priors for Phi
#' pP <- matrix(data = 0, nrow = ((nlags * ncol(pA)) + 1), ncol = ncol(pA))
#' pP[1:nrow(pA), 1:ncol(pA)] <-
#'   diag(x = 1, nrow = nrow(pA), ncol = ncol(pA))
#' 
#' # Confidence in the priors for Phi
#' x1 <- 
#'   matrix(data = NA, nrow = (nrow(y) - nlags), 
#'          ncol = (ncol(y) * nlags))
#' for (k in 1:nlags) {
#'   x1[, (ncol(y) * (k - 1) + 1):(ncol(y) * k)] <-
#'     y[(nlags - k + 1):(nrow(y) - k),]
#' }
#' x1 <- cbind(x1, 1)
#' colnames(x1) <- 
#'   c(paste(rep(colnames(y), nlags),
#'           "_L",
#'           sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(y)),
#'                decreasing = FALSE),
#'           sep = ""),
#'     "cons")
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
#' # Confidence in long-run restriction priors
#' pR_sig <-
#'   array(data = 0,
#'         dim = c(((nlags * ncol(y)) + 1),
#'                 ((nlags * ncol(y)) + 1),
#'                 ncol(y)))
#' Ri <-
#'   cbind(kronecker(matrix(data = 1, nrow = 1, ncol = nlags),
#'                   matrix(data = c(1, 0), nrow = 1)),
#'         0)
#' pR_sig[, , 2] <- (t(Ri) %*% Ri) / 0.1
#' 
#' # Confidence in priors for D
#' kappa1 <- matrix(data = 2, nrow = 1, ncol = ncol(y))
#' 
#' # Set graphical parameters
#' par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
#'     mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
#' 
#' # Estimate the parameters of the model
#' results1 <- 
#'   BH_SBVAR(y = y, nlags = nlags, pA = pA, pP = pP, pP_sig = pP_sig,
#'            pR_sig = pR_sig, kappa1 = kappa1, itr = itr, burn = burn,
#'            thin = thin, cri = cri)
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
  test <- plot_funs_args_check(results = results, xlab = xlab, ylab = ylab)
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
      if ((is.na(pA[i, j, 1])) && (!is.na(A_titles[i, j]))) {
        stop(paste("A_titles: A_titles[", i, ", ", j, "] should be empty since pA[", i, ", ", j, ", ", 1, "] is empty.", sep = ""))
      }
      if ((!is.na(pA[i, j, 1])) && (is.na(A_titles[i, j]))) {
        stop(paste("A_titles: A_titles[", i, ", ", j, "] is missing.", sep = ""))
      }
      if ((is.na(pH[i, j, 1])) && (!is.na(H_titles[i, j]))) {
        stop(paste("H_titles: H_titles[", i, ", ", j, "] should be empty since pH[", i, ", ", j, ", ", 1, "] is empty.", sep = ""))
      }
      if ((!is.na(pH[i, j, 1])) && (is.na(H_titles[i, j]))) {
        stop(paste("H_titles: H_titles[", i, ", ", j, "] is missing.", sep = ""))
      }
    }
  }
  
  nticks0 <- 3
  
  list1 <- list(A_den = A_den, H_den = H_den)
  list2 <- list(pA = pA, pH = pH)
  hori <- matrix(data = NA_real_, nrow = dim(A_den$hori)[3], ncol = 1)
  vert <- matrix(data = NA_real_, nrow = dim(A_den$vert)[3], ncol = 1)
  den1 <- matrix(data = NA_real_, nrow = dim(A_den$vert)[3], ncol = 2)
  prior_den <- matrix(data = NA_real_, nrow = 501, ncol = 2)
  max_distance <- matrix(data = 0, nrow = dim(pA)[1], ncol = 1)
  crr_lb <- 0.05
  for (k in 1:length(list1)) {
    if (names(list2[k]) == "pA") {
      elast <- -1
    } else {
      elast <- 1
    }
    hori[, ] <- NA_real_
    vert[, ] <- NA_real_
    den1[, ] <- NA_real_
    prior_den[, ] <- NA_real_
    max_distance[, ] <- 0
    distance <- 0
    for (i in 1:(dim(list2[[k]])[1])) { #equations are by column
      for (j in 1:(dim(list2[[k]])[2])) {
        hori[, ] <- list1[[k]]$hori[i, j, ]
        vert[, ] <- list1[[k]]$vert[i, j, ]
        if (any(!is.na(hori))) {
          distance <- ceiling(max(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1], na.rm = TRUE) - min(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1], na.rm = TRUE)) + 1
        } else {
          distance <- 0
        }
        if (distance > max_distance[i, 1]) {
          max_distance[i, 1] <- distance
        }
      }
    }
    
    for (i in 1:(dim(list2[[k]])[1])) { #equations are by column
      for (j in 1:(dim(list2[[k]])[2])) {
        hori[, ] <- list1[[k]]$hori[i, j, ]
        vert[, ] <- list1[[k]]$vert[i, j, ]
        if (!is.na(list2[[k]][i, j, 1])) {
          if (is.na(list2[[k]][i, j, 2])) {
            ub <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) + (max_distance[i, 1] * 0.5)
            lb <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) - (max_distance[i, 1] * 0.5)
          } else if (list2[[k]][i, j, 2] == 1) {
            if ((round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0) - (max_distance[i, 1] * 0.5)) < 1) {
              if (names(list2[k]) == "pA") {
                ub <- 0
                lb <- (-1) * max_distance[i, 1]
              } else {
                ub <- max_distance[i, 1]
                lb <- 0
              }
            } else {
              ub <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) + (max_distance[i, 1] * 0.5)
              lb <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) - (max_distance[i, 1] * 0.5)
            }
          } else if (list2[[k]][i,j,2] == (-1)) {
            if ((round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0) + (max_distance[i, 1] * 0.5)) > (-1)) {
              if (names(list2[k]) == "pA") {
                ub <- max_distance[i, 1]
                lb <- 0
              } else {
                ub <- 0
                lb <- (-1) * max_distance[i, 1]
              }
            } else {
              ub <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) + (max_distance[i, 1] * 0.5)
              lb <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) - (max_distance[i, 1] * 0.5)
            }
          }
          den1[, ] <- cbind(hori[, 1], vert[, 1])
          den_plot(list2 = list2, den1 = den1, elast = elast, lb = lb, ub = ub, nticks0 = nticks0, A_titles = A_titles, H_titles = H_titles, xlab = xlab, ylab = ylab, k = k, j = j, i = i)
          
          prior_den[, 1] <- seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500))
          
          if (list2[[k]][i, j, 1] == 0) {
            if (is.na(list2[[k]][i, j, 2])) {
              
              for (h in 1:nrow(prior_den)) {
                prior_den[h, 2] <- prior_t(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              }

            } else if (list2[[k]][i, j, 2] == 1) {
              
              for (h in 1:nrow(prior_den)) {
                prior_den[h, 2] <- prior_t_p(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              }

            } else if (list2[[k]][i, j, 2] == (-1)) {
              
              for (h in 1:nrow(prior_den)) {
                prior_den[h, 2] <- prior_t_n(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              }

            }
          } else if (list2[[k]][i, j, 1] == 1) {
            
            for (h in 1:nrow(prior_den)) {
              prior_den[h, 2] <- prior_nonc_t(prior_den[h, 1], list2[[k]][i, j, 3], list2[[k]][i, j, 4], list2[[k]][i, j, 5], list2[[k]][i, j, 6])
            }

          } else if (list2[[k]][i, j, 1] == 2) {
            
            for (h in 1:nrow(prior_den)) {
              if (((list2[[k]][i, j, 2] == 1) & (prior_den[h, 1] >= 1)) | ((list2[[k]][i, j, 2] == -1) & (prior_den[h, 1] <= -1))) {
                prior_den[h, 2] <- prior_ibeta(abs(prior_den[h, 1]), list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              } else {
                prior_den[h, 2] <- 0
              }
            }
            
          } else if (list2[[k]][i, j, 1] == 3) {
            
            for (h in 1:nrow(prior_den)) {
              
              if (((list2[[k]][i, j, 2] == 1) & (prior_den[h, 1] >= 0) & (prior_den[h, 1] <= 1)) | ((list2[[k]][i, j, 2] == -1) & (prior_den[h, 1] <= 0) & (prior_den[h, 1] >= -1))) {
                prior_den[h, 2] <- prior_beta(abs(prior_den[h, 1]), list2[[k]][i, j, 4], list2[[k]][i, j, 5])
              } else {
                prior_den[h, 2] <- 0
              }
            }
            
          }
          
          graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
          
        }
      }
    }
  }
  
  list2 <- list(pdetA = pdetA)
  hori[, ] <- NA_real_
  vert[, ] <- NA_real_
  den1[, ] <- NA_real_
  prior_den[, ] <- NA_real_
  max_distance[, ] <- 0
  elast <- 1
  if (!is.na(pdetA[1, 1, 1])) {
    hori[, ] <- detA_den$hori[1, 1, ]
    vert[, ] <- detA_den$vert[1, 1, ]
    max_distance[1 ,1] <- ceiling(max(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1], na.rm = TRUE) - min(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1], na.rm = TRUE)) + 1
    
    if (is.na(list2[[1]][1, 1, 2])) {
      ub <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) + (max_distance[1 ,1] * 0.5)
      lb <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) - (max_distance[1 ,1] * 0.5)
    } else if (list2[[1]][1, 1, 2] == 1) {
      if ((round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0) - (max_distance[1 ,1] * 0.5)) < 0) {
        ub <- max_distance[1 ,1]
        lb <- 0
      } else {
        ub <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) + (max_distance[1 ,1] * 0.5)
        lb <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) - (max_distance[1 ,1] * 0.5)
      }
    } else if (list2[[1]][1, 1, 2] == (-1)) {
      if ((round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0) + (max_distance[1 ,1] * 0.5)) > 0) {
        ub <- 0
        lb <- (-1) * max_distance[1 ,1]
      } else {
        ub <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) + (max_distance[1 ,1] * 0.5)
        lb <- (elast * round(x = mean(hori[which(vert > (max(vert, na.rm = TRUE) * crr_lb)), 1]), digits = 0)) - (max_distance[1 ,1] * 0.5)
      }
    }
    
    den1[, ] <- cbind(hori[, 1], vert[, 1])
    den_plot(list2 = list2, den1 = den1, elast = elast, lb = lb, ub = ub, nticks0 = nticks0, A_titles = A_titles, H_titles = H_titles, xlab = xlab, ylab = ylab, k = 1, j = 1, i = 1)
    
    prior_den[, 1] <- seq(from = (elast * lb), to = (elast * ub), by = (elast * (ub - lb) / 500))
    
    if (pdetA[1, 1, 1] == 0) {
      if (is.na(pdetA[1, 1, 2])) {
        
        for (h in 1:nrow(prior_den)) {
          prior_den[h, 2] <- prior_t(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        }

      } else if (pdetA[1, 1, 2] == 1) {
        
        for (h in 1:nrow(prior_den)) {
          prior_den[h, 2] <- prior_t_p(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        }

      } else if (pdetA[1, 1, 2] == (-1)) {
        
        for (h in 1:nrow(prior_den)) {
          prior_den[h, 2] <- prior_t_n(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        }

      }
    } else if (pdetA[1, 1, 1] == 1) {
      
      for (h in 1:nrow(prior_den)) {
        prior_den[h, 2] <- prior_nonc_t(prior_den[h, 1], list2[[1]][1, 1, 3], list2[[1]][1, 1, 4], list2[[1]][1, 1, 5], list2[[1]][1, 1, 6])
      }

    } else if (pdetA[1, 1, 1] == 2) {
      for (h in 1:nrow(prior_den)) {
        if (((list2[[1]][1, 1, 2] == 1) & (prior_den[h, 1] >= 1)) | ((list2[[1]][1, 1, 2] == -1) & (prior_den[h, 1] <= -1))) {
          prior_den[h, 2] <- prior_ibeta(abs(prior_den[h, 1]), list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        } else {
          prior_den[h, 2] <- 0
        }
      }
    } else if (pdetA[1, 1, 1] == 3) {
            
      for (h in 1:nrow(prior_den)) {
        if (((list2[[1]][1, 1, 2] == 1) & (prior_den[h, 1] >= 0) & (prior_den[h, 1] <= 1)) | ((list2[[1]][1, 1, 2] == -1) & (prior_den[h, 1] <= 0) & (prior_den[h, 1] >= -1))) {
          prior_den[h, 2] <- prior_beta(abs(prior_den[h, 1]), list2[[1]][1, 1, 4], list2[[1]][1, 1, 5])
        } else {
          prior_den[h, 2] <- 0
        }
      }
            
    }
    
    graphics::lines(x = (elast * prior_den[, 1]), y = prior_den[, 2], type = "l", col = "red")
    
  }
}

# Check results from BH_SBVAR.
#' @keywords internal
BH_SBVAR_results_check <- function(results) {
  if (!is.list(results)) {
    return("results: Must be a list.")
  }
  if (!is.matrix(results$y)) {
    return("results: Must be a list containing a matrix 'y'.")
  }
  if (!is.array(results$A_chain) | !is.array(results$B_chain) | !is.array(results$D_chain)) {
    return("results: Must be a list containing an array 'A_chain', 'B_chain', and 'D_chain'.")
  }
  if ((length(dim(results$A_chain)) != 3) | (length(dim(results$B_chain)) != 3) | (length(dim(results$D_chain)) != 3)) {
    return("results: 'A_chain', 'B_chain', and 'D_chain' must have 3 dimensions.")
  }
  if ((dim(results$A_chain)[3] != dim(results$B_chain)[3]) | (dim(results$A_chain)[3] != dim(results$D_chain)[3])) {
    return("results: The number of slices of the third dimension of 'A_chain' must equal the number of slices of the third dimension of 'B_chain' (or 'D_chain').")
  }
  if ((ncol(results$y) != ncol(results$A_chain)) | (ncol(results$y) != nrow(results$A_chain))) {
    return("results: The number of columns and the number of rows of 'A_chain' must be equal to the number of columns in 'y'.")
  }
  if ((ncol(results$y) != ncol(results$D_chain)) | (ncol(results$y) != nrow(results$D_chain))) {
    return("results: The number of columns (or rows) of 'D_chain' must be equal to the number of columns in 'y'.")
  }
  if ((ncol(results$y) != ncol(results$D_chain)) | (ncol(results$y) != nrow(results$D_chain))) {
    return("results: The number of columns (or rows) of 'D_chain' must be equal to the number of columns in 'y'.")
  }
  if (all(!is.numeric(results$nlags)) || (length(results$nlags) > 1)) {
    return("results: 'nlags' should be a single number.")
  }
  if ((results$nlags != round(x = results$nlags, digits = 0)) | (results$nlags < 2)) {
    return("results: 'nlags' should be a positive integer greater than 2.")
  }
  if ((ncol(results$y) != ncol(results$B_chain)) || (((results$nlags * ncol(results$y)) + 1) != nrow(results$B_chain))) {
    return(paste0("results: 'B_chain' should have ", ((results$nlags * ncol(results$y)) + 1), " ((results$nlags * ncol(results$y)) + 1) rows and ", ncol(results$y), " (ncol(results$y)) columns."))
  }
  return("pass")
}

# Check arguments from the IRF_Plots, HD_Plots, Dist_Plots functions.
#' @keywords internal
plot_funs_args_check <- function(results, xlab, ylab) {
  if (is.array(results) && (length(dim(results)) == 3)) {
    if (length(results) == 0) {
      return(paste("results: Must be an array obtained from running IRF(), HD(), or FEVD().", sep = ""))
    }
    if ((any(!is.finite(results))) || (dim(results)[1] < 4) || (dim(results)[3] != 3) || ((dim(results)[2] / sqrt(dim(results)[2])) != (sqrt(dim(results)[2])))) {
      return(paste("results: Results from IRF(), or FEVD() are not present.", sep = ""))
    }
  } else if (is.list(results)) {
    if (any(names(results) == "HD")) {
      if (!is.array(results$HD) || (length(results$HD) == 0)) {
        return(paste("results: Must be an array obtained from running HD().", sep = ""))
      }
      if ((any(!is.finite(results$HD))) || (dim(results$HD)[1] < 4) || (dim(results$HD)[3] != 3) || ((dim(results$HD)[2] / sqrt(dim(results$HD)[2])) != (sqrt(dim(results$HD)[2])))) {
        return(paste("results: Results from HD() are not present.", sep = ""))
      }
    } else {
      if ((is.null(results$A)) || (!is.array(results$A)) || (any(!is.finite(results$A))) || (dim(results$A)[1] != dim(results$A)[2]) || (dim(results$A)[3] != 3) || (dim(results$A)[2] != ncol(results$y)) || (dim(results$A)[2] < 2)) {
        return(paste("results: A from BH_SBVAR() is not present.", sep = ""))
      }
    }
  } else {
    return(paste0("results: Results must be an array or a list containing an array"))
  }
  
  if ((!is.null(xlab)) && ((!is.character(xlab)) || (length(xlab) != 1))) {
    return(paste("xlab: Must be a character vector containing the label for the horizontal axis.", sep = ""))
  }
  if ((!is.null(ylab)) && ((!is.character(ylab)) || (length(ylab) != 1))) {
    return(paste("ylab: Must be a character vector containing the label for the vertical axis.", sep = ""))
  }
  return("pass")
}


