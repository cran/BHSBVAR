#' Forecast Error Variance Decompositions
#' 
#' Forecast Error Variance Decompositions
#' @author Paul Richardson
#' @export
#' @import Rcpp
#' @name FEVD
#' @param results List containing the results from running BH_SBVAR().
#' @param h Integer specifying the time horizon for computing impulse responses (default = 12).
#' @param acc Boolean indicating whether accumulated impulse responses are to be returned (default = TRUE).
#' @param cri credibility intervals for the estimates to be returned (default = 0.95). A value of 0.95 will return 95\% credibility intervals. A value of 0.90 will return 90\% credibility intervals.
#' @details Computes forecast error variance decomposition estimates.
#' @return An array containing forecast error variance decomposition estimates.
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
#' fevd <- FEVD(results = results1, h = h, acc = acc, cri = cri)
#' 
#' # Plot impulse responses
#' varnames <- colnames(USLMData)[2:3]
#' shocknames <- c("Labor Demand","Labor Supply")
#' fevd_results <- 
#'   FEVD_Plots(results = fevd, varnames = varnames,
#'             shocknames = shocknames)
FEVD <- function(results, h = 12, acc = TRUE, cri = 0.95) {
  
  test <- BH_SBVAR_results_check(results = results)
  if (test != "pass") {
    stop(test)
  }
  
  if ((all(!is.numeric(h))) || (length(h) > 1)) {
    stop("h: Should be a single number.")
  }
  if ((h != round(x = h, digits = 0)) | (h < 3)) {
    stop("h: Should be a positive integer greater than 3.")
  }
  if ((all(!is.numeric(cri))) || (length(cri) > 1)) {
    "cri: Should be a single number."
  }
  if ((cri > 1) | (cri < 0.5)) {
    stop("cri: Should be a positive value between 1 and 0.5.")
  }
  if ((!is.logical(acc)) || (is.na(acc))) {
    stop(paste("acc_irf: Must be logical 'TRUE' or 'FALSE'.", sep = ""))
  }
  
  varnames <- colnames(results$y)
  ci <- (1.0 - ((1.0 - cri) / 2.0))
  nvar <- dim(results$A_chain[, , ])[2]
  nlags <- results$nlags
  nsli <- dim(results$A_chain[, , ])[3]
  
  fevd <- fevd_estimates(A_chain = results$A_chain[, , ], B_chain = results$B_chain[, , ], D_chain = results$D_chain[, , ], nlags = nlags, h = h, acc = acc, ci = ci)
  
  dimnames(fevd) <- list(NULL, paste0("Res_", varnames, paste0("_Shk_", varnames[c(sort(x = rep(x = c(1:nvar), times = nvar)))])), paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
  
  return(fevd)
}

#' Plot Forecast Error Variance Decompositions
#' 
#' Plot Forecast Error Variance Decompositions
#' @author Paul Richardson
#' @export
#' @name FEVD_Plots
#' @param results List containing the results from running BH_SBVAR().
#' @param varnames Character vector containing the names of the endogenous variables.
#' @param shocknames Character vector containing the names of the shocks.
#' @param xlab Character label for the horizontal axis of impulse response plots (default = NULL). Default produces plots without a label for the horizontal axis.
#' @param ylab Character label for the vertical axis of impulse response plots (default = NULL). Default produces plots without a label for the vertical axis.
#' @param rel Boolean indicating whether to display forecast error variance explained by the shock as a percent of total forecast error variance (default = TRUE).
#' @details  Plots forecast error variance decompositions and returns a list containing the actual processed data used to create the plots.
#' @return A list containing forecast error variance decompositions.
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
#' fevd <- FEVD(results = results1, h = h, acc = acc, cri = cri)
#' 
#' # Plot impulse responses
#' varnames <- colnames(USLMData)[2:3]
#' shocknames <- c("Labor Demand","Labor Supply")
#' fevd_results <- 
#'   FEVD_Plots(results = fevd, varnames = varnames,
#'             shocknames = shocknames)
FEVD_Plots <- function(results, varnames, shocknames = NULL, xlab = NULL, ylab = NULL, rel = TRUE) {
  
  #test arguments
  test <- plot_funs_args_check(results = results, xlab = xlab, ylab = ylab)
  if (test != "pass") {
    stop(test)
  }
  if (!is.vector(varnames) || (!is.character(varnames)) || (length(varnames) != sqrt(dim(results)[2]))) {
    stop(paste("varnames: Must be a character vector containing the names of the endogenous variables.", sep = ""))
  }
  if (is.null(shocknames)) {
    shocknames <- varnames
  }
  if (!is.vector(shocknames) || (!is.character(shocknames)) || (length(shocknames) != sqrt(dim(results)[2]))) {
    stop(paste("shocknames: Must be a character vector containing the names of the shocks.", sep = ""))
  }
  
  if (is.null(xlab)) {
    xlab <- ""
  }
  if (is.null(ylab)) {
    ylab <- ""
  }
  
  if (!isTRUE(rel) & !isFALSE(rel)) {
    stop("rel: Must be TRUE or FALSE.")
  }
  
  fevd <- results
  nvar <- round(sqrt(dim(results)[2]))
  xticks <- floor(dim(fevd)[1] / 4)
  fevdcolnames <- dimnames(fevd)[[2]]
  
  colors1 <- c("black", "red", "black", "red")
  lty1 <- c(1, 2, 1, 2)
  if (rel) {
    colors1[] <- rep(x = "black", times = length(colors1))
    lty1[] <- rep(x = 1, times = length(lty1))
    ind1 <- rep(x = 0, times = sqrt(dim(fevd)[2]))
    for (i in 1:sqrt(dim(fevd)[2])) {
      ind1[] <- seq(from = i, to = dim(fevd)[2], by = sqrt(dim(fevd)[2]))
      fevd[, ind1, "50%"] <- fevd[, ind1, "50%"] * matrix(data = rep(x = (100 / rowSums(fevd[, ind1, "50%"])), times = nvar), ncol = nvar, dimnames = list(NULL, fevdcolnames[ind1]))
    }
    fevd[, , c(1, 3)] <- fevd[, , "50%"]
  }
  dimnames(fevd)[[3]] <- rep(x = "50%", times = 3)
  
  #store results from fevd responses
  fevd_results <- vector(mode = "list", length = (nvar * nvar))
  for (j in 1:nvar) {
    for (i in 1:nvar) {
      #fevd responses
      names(fevd_results)[((nvar * (j - 1)) + i)] <- fevdcolnames[((nvar * (j - 1)) + i)]
      fevd_results[[((nvar * (j - 1)) + i)]] <- fevd[, ((nvar * (j - 1)) + i), ]
      #fevd response plots
      mat_ts <- stats::ts(cbind(0, fevd_results[[((nvar * (j - 1)) + i)]]))
      colnames(mat_ts) <- c("Series1", "Series2", "Series3", "Series4")
      stats::ts.plot(mat_ts, col = colors1, gpars = list(xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "r", xaxt = "n", lty = lty1))
      graphics::title(main = paste("Response of ", varnames[i], " to ", shocknames[j], sep = ""), col.main = "black")
      graphics::axis(side = 1, at = seq(from = 1, to = nrow(mat_ts), by = xticks), labels = seq(from = 0, to = (nrow(mat_ts) - 1),by = xticks))
    }
  }
  return(fevd_results)
}
