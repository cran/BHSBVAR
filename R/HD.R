#' Historical Decompositions
#' 
#' Historical Decompositions
#' @author Paul Richardson
#' @export
#' @import Rcpp
#' @name HD
#' @param results List containing the results from running BH_SBVAR().
#' @param cri credibility intervals for the estimates to be returned (default = 0.95). A value of 0.95 will return 95\% credibility intervals. A value of 0.90 will return 90\% credibility intervals.
#' @details Computes historical decomposition estimates.
#' @return An array containing historical decomposition estimates.
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
#' hd <- HD(results = results1, cri = cri)
#' 
#' # Plot historical decompositions
#' varnames <- colnames(USLMData)[2:3]
#' shocknames <- c("Labor Demand","Labor Supply")
#' freq <- 4
#' start_date <- 
#'   c(floor(USLMData[(nlags + 1), 1]),
#'     round(((USLMData[(nlags + 1), 1] %% 1) * freq), digits = 0))
#' hd_results <- 
#'   HD_Plots(results  = hd, varnames = varnames,
#'            shocknames = shocknames,
#'            freq = freq, start_date = start_date)
HD <- function(results, cri = 0.95) {
  
  test <- BH_SBVAR_results_check(results = results)
  if (test != "pass") {
    stop(test)
  }
  
  if ((all(!is.numeric(cri))) || (length(cri) > 1)) {
    "cri: Should be a single number."
  }
  if ((cri > 1) | (cri < 0.5)) {
    stop("cri: Should be a positive value between 1 and 0.5.")
  }
  
  varnames <- colnames(results$y)
  ci <- (1.0 - ((1.0 - cri) / 2.0))
  nvar <- dim(results$A_chain[, , ])[2]
  nlags <- results$nlags
  nsli <- dim(results$A_chain[, , ])[3]
  
  hd <- list("y" = NULL, "HD" = NULL)
  
  hd$HD <- hd_estimates(A_chain = results$A_chain, B_chain = results$B_chain, y1 = results$y, x1 = results$x, nlags = nlags, ci = ci)
  
  dimnames(hd$HD) <- list(NULL, paste0("Res_", varnames, paste0("_Shk_", varnames[c(sort(x = rep(x = c(1:nvar), times = nvar)))])), paste0(c(((1 - ci) * 100), 50, (ci * 100)),"%"))
  
  hd$y <- results$y
  
  return(hd)
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
#' @param start_date Numeric vector indicating the date of the first observation of the endogenous variables included in the model.
#' @details Plots historical decompositions and returns a list containing the actual processed data used to create the plots.
#' @return A list containing historical decompositions.
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
#' hd <- HD(results = results1, cri = cri)
#' 
#' # Plot historical decompositions
#' varnames <- colnames(USLMData)[2:3]
#' shocknames <- c("Labor Demand","Labor Supply")
#' freq <- 4
#' start_date <- 
#'   c(floor(USLMData[(nlags + 1), 1]),
#'     round(((USLMData[(nlags + 1), 1] %% 1) * freq), digits = 0))
#' hd_results <- 
#'   HD_Plots(results  = hd, varnames = varnames,
#'            shocknames = shocknames,
#'            freq = freq, start_date = start_date)
HD_Plots <- function(results, varnames, shocknames = NULL, xlab = NULL, ylab = NULL, freq, start_date) {
  
  #test arguments
  test <- plot_funs_args_check(results = results, xlab = xlab, ylab = ylab)
  if (test != "pass") {
    stop(test)
  }
  if (!is.vector(varnames) || (!is.character(varnames)) || (length(varnames) != sqrt(dim(results$HD)[2])) || (length(varnames) != ncol(results$y))) {
    stop(paste("varnames: Must be a character vector containing the names of the endogenous variables.", sep = ""))
  }
  if (is.null(shocknames)) {
    shocknames <- varnames
  }
  if (!is.vector(shocknames) || (!is.character(shocknames)) || (length(shocknames) != sqrt(dim(results$HD)[2])) || (length(shocknames) != ncol(results$y))) {
    stop(paste("shocknames: Must be a character vector containing the names of the shocks.", sep = ""))
  }
  if ((!is.numeric(freq)) || (!is.finite(freq)) || (length(freq) != 1) || ((freq %% 1) != 0) || (freq < 1)) {
    stop("freq: Must be a finite whole number grater than 0.")
  }
  if ((!is.numeric(start_date)) || (any(!is.finite(start_date))) || (length(start_date) != 2) || (any((start_date %% 1) != 0)) || (any(start_date < 0))) {
    stop("start_date: Must be a numeric vector containing finite whole numbers greater than or equal to 0.")
  }
  
  if (is.null(xlab)) {
    xlab <- ""
  }
  if (is.null(ylab)) {
    ylab <- ""
  }
  
  y <- results$y
  hd <- results$HD
  nvar <- dim(results$y)[2]
  
  #store results from histroical decompositions
  hd_results <- vector(mode = "list", length = (nvar * nvar))
  for (j in 1:nvar) {
    for (i in 1:nvar) {
      #historical decomposition
      names(hd_results)[((nvar * (j - 1)) + i)] <- dimnames(hd)[[2]][((nvar * (j - 1)) + i)]
      hd_results[[((nvar * (j - 1)) + i)]] <- hd[, ((nvar * (j - 1)) + i), ]
      #historical decomposition plots
      mat_ts <- stats::ts(cbind(0, y[, i], hd_results[[((nvar * (j - 1)) + i)]]), frequency = freq, start = start_date)
      colnames(mat_ts) <- c("Series1", "Series2", "Series3", "Series4", "Series5")
      stats::ts.plot(mat_ts, col = c("black", "black", "red", "red", "red"), gpars = list(xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "r", lty = c(1, 1, 2, 1, 2)))
      graphics::title(main = paste("Contribution of ", shocknames[j], " Shocks on ", varnames[i], sep = ""), col.main = "black")
    }
  }
  return(hd_results)
}
