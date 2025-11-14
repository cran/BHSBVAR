## ----Setup, include=TRUE, echo=FALSE------------------------------------------
knitr::opts_knit$set(
  concordance = TRUE
  )
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Data, include=TRUE, echo=TRUE, background='NA', highlight=FALSE----------
rm(list = ls())
library(BHSBVAR)
set.seed(123)
data(USLMData)
y0 <- matrix(data = c(USLMData$Wage, USLMData$Employment), ncol = 2)
y <- y0 - (matrix(data = 1, nrow = nrow(y0), ncol = ncol(y0)) %*% 
             diag(x = colMeans(x = y0, na.rm = FALSE, dims = 1)))
colnames(y) <- c("Wage", "Employment")

## ----Inputs, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
nlags <- 8
itr <- 200000
burn <- 0
thin <- 20
acc <- TRUE
h <- 20
cri <- 0.95

## ----pA, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
pA <- array(data = NA, dim = c(ncol(y), ncol(y), 8))
pA[, , 1] <- c(0, NA, 0, NA)
pA[, , 2] <- c(1, NA, -1, NA)
pA[, , 3] <- c(0.6, 1, -0.6, 1)
pA[, , 4] <- c(0.6, NA, 0.6, NA)
pA[, , 5] <- c(3, NA, 3, NA)
pA[, , 6] <- c(NA, NA, NA, NA)
pA[, , 7] <- c(NA, NA, 1, NA)
pA[, , 8] <- c(2, NA, 2, NA)

## ----pP_pP_sig, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
pP <- matrix(data = 0, nrow = ((nlags * ncol(pA)) + 1), ncol = ncol(pA))
pP[1:nrow(pA), 1:ncol(pA)] <-
  diag(x = 1, nrow = nrow(pA), ncol = ncol(pA))
x1 <- 
  matrix(data = NA, nrow = (nrow(y) - nlags), 
         ncol = (ncol(y) * nlags))
for (k in 1:nlags) {
  x1[, (ncol(y) * (k - 1) + 1):(ncol(y) * k)] <-
    y[(nlags - k + 1):(nrow(y) - k),]
}
x1 <- cbind(x1, 1)
colnames(x1) <- 
  c(paste(rep(colnames(y), nlags),
          "_L",
          sort(rep(seq(from = 1, to = nlags, by = 1), times = ncol(y)),
               decreasing = FALSE),
          sep = ""),
    "cons")
y1 <- y[(nlags + 1):nrow(y),]
ee <- matrix(data = NA, nrow = nrow(y1), ncol = ncol(y1))
for (i in 1:ncol(y1)) {
  xx <- cbind(x1[, seq(from = i, to = (ncol(x1) - 1), by = ncol(y1))], 1)
  yy <- matrix(data = y1[, i], ncol = 1)
  phi <- solve(t(xx) %*% xx, t(xx) %*% yy)
  ee[, i] <- yy - (xx %*% phi)
}
somega <- (t(ee) %*% ee) / nrow(ee)
lambda0 <- 0.2
lambda1 <- 1
lambda3 <- 100
v1 <- matrix(data = (1:nlags), nrow = nlags, ncol = 1)
v1 <- v1^((-2) * lambda1)
v2 <- matrix(data = diag(solve(diag(diag(somega)))), ncol = 1)
v3 <- kronecker(v1, v2)
v3 <- (lambda0^2) * rbind(v3, (lambda3^2))
v3 <- 1 / v3
pP_sig <- diag(x = c(v3), nrow = nrow(v3), ncol = nrow(v3))

## ----pR_pR_sig_kappa1, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
pR_sig <-
  array(data = 0,
        dim = c(((nlags * ncol(y)) + 1),
                ((nlags * ncol(y)) + 1),
                ncol(y)))
Ri <-
  cbind(kronecker(matrix(data = 1, nrow = 1, ncol = nlags),
                  matrix(data = c(1, 0), nrow = 1)),
        0)
pR_sig[, , 2] <- (t(Ri) %*% Ri) / 0.1
kappa1 <- matrix(data = 2, nrow = 1, ncol = ncol(y))

## ----model, include=TRUE, echo=TRUE, fig.keep='all', fig.show='asis', fig.width=7, fig.align='center', fig.asp=0.4, background='NA', highlight=FALSE, fig.cap="Line and Autocorrelation Diagnostic Plots of the Posterior Estimates"----
par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
    mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
results1 <- 
  BH_SBVAR(y = y, nlags = nlags, pA = pA, pP = pP, pP_sig = pP_sig,
           pR_sig = pR_sig, kappa1 = kappa1, itr = itr,
           burn = burn, thin = thin, cri = cri)

## ----irf-plots, include=TRUE, echo=TRUE, fig.keep='all', fig.show='asis', fig.width=7, fig.align='center', fig.asp=0.4, background='NA', highlight=FALSE, fig.cap="Impulse Responses"----
irf <- IRF(results = results1, h = h, acc = acc, cri = cri)
varnames <- colnames(USLMData)[2:3]
shocknames <- c("Labor Demand","Labor Supply")
par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
    mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
irf_results <- 
  IRF_Plots(results = irf, varnames = varnames,
            shocknames = shocknames)

## ----fevd-plots, include=TRUE, echo=TRUE, fig.keep='all', fig.show='asis', fig.width=7, fig.align='center', fig.asp=0.4, background='NA', highlight=FALSE, fig.cap="Forecast Error Variance Decompositions"----
fevd <- FEVD(results = results1, h = h, acc = acc, cri = cri)
varnames <- colnames(USLMData)[2:3]
shocknames <- c("Labor Demand","Labor Supply")
par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
    mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
fevd_results <- 
  FEVD_Plots(results = fevd, varnames = varnames,
            shocknames = shocknames)

## ----hd-plots, include=TRUE, echo=TRUE, fig.keep='all', fig.show='asis', fig.width=7, fig.align='center', fig.asp=0.4, background='NA', highlight=FALSE, fig.cap="Historical Decompositions"----
hd <- HD(results = results1, cri = cri)
varnames <- colnames(USLMData)[2:3]
shocknames <- c("Labor Demand","Labor Supply")
freq <- 4
start_date <- 
  c(floor(USLMData[(nlags + 1), 1]),
    (floor(((USLMData[(nlags + 1), 1] %% 1) * freq)) + 1))
par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
    mfrow = c(2, 2), mar = c(2, 2.2, 2, 1), las = 1)
hd_results <- 
  HD_Plots(results  = hd, varnames = varnames,
           shocknames = shocknames,
           freq = freq, start_date = start_date)

## ----dist-plots, include=TRUE, echo=TRUE, fig.keep='all', fig.show='asis', fig.width=7, fig.align='center', fig.asp=0.4, background='NA', highlight=FALSE, fig.cap="Posterior and Prior Distribution Plots"----
A_titles <- 
  matrix(data = NA_character_, nrow = dim(pA)[1], ncol = dim(pA)[2])
A_titles[1, 1] <- "Wage Elasticity of Labor Demand"
A_titles[1, 2] <- "Wage Elasticity of Labor Supply"
par(cex.axis = 0.8, cex.main = 1, font.main = 1, family = "serif",
    mfrow = c(1, 2), mar = c(2, 2.2, 2, 1), las = 1)
Dist_Plots(results = results1, A_titles = A_titles)

## ----Density1, eval = FALSE, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
# density <-
#   dt(x = ((a1 - p1) / sigma1), df = nu, ncp = 0, log = FALSE) / sigma1

## ----Density2, eval = FALSE, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
# density <-
#   dt(x = ((a1 - p1) / sigma1), df = nu, ncp = lam1, log = FALSE) / sigma1

## ----Density3, eval = FALSE, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
# density <-
#   dt(x = ((a1 - p1) / sigma1), df = nu, ncp = 0, log = FALSE) /
#   (sigma1 *
#      (1 - pt(q = ((-p1) / sigma1), df = nu, ncp = 0, lower.tail = TRUE,
#              log.p = FALSE)))

## ----Density5, eval = FALSE, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
# density <- dt(x = ((a1 - p1) / sigma1), df = nu, ncp = 0, log = FALSE) /
#   (sigma1 * pt(q = ((-p1) / sigma1), df = nu, ncp = 0, lower.tail = TRUE,
#                log.p = FALSE))

## ----Density6, eval = FALSE, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
# density <- 0
# if (a1 >= 1) {
#   density <- exp(
#     ((sh2 - 1) * log((a1 - 1))) +
#       (((-1) * (sh2 + sh1)) * log((1 + (a1 - 1)))) -
#       log(beta(sh2, sh1))
#     )
# }

## ----Density4, eval = FALSE, include = TRUE, echo = TRUE, background = 'NA', highlight = FALSE----
# density <- dbeta(x = a1, shape1 = sh1, shape2 = sh2, ncp = 0, log = FALSE)

