
#' @docType data
#' @title U.S. Labor Market Data
#' @description Quarterly U.S. labor market time-series data. These data are the data used in Baumeister and Hamilton (2015).
#' @source Dr. Christiane Baumeister's website \href{https://sites.google.com/site/cjsbaumeister/}{https://sites.google.com/site/cjsbaumeister/}.
#' @source Dr. James D. Hamilton's website \href{https://econweb.ucsd.edu/~jhamilton/}{https://econweb.ucsd.edu/~jhamilton/}.
#' @usage data(USLMData)
#' @format Data frame object that includes "Date", "Wage", and "Employment" variables. These data are the percent change in U.S. real wage and employment and were created by taking the difference of the natural log of U.S. real wage and employment levels and multiplying by 100.
#' @keywords datasets
#' @references Baumeister, C., & Hamilton, J.D. (2015). Sign restrictions, structural vector autoregressions, and useful prior information. \emph{Econometrica}, 83(5), 1963-1999.
#' @examples
#' data(USLMData)
 "USLMData"