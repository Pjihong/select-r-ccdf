#' Sancheong Annual Maximum Daily Rainfall (r-Largest Order Statistics)
#'
#' Annual r-largest order statistics of daily rainfall recorded at
#' Sancheong weather station (산청 기상관측소), South Korea.
#' Data span 51 years (1972--2022) with up to 20 largest values per year.
#'
#' @format A data frame with 51 rows and 21 columns:
#' \describe{
#'   \item{year}{Observation year (1972--2022).}
#'   \item{X1}{Annual maximum daily rainfall (mm), 1st largest.}
#'   \item{X2}{2nd largest daily rainfall (mm).}
#'   \item{...}{...}
#'   \item{X20}{20th largest daily rainfall (mm).}
#' }
#' The 1972 record contains only the annual maximum (X1);
#' all other r-values for that year are \code{NA}.
#'
#' @source Korea Meteorological Administration (KMA).
#'
#' @examples
#' data(sancheong)
#' head(sancheong)
#' dim(sancheong)   # 51 x 21
#'
#' # Use with rsel.rgev (stationary model, r up to 10)
#' \dontrun{
#' xdat <- as.matrix(sancheong[, 2:11])   # columns X1~X10
#' result <- rsel.rgev(xdat, method = "ed", sigL = 0.05)
#' cat("Selected r =", result$r.sel, "\n")
#' }
#'
#' # Use with rsel.rgev11 (nonstationary model)
#' \dontrun{
#' xdat <- as.matrix(sancheong[, 2:11])
#' result <- rsel.rgev11(xdat, model = "rgev11",
#'                       method = "ed", sigL = 0.05)
#' cat("Selected r =", result$r.sel, "\n")
#' cat("MLE =", result$mle, "\n")
#' }
"sancheong"
