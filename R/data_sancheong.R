#' Sancheong Station Annual r-Largest Daily Rainfall Data
#'
#' Annual r-largest daily rainfall observations recorded at Sancheong
#' (산청) meteorological station, South Korea. The data consist of the
#' 20 largest daily rainfall amounts per year, used for extreme value
#' analysis with r-largest order statistics models.
#'
#' @format A data frame with 51 rows and 21 columns:
#' \describe{
#'   \item{year}{Observation year (1972--2022).}
#'   \item{X1}{Annual maximum daily rainfall (mm) — 1st largest.}
#'   \item{X2}{2nd largest daily rainfall (mm).}
#'   \item{X3}{3rd largest daily rainfall (mm).}
#'   \item{...}{...}
#'   \item{X20}{20th largest daily rainfall (mm).}
#' }
#' Note: 1972 has only the annual maximum (X1); X2--X20 are \code{NA}
#' for that year.
#'
#' @details
#' Station: Sancheong (산청), South Korea (ASOS 285)\cr
#' Period : 1972--2022 (51 years)\cr
#' Unit   : millimetres (mm)\cr
#' Source : Korea Meteorological Administration (KMA)\cr
#'
#' The record maximum is 332.5 mm, observed in 1998.
#' The dataset is intended for use with \code{\link{rsel.rgev}} and
#' \code{\link{rsel.rgev11}} to select the optimal number r of
#' r-largest order statistics.
#'
#' @examples
#' data(sancheong)
#' head(sancheong)
#'
#' # Extract the r-largest matrix (exclude year column)
#' xdat <- as.matrix(sancheong[, -1])
#'
#' \dontrun{
#' # Select optimal r using energy-distance method
#' result <- rsel.rgev(xdat, method = "ed", sigL = 0.05)
#' cat("Selected r =", result$r.sel, "\n")
#' }
#'
#' @source Korea Meteorological Administration (KMA),
#'   \url{https://data.kma.go.kr}
"sancheong"
