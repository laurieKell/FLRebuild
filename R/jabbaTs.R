#' Extract JABBA Time Series as FLQuants
#'
#' @description
#' Converts selected slices of a JABBA \code{fit$timeseries} array into an
#' \code{FLQuants} object for FLR plotting and process-error diagnostics.
#'
#' @param fit A JABBA fit with a 3-d \code{timeseries} array
#'   (year x statistic x quantity). Also accepts wrappers with
#'   \code{fit$fit$timeseries}.
#' @param quant Character. Second-dimension statistic. Default \code{"mu"}.
#' @param vars Named character vector mapping output names to JABBA quantities.
#'   Default maps \code{B}, \code{F}, \code{BBmsy}, \code{FFmsy}, \code{procB},
#'   \code{SPt} to \code{stock}, \code{harvest}, \code{bbmsy}, \code{ffmsy},
#'   \code{pe}, \code{sprod}. Missing names are skipped with a warning.
#'
#' @return An \code{FLQuants} of annual series.
#'
#' @details
#' JABBA stores process error as \code{procB}. This exposes that series
#' directly; see \code{\link{jabbaPE}} for production-function residuals.
#'
#' @examples
#' \dontrun{
#' fit <- runJABBA(stk, method = "ices", quick = TRUE)
#' ts  <- jabbaTs(fit)
#' plotPe(ts)
#' }
#'
#' @seealso \code{\link{plotPe}}, \code{\link{jabbaPE}}, \code{\link{rod}}
#' @export
#' @importFrom FLCore FLQuants as.FLQuant
jabbaTs <- function(fit,
                    quant = "mu",
                    vars = c(stock = "B",
                             harvest = "F",
                             bbmsy = "BBmsy",
                             ffmsy = "FFmsy",
                             pe = "procB",
                             sprod = "SPt")) {

  if (is.null(fit))
    stop("'fit' cannot be NULL")

  if (is.null(fit$timeseries) && !is.null(fit$fit$timeseries))
    fit <- fit$fit

  ts <- fit$timeseries
  if (is.null(ts))
    stop("JABBA fit object must have a 'timeseries' component")
  if (length(dim(ts)) != 3)
    stop("'fit$timeseries' must be a 3-d array (year x statistic x quantity)")

  dn <- dimnames(ts)
  years <- dn[[1]]
  stats <- dn[[2]]
  quants <- dn[[3]]

  if (!quant %in% stats)
    stop("Statistic '", quant, "' not found in fit$timeseries. Available: ",
         paste(stats, collapse = ", "))

  out <- list()
  for (nm in names(vars)) {
    jnm <- vars[[nm]]
    if (!jnm %in% quants) {
      warning("Quantity '", jnm, "' (as '", nm,
              "') not found in fit$timeseries; skipping")
      next
    }
    out[[nm]] <- FLCore::as.FLQuant(
      data.frame(year = years, data = as.numeric(ts[, quant, jnm]))
    )
  }

  if (!length(out))
    stop("No requested timeseries quantities were found in fit$timeseries")

  do.call(FLCore::FLQuants, out)
}
