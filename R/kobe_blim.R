# Kobe status and implied B_lim(t*) from a moment grid.
# Kept separate from mvln.R so a paper can source this file without
# replacing getRuns / ssLoad in the search path.

#' Kobe status from one or more moment-grid cells
#'
#' Median \eqn{B/B_{\mathrm{MSY}}} and \eqn{F/F_{\mathrm{MSY}}} are
#' \code{exp(muB)} and \code{exp(muF)}; deterministic \eqn{F = 0}
#' (\code{muF = -Inf}) is reported as zero. Quadrant probabilities use the
#' same MVLN / posterior moments as \code{\link{makeK2SM}}.
#'
#' @param grid A \code{\link{kobeMomentGrid}} / \code{\link{kobeMomentGridJabba}}
#'   \code{data.frame}.
#' @param year Optional year filter.
#' @param catch Optional catch filter (when \code{grid} has \code{catch}).
#' @return A \code{data.frame} with medians and \code{pNotOverfished},
#'   \code{pNoOverfishing}, \code{pGreen}.
#' @seealso \code{\link{makeK2SM}}, \code{\link{kobeGreen}}
#' @export
kobeStatusAt <- function(grid, year = NULL, catch = NULL) {

  need <- c("year", "muB", "muF", "sigmaB", "sigmaF", "rhoLog")
  if (!is.data.frame(grid) || !all(need %in% names(grid)))
    stop("'grid' must be kobeMomentGrid / kobeMomentGridJabba output")

  g <- grid
  if (!is.null(year)) g <- g[g$year %in% year, , drop = FALSE]
  if (!is.null(catch)) {
    if (!"catch" %in% names(g))
      stop("'catch' was supplied but 'grid' has no catch column")
    g <- g[g$catch %in% catch, , drop = FALSE]
  }
  if (!nrow(g))
    stop("no moment-grid rows match year/catch")

  zeroF <- is.infinite(g$muF) & g$muF < 0 & is.finite(g$sigmaF) & g$sigmaF == 0
  pB <- ifelse(!is.finite(g$muB) | !is.finite(g$sigmaB) | g$sigmaB <= 0,
               NA_real_, 1 - stats::pnorm(0, g$muB, g$sigmaB))
  pF <- ifelse(zeroF, 1,
               ifelse(!is.finite(g$muF) | !is.finite(g$sigmaF) | g$sigmaF <= 0,
                      NA_real_, stats::pnorm(0, g$muF, g$sigmaF)))
  out <- data.frame(
    year = g$year,
    bbmsy = exp(g$muB),
    ffmsy = ifelse(zeroF, 0, exp(g$muF)),
    pNotOverfished = pB,
    pNoOverfishing = pF,
    pGreen = mapply(kobeGreen, g$muB, g$muF, g$sigmaB, g$sigmaF, g$rhoLog,
                    MoreArgs = list(quadrant = "green")),
    stringsAsFactors = FALSE)
  if ("catch" %in% names(g)) out$catch <- g$catch
  if ("run" %in% names(g))   out$run   <- g$run
  out
}

#' Implied \eqn{B_{\mathrm{lim}}(t^{*})} along a constant-catch median path
#'
#' Takes the median relative biomass \eqn{e^{\mu_B}} at one catch (typically
#' TAC \eqn{= 0}, i.e.\ \eqn{F = 0}), finds the last year the median remains
#' at or above \code{threshold} (default \eqn{B_{\mathrm{MSY}}}), and inverts:
#' \eqn{B_{\mathrm{lim}}(t^{*})/B_{\mathrm{MSY}}} is the median relative
#' biomass \code{tstar} years before that crossing. This is the trajectory
#' analogue of \code{\link{makeRebuildRefs}} when a depletion grid is not
#' available, as in Kobe projections of a single stock.
#'
#' @param grid A moment grid with \code{year}, \code{muB}, and (usually)
#'   \code{catch}.
#' @param catch Catch level to follow; default \code{0}.
#' @param tstar Numeric horizon(s) in years.
#' @param years Optional year subset; default all years in \code{grid}.
#' @param threshold Relative-biomass threshold; default \code{1} (\eqn{B_{\mathrm{MSY}}}).
#' @return A \code{data.frame} with one row per \code{tstar}: \code{yearHit},
#'   \code{T_years} from the first year in the series, \code{blim_bbmsy}.
#' @seealso \code{\link{kobeRebuildCurve}}, \code{\link{makeRebuildRefs}}
#' @export
kobeImpliedBlim <- function(grid, catch = 0, tstar, years = NULL,
                            threshold = 1) {

  if (missing(tstar) || !length(tstar) || any(!is.finite(tstar) | tstar < 0))
    stop("'tstar' must be non-negative finite year-horizon(s)")
  if (!is.data.frame(grid) || !all(c("year", "muB") %in% names(grid)))
    stop("'grid' must have columns year and muB")

  g <- grid
  if (!is.null(years)) g <- g[g$year %in% years, , drop = FALSE]
  if ("catch" %in% names(g))
    g <- g[g$catch == catch, , drop = FALSE]
  g <- g[order(g$year), , drop = FALSE]
  y <- g$year
  b <- exp(g$muB)
  ok <- is.finite(y) & is.finite(b)
  y <- y[ok]; b <- b[ok]

  empty <- function()
    data.frame(catch = catch, startYear = if (length(y)) min(y) else NA_real_,
               yearHit = NA_real_, T_years = NA_real_, attained = FALSE,
               tstar = tstar, blim_bbmsy = NA_real_,
               stringsAsFactors = FALSE)

  if (length(y) < 2L) return(empty())

  if (b[length(b)] < threshold) {
    out <- empty()
    out$startYear <- min(y)
    return(out)
  }

  cross <- numeric(0)
  for (i in seq_len(length(y) - 1L)) {
    a <- b[i] - threshold
    c1 <- b[i + 1L] - threshold
    if (a * c1 <= 0 && b[i] != b[i + 1L])
      cross <- c(cross,
                 y[i] + (y[i + 1L] - y[i]) * (threshold - b[i]) /
                   (b[i + 1L] - b[i]))
  }
  yearHit <- if (!length(cross)) min(y) else max(cross)
  start <- min(y)
  blim <- vapply(tstar, function(ts) {
    as.numeric(stats::approx(y, b, xout = yearHit - ts, rule = 1)$y)
  }, numeric(1))

  data.frame(catch = catch, startYear = start, yearHit = yearHit,
             T_years = yearHit - start, attained = TRUE,
             tstar = tstar, blim_bbmsy = blim,
             stringsAsFactors = FALSE)
}
