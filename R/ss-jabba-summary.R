.fitType <- function(x) {
  if (!is.null(x$derived_quants)) "ss3"
  else if (!is.null(x$refpts_posterior)) "jabba"
  else NA_character_
}

.byScenario <- function(object, fn, ...) {
  res <- lapply(names(object), function(s) {
    d <- fn(object[[s]], ...)
    if (nrow(d)) data.frame(scenario = s, d) else NULL
  })
  do.call(rbind, res)
}

.ptRow <- function(msy, fmsy, bmsy, b0, tb0, m = NULL, r = NULL) {
  shape <- bmsy / b0
  if (is.null(m)) m <- shape2p(shape) + 1
  if (is.null(r)) r <- r_from_fmsy_shape(fmsy, shape, m - 1)
  data.frame(msy = msy, fmsy = fmsy, bmsy = bmsy, b0 = b0, tb0 = tb0,
             shape = shape, m = m, r = r, td = log(2) / r)
}

#' @rdname refptsSummary
#' @details For SS3, FMSY is the annual apical F, BMSY and B0 are SSB
#'   (\code{SSB_Virgin}) in SS3 units and \code{tb0} is the virgin total
#'   biomass. SS3 has no intrinsic rate, so \code{r} is the Pella-Tomlinson
#'   equivalent from \code{\link{r_from_fmsy_shape}}, with the shape solved
#'   from BMSY/B0 (\eqn{r = m F_{MSY}}). For JABBA, values are posterior
#'   medians and \code{r} and \code{m} are the estimates.
#' @export
setMethod("refptsSummary", signature(object = "list"), function(object, ...) {
  type <- .fitType(object)
  if (is.na(type)) return(.byScenario(object, refptsSummary, ...))
  if (type == "ss3") {
    dq <- object$derived_quants
    v  <- function(k) if (k %in% rownames(dq)) dq[k, "Value"] else NA_real_
    ts <- object$timeseries
    return(.ptRow(v("Dead_Catch_MSY"), v("annF_MSY"), v("SSB_MSY"), v("SSB_Virgin"),
                  ts$Bio_all[ts$Era == "VIRG"][1]))
  }
  q <- sapply(object$refpts_posterior, stats::median)
  p <- sapply(object$pars_posterior, stats::median)
  .ptRow(q[["MSY"]], q[["Fmsy"]], q[["Bmsy"]], p[["K"]], p[["K"]], m = p[["m"]], r = p[["r"]])
})

#' @rdname refptsSummary
#' @details For a \code{derived_posteriors.sso} file, values are posterior
#'   medians; \code{-mceval} does not report total biomass, so \code{tb0} is
#'   \code{NA}.
#' @export
setMethod("refptsSummary", signature(object = "character"), function(object, ...) {
  d <- mcmcTable(object[1])
  v <- function(k) if (k %in% names(d)) stats::median(as.numeric(d[[k]]), na.rm = TRUE) else NA_real_
  .ptRow(v("Dead_Catch_MSY"), v("annF_MSY"), v("SSB_MSY"), v("SSB_Virgin"), NA_real_)
})

#' @rdname surplusError
#' @details For SS3, \eqn{B} is \code{Bio_all}, \eqn{S} is SSB and \eqn{P} the
#'   equilibrium yield curve from \code{\link{curveSS}}. For JABBA, \eqn{B = S}
#'   is total biomass and \eqn{P} the fitted Pella-Tomlinson curve
#'   (\code{pfunc}).
#' @export
setMethod("surplusError", signature(object = "list"), function(object, ...) {
  type <- .fitType(object)
  if (is.na(type)) return(.byScenario(object, surplusError, ...))
  if (type == "ss3") {
    cs <- curveSS(object)
    ts <- cs$tseries[order(cs$tseries$year), ]
    d  <- data.frame(year = ts$year, biomass = ts$biomass, stock = ts$ssb, catch = ts$yield,
                     pf = ts$pf, bmsy = cs$refpts$bmsy, msy = cs$refpts$msy)
  } else {
    yr <- as.numeric(dimnames(object$timeseries)[[1]])
    b  <- object$timeseries[, "mu", "B"]
    ct <- object$inputseries$catch
    d  <- data.frame(year = yr, biomass = b, stock = b, catch = ct[match(yr, ct[, 1]), 2],
                     pf = stats::splinefun(object$pfunc$SB_i, object$pfunc$SP, method = "natural")(b),
                     bmsy = object$refpts$bmsy[1], msy = object$pfunc$MSY[1])
  }
  d$sp <- c(d$biomass[-1], NA) - d$biomass + d$catch
  d$x  <- d$stock / d$bmsy
  d$pe <- (d$sp - d$pf) / d$msy
  d
})
