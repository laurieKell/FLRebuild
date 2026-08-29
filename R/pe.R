# =============================================================================
# Process error: location of η in Pella–Tomlinson / surplus-production dynamics
# =============================================================================
# Three placements (catch C_y replaces F_y B_y), with
#   SP_y = (r/p) B_y (1 - (B_y/K)^p)   for Pella–Tomlinson
# or equilibrium yield from FLBRP for production = "FLBRP":
#
#   after catch:        B_{y+1} = (B_y + SP_y - C_y) e^{η_y}
#   before catch:       B_{y+1} = (B_y + SP_y) e^{η_y} - C_y
#   on productivity:    B_{y+1} = B_y + SP_y e^{η_y} - C_y
#
# Solving for the multiplicative residual e^{η_y}:
#   after:         e^η = B_{y+1} / (B_y + SP_y - C_y)
#   before:        e^η = (B_{y+1} + C_y) / (B_y + SP_y)
#   productivity:  e^η = (B_{y+1} - B_y + C_y) / SP_y

#' Pella–Tomlinson surplus production at biomass
#' SP = (r/p) B (1 - (B/K)^p)
#' @noRd
.ptSpAtB <- function(B, params) {
  pn <- dimnames(params)$params
  r <- params["r"]
  p <- params["p"]
  k <- if ("k" %in% pn) params["k"] else if ("virgin" %in% pn)
    params["virgin"]
  else
    stop(".ptSpAtB: params need r, p and k (or virgin)", call. = FALSE)
  ## same form as mpb pellat / production(FLQuant, FLPar)
  ((r %/% p) %*% B) %*% (1 - exp(log(B %/% k) %*% p))
}

#' Align B_{y+1} with year y (last year NA)
#' @noRd
.peBnext <- function(B) {
  Bnext <- B %=% NA_real_
  ny <- dim(B)[2L]
  if (ny > 1L)
    Bnext[, seq_len(ny - 1L)] <- B[, seq.int(2L, ny)]
  Bnext
}

#' Multiplicative / log / relative PE from B, C, SP and error location
#' @noRd
.peFromDynamics <- function(B, C, SP, location = "after", scale = "log") {
  location <- match.arg(location, c("after", "before", "productivity"))
  scale    <- match.arg(scale, c("log", "mult", "relative", "legacy", "raw"))
  Bnext <- .peBnext(B)

  mult <- switch(location,
    after = Bnext %/% (B %+% SP %-% C),
    before = (Bnext %+% C) %/% (B %+% SP),
    productivity = (Bnext %-% B %+% C) %/% SP)

  if (scale %in% c("relative", "raw")) {
    if (identical(location, "after"))
      return((Bnext %-% B %+% C %-% SP) %/% B)
    return(mult %-% 1)
  }
  switch(scale,
    mult = mult,
    log  = {
      m <- c(mult)
      eta <- rep(NA_real_, length(m))
      ok <- is.finite(m) & m > 0
      eta[ok] <- log(m[ok])
      FLQuant(eta, dimnames = dimnames(mult))
    },
    legacy = {
      if (!identical(location, "after"))
        stop("scale = 'legacy' is only defined for location = 'after'",
             call. = FALSE)
      (B %-% Bnext %-% C %+% SP) %/% B
    })
}

#' PT params from eq (FLPar, PellaTomlinson, or FLBRP via pellatParams)
#' @noRd
.ptParamsFromEq <- function(eq) {
  if (is(eq, "FLPar"))
    return(eq)
  if (is(eq, "PellaTomlinson"))
    return(eq@params)
  if (!is(eq, "FLBRP"))
    stop("For production = 'pellat', 'eq' must be FLBRP, FLPar or PellaTomlinson",
         call. = FALSE)
  fn <- get0("pellatParams", envir = asNamespace("FLRebuild"), inherits = FALSE)
  if (is.null(fn))
    stop("pellatParams not found — load or reinstall FLRebuild (devtools::load_all())",
         call. = FALSE)
  fn(eq)
}

#' Equilibrium yield (SP) at observed biomass from an FLBRP
#' @noRd
.flbrpSp <- function(stk, eq, stock = FLCore::ssb) {
  fbar(eq) <- FLQuant(seq(0, 1, length.out = 201)) *
    computeRefpts(eq)["crash", "harvest"]
  mf <- model.frame(FLQuants(eq,
                             stock = function(x) stock(x),
                             catch = function(x) catch(x)),
                    drop = TRUE)
  dat <- with(mf, approx(stock, catch, xout = c(stock(stk))))
  FLQuant(dat$y, dimnames = dimnames(stock(stk)))
}

#' Surplus production series for pe()
#' @noRd
.peSurplus <- function(object, eq, stockFn, production) {
  production <- match.arg(production, c("FLBRP", "pellat"))
  B <- stockFn(object)
  if (identical(production, "FLBRP")) {
    if (!is(eq, "FLBRP"))
      stop("production = 'FLBRP' requires an FLBRP in 'eq'", call. = FALSE)
    return(.flbrpSp(object, eq, stock = stockFn))
  }
  .ptSpAtB(B, .ptParamsFromEq(eq))
}

#' @rdname pe
#' @param stock Function returning the biomass metric (default \code{FLCore::ssb}).
#' @param location Where process error enters the update (catch replaces
#'   \eqn{F_y B_y}):
#'   \itemize{
#'     \item \code{"after"}: \eqn{B_{y+1}=(B_y+\mathrm{SP}_y-C_y)e^{\eta_y}}
#'     \item \code{"before"}: \eqn{B_{y+1}=(B_y+\mathrm{SP}_y)e^{\eta_y}-C_y}
#'     \item \code{"productivity"}: \eqn{B_{y+1}=B_y+\mathrm{SP}_y e^{\eta_y}-C_y}
#'   }
#'   For Pella–Tomlinson,
#'   \eqn{\mathrm{SP}_y=(r/p)B_y\bigl(1-(B_y/K)^p\bigr)}.
#' @param scale Residual scale: \code{"log"} (default; \eqn{\eta_y}),
#'   \code{"mult"} (\eqn{e^{\eta_y}}),
#'   \code{"legacy"} (historical
#'   \code{(B_y-B_{y+1}-C_y+\mathrm{SP}_y)/B_y} for \code{location="after"}),
#'   \code{"relative"} / \code{"raw"} (\eqn{(B_{y+1}-B_y+C_y-\mathrm{SP}_y)/B_y}
#'   for \code{"after"}; \eqn{e^{\eta_y}-1} otherwise).
#' @param production Surplus production source: \code{"FLBRP"} (default;
#'   equilibrium yield at observed biomass from \code{eq}) or \code{"pellat"}
#'   (Pella–Tomlinson via \code{\link{pellatParams}}). For an \code{FLBRP}
#'   \code{eq}, the default \code{"FLBRP"} matches \code{is(eq)[1]}.
#'
#' @return An \code{FLQuant} of process-error residuals by year.
#'
#' @examples
#' \dontrun{
#' data(ple4)
#' eq <- brp(FLBRP(ple4))
#' ## default: after catch, FLBRP production, log residuals
#' pe(ple4, eq)
#' ## Pella–Tomlinson, three locations
#' pe(ple4, eq, production = "pellat", location = "after")
#' pe(ple4, eq, production = "pellat", location = "before")
#' pe(ple4, eq, production = "pellat", location = "productivity")
#' }
#' @export
setMethod("pe", signature(object = "FLStock", eq = "FLBRP"),
          function(object, eq, stock = FLCore::ssb,
                   location = "after",
                   scale = "log",
                   production = "FLBRP") {
            stockFn <- stock
            ## Preserve exact historical formula for the original defaults
            if (identical(location, "after") && identical(scale, "legacy") &&
                identical(production, "FLBRP")) {
              return((stockFn(object) %-%
                        window(stockFn(object)[, -1],
                               end = dims(object)$maxyear + 1) -
                        catch(object) %+% .flbrpSp(object, eq, stock = stockFn)) %/%
                       stockFn(object))
            }
            B  <- stockFn(object)
            C  <- catch(object)
            SP <- .peSurplus(object, eq, stockFn, production)
            .peFromDynamics(B, C, SP, location = location, scale = scale)
          })

#' @rdname pe
#' @export
setMethod("pe", signature(object = "FLStock", eq = "FLPar"),
          function(object, eq, stock = FLCore::ssb,
                   location = "after",
                   scale = "log",
                   production = "pellat") {
            if (!identical(production, "pellat"))
              stop("FLPar method requires production = 'pellat'", call. = FALSE)
            stockFn <- stock
            B  <- stockFn(object)
            C  <- catch(object)
            SP <- .ptSpAtB(B, eq)
            .peFromDynamics(B, C, SP, location = location, scale = scale)
          })
