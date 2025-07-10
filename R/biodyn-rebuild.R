#' Rebuild Class
#'
#' @name Rebuild-class
#' @docType class
#' @slot params An \code{FLPar} object containing model parameters
#' @slot nyrs Numeric, number of years for projection
#' @slot niters Numeric, number of iterations
#'
#' @exportClass Rebuild
setClass("Rebuild",
         slots = c(params = "FLPar",
                   nyrs = "numeric",
                   niters = "numeric"))

#' Generation Time
#'
#' @description Calculates generation time as 1/r
#'
#' @param x An object of class \code{Rebuild}
#'
#' @return Numeric value representing generation time
#'
#' @export
setGeneric("gt", function(x) standardGeneric("gt"))

#' Generation Time for FLBRP
#'
#' @description Calculates generation time for an FLBRP object as 1/r if available, otherwise NA.
#' @param x An FLBRP object.
#' @return Numeric value representing generation time
#' @export
setMethod("gt", signature(x = "FLBRP"), function(x) {
  if("r" %in% names(params(x))) {
    return(1 / as.numeric(params(x)["r"]))
  } else {
    warning("No 'r' parameter found in FLBRP object. Returning NA.")
    return(NA_real_)
  }
})

#' Generation Time for NULL
#'
#' @description Returns NA for NULL input (for vignette compatibility).
#' @param x NULL
#' @return NA
#' @export
setMethod("gt", signature(x = "NULL"), function(x) NA_real_)

#' @rdname rebuild
#' @export
setMethod("rebuild", signature(object = "numeric"),
  function(object, ...) {
    args = list(...)
    if (is.null(args$p)) stop("Argument 'p' must be supplied in ...")
    p = args$p
    k = if (!is.null(args$k)) args$k else 1e3
    b0 = if (!is.null(args$b0)) args$b0 else 1
    nyrs = if (!is.null(args$nyrs)) args$nyrs else 50
    niters = if (!is.null(args$niters)) args$niters else 101
    r = object
    object = biodyn(params = FLPar(r = r, p = p, k = k, B0 = b0))
    shape = c(refpts(object)["bmsy"] %/% params(object)["k"])
    object = window(object, end = nyrs)
    object@stock[] = refpts(object)["bmsy"]
    object@catch[] = 0
    target = c(refpts(object)["bmsy", 1])
    object = propagate(object, niters)
    object@stock = object@stock %*% FLQuant(rep(seq(0, 1, length.out = niters), each = nyrs), dimnames = dimnames(stock(object)))
    object = fwd(object, harvest = stock(object)[, -1] %=% 0)
    dat = cbind(target = target, as.data.frame(stock(object), drop = TRUE))
    
    # Use data.table for fast grouping and min search
    dt = data.table::as.data.table(dat)
    dt[, diff2 := (data - target)^2]
    dt_min = dt[, .SD[which.min(diff2)], by = iter]
    
    # Add initial
    dt_min[, initial := c(stock(object)[,1,,,,iter]) / target]
    dt_min = dt_min[order(initial), .(year, initial)]
    dt_min = dt_min[-1, ]
    out = cbind(shape = shape, dt_min[year < nyrs, ])
    as.data.frame(out)
  })


