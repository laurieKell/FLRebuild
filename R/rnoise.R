#' Autocorrelated recruitment residuals
#'
#' Methods for \code{\link[FLCore]{rnoise}} when \code{n} is an \code{FLQuant}
#' of residuals rather than a number of iterations. An AR(1) is fitted to the
#' residuals after they have been standardised by the regime mean and standard
#' deviation from \code{\link{recDevs}}.
#'
#' @param n \code{FLQuant} of residuals, typically log stock--recruit
#'   residuals. The two-argument method requires a single iteration; that
#'   series supplies the autocorrelation.
#' @param len \code{FLQuant} template for the years and iterations to fill.
#'   If omitted, the result has the same dimensions as \code{n}.
#' @param ... Passed along by the \pkg{FLCore} generic. These methods do not
#'   use further arguments.
#'
#' @return An \code{FLQuant}.
#' \describe{
#'   \item{\code{len} missing}{Each iteration is replaced by an AR(1) draw
#'     and the regime mean is added back.}
#'   \item{\code{len} supplied}{An AR(1) series with the dimensions of
#'     \code{len}, on the standardised-residual scale. Callers that want
#'     multiplicative deviations exponentiate the result.}
#' }
#'
#' @details
#' Standardisation uses \code{\link{recDevs}}. The AR(1) coefficient and
#' innovation standard deviation are estimated with \code{stats::arima} on
#' each iteration, and the draws themselves use the numeric method
#' \code{\link[FLCore]{rnoise}}.
#'
#' @seealso \code{\link{recDevs}}, \code{\link[FLCore]{rnoise}}
#'
#' @examples
#' \dontrun{
#' data(nsher, package = "FLCore")
#' rsdl <- recDevs(nsher)$residuals
#' rnoise(rsdl)
#' tmpl <- FLQuant(0, dimnames = list(year = 2000:2020, iter = 1:5))
#' rnoise(iter(rsdl, 1), tmpl)
#' }
#'
#' @name rnoise
#' @aliases rnoise,FLQuant,missing-method rnoise,FLQuant,FLQuant-method
#' @export
setMethod("rnoise", signature(n = "FLQuant", len = "missing"),
          function(n, ...) {
            rs1 <- .recDevsFromQuant(n)
            rs2 <- n %-% rs1$mn
            rs2 <- rs2 %/% rs1$sd

            ft <- plyr::alply(rs2, 6, arima, order = c(1, 0, 0))
            sd <- plyr::laply(ft, function(x) unname(sqrt(x$sigma2)))
            b  <- plyr::laply(ft, function(x) unname(coef(x)["ar1"]))

            rtn <- n
            for (i in seq_along(sd))
              FLCore::iter(rtn, i) <- rnoise(n = 1,
                                             FLCore::iter(n, 1) %=% 0,
                                             sd = sd[i], b = b[i])
            rtn %+% rs1$mn
          })

#' @rdname rnoise
#' @export
setMethod("rnoise", signature(n = "FLQuant", len = "FLQuant"),
          function(n, len, ...) {
            if (dim(n)[6] != 1)
              stop("n must have one iteration; len sets the iterations of the result")

            rs1 <- .recDevsFromQuant(n)
            rs2 <- n %-% rs1$mn[, 1]
            rs2 <- rs2 %/% rs1$sd[, 1]

            ft <- plyr::alply(rs2, 6, arima, order = c(1, 0, 0))
            sd <- plyr::laply(ft, function(x) unname(sqrt(x$sigma2)))
            b  <- plyr::laply(ft, function(x) unname(coef(x)["ar1"]))

            rtn <- len
            for (i in seq(dim(len)[6]))
              FLCore::iter(rtn, i) <- rnoise(n = 1,
                                             FLCore::iter(len, 1) %=% 0,
                                             sd = sd, b = b)
            rtn
          })
