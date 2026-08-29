#' Calculate Exploitable Biomass
#'
#' @description
#' Calculates the exploitable biomass from an FLStock object using selectivity-weighted
#' catch weights and stock numbers.
#'
#' @param object An object of class FLStock
#'
#' @return An FLQuant object containing the exploitable biomass time series
#'
#' @details
#' Exploitable biomass is based on weighting catch weights by selectivity normalised 
#' by peak selectivity, then multiplying by stock numbers and summing across ages
#'
#' @examples
#' \dontrun{
#' data(ple4)
#' eb <- ebiomass(ple4)
#' }
#'
#' @export
setGeneric("ebiomass", function(object) standardGeneric("ebiomass"))

#' @rdname ebiomass
#' @export
setMethod("ebiomass", signature(object="FLStock"),
          function(object) {
            sel   <- FLCore::harvest(object)
            wt    <- catch.wt(object) %*% sel %/% FLCore::fapex(sel)
            eb.wt <- qmax(wt, 0.000001)
            
            apply(eb.wt %*% stock.n(object), 2:6, sum)
          })

setMethod("ebiomass", signature(object="FLBRP"),
          function(object) {
            sel   <- FLCore::harvest(object)
            wt    <- catch.wt(object) %*% sel %/% FLCore::fapex(sel)
            eb.wt <- qmax(wt, 0.000001)
            
            apply(eb.wt %*% stock.n(object), 2:6, sum)
          })
