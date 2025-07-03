#' Rebuild a fish population
#' 
#' @description Projects rebuilding trajectories from different initial SSB levels
#'
#' @param object An object representing the population
#' @param targetF Target fishing mortality during rebuilding
#' @param targetSSB Target spawning stock biomass 
#' @param nInitial Number of initial SSB levels
#' @param growthRate Growth rate for depletion sequence
#' @param minVal Minimum depletion value
#' @param maxVal Maximum depletion value 
#' @param burnin Number of years for burn-in period
#' @param truncate Whether to remove burn-in period
#' @return An object with rebuilding trajectories
#' @export
setGeneric("rebuild", function(object, ...) standardGeneric("rebuild"))

#' Calculate rebuilding time
#'
#' @param object An object containing rebuilding trajectories
#' @param ... Additional arguments
#' @return A data frame with columns year and initial
#' @export
setGeneric("rebuildTime", function(object, ...) standardGeneric("rebuildTime"))

#' Calculate Blim Reference Point
#' 
#' @param object An FLBRP object
#' @param ratio Ratio of virgin recruitment to use (default 0.3)
#' @return An FLPar object containing Blim reference points
setGeneric("blim", function(object, ...) standardGeneric("blim"))


#' Calculate MSY and Virgin State Metrics
#'
#' @description
#' Generic for calculating key metrics (F, SSB, catch, ebiomass) at virgin and MSY states for an FLBRP object.
#'
#' @param object An object of class FLBRP.
#' @return A named vector containing metrics for virgin and MSY states.
#' @note No method is currently implemented for this generic.
#' @export
setGeneric("msyVirgin", function(object, ...) standardGeneric("msyVirgin"))

#' Calculate Process Error for FLBRP Object

#' Calculate the reference age for a FLBRP object.
#'
#' This function calculates the reference age for a FLBRP object.
#'
#' @param object FLBRP object.
#' @param ref Reference point, e.g., "msy" (default) or "f0.1".
#' @param p Probability threshold (default = 0.9).
#'
#' @return An FLQuant object containing reference ages.
#'
#' @export
setGeneric("abiAge", function(object, ref = "msy", p = 0.9) {
  standardGeneric("abiAge")
})


#' Calculate P obs for a FLStock object.
#'
#' This function calculates P obs for a FLStock object.
#'
#' @param object A FLStock object.
#' @param age Reference ages obtained from abiAge.
#'
#' @return An FLQuant object containing P obs.
#'
#' @export
setGeneric("abi", function(object, age, ...) {
  standardGeneric("abi")})


#' Calculate P(N) at FMSY for a FLBRP object.
#'
#' This function calculates P(N) at FMSY for a FLBRP object.
#'
#' @param y A FLBRP object.
#' @param ref Reference point, e.g., "msy" (default) or "f0.1".
#' @param p Probability threshold (default = 0.9).
#'
#' @return An FLQuant object containing P(N) at FMSY.
#'
#' @export
setGeneric("abiMsy", function(object, ref = "msy", p = 0.9) {
  standardGeneric("abiMsy")
})

