#' Calculate recovery time to BMSY under Pella-Tomlinson dynamics
#'
#' @description Numerically integrates the Pella-Tomlinson production function to estimate the time required 
#' for a population to recover from an initial biomass to BMSY, assuming zero fishing.
#'
#' @param object Initial biomass (must be less than BMSY) or FLPar object with parameters
#' @param r Intrinsic rate of population increase (if object is numeric)
#' @param p Shape parameter of the Pella-Tomlinson model (if object is numeric)
#' @param bmsy Biomass at maximum sustainable yield (default=1)
#' @param ... Additional arguments
#' @return Recovery time in years
#' @export
#' @examples
#' rTime(0.3, r=0.2, p=1.5)
#' 
#' # With FLPar
#' library(FLCore)
#' params <- FLPar(initial=0.3, r=0.2, p=1.5, bmsy=1)
#' rTime(params)
setGeneric("rTime", function(object, ...) standardGeneric("rTime"))

#' @rdname rTime
#' @export
setMethod("rTime", signature(object="numeric"),
  function(object, r, p, bmsy=1, ...) {
    if (missing(r) || missing(p)) 
      stop("r and p must be provided when object is numeric")
    
    initial <- object
    virgin <- (p+1)^(1/p) * bmsy
    
    integrand <- function(x) 1/(r * x * (1-(x/virgin)^p))
    res <- integrate(integrand, lower=initial, upper=bmsy)
    
    return(res$value)
  })

#' @rdname rTime
#' @export
setMethod("rTime", signature(object="FLPar"),
  function(object, ...) {
    if (!all(c("initial", "r", "p") %in% names(object)))
      stop("FLPar must contain 'initial', 'r', and 'p' parameters")
    
    bmsy <- if ("bmsy" %in% names(object)) object["bmsy"] else FLPar(rep(1, dims(object)$iter))
    
    n_iter <- dims(object)$iter
    result <- numeric(n_iter)
    
    for (i in 1:n_iter) {
      result[i] <- rTime(as.numeric(object["initial", i]), 
                        as.numeric(object["r", i]), 
                        as.numeric(object["p", i]), 
                        as.numeric(bmsy[i]))
    }
    
    return(FLPar(result, dimnames=list(iter=1:n_iter)))
  })

#' Find initial biomass for a given recovery time to BMSY
#'
#' @description For a specified number of years, returns the biomass from which it takes exactly that time
#' to recover to BMSY under Pella-Tomlinson dynamics with zero fishing.
#'
#' @param object Desired recovery time in years or FLPar object with parameters
#' @param r Intrinsic rate of population increase (if object is numeric)
#' @param p Shape parameter of the Pella-Tomlinson model (if object is numeric)
#' @param bmsy Biomass at maximum sustainable yield (default=1)
#' @param ... Additional arguments
#' @return Initial biomass
#' @export
#' @examples
#' brebuild(3, r=0.2, p=1.5)
#' 
#' # With FLPar
#' params <- FLPar(nyrs=3, r=0.2, p=1.5, bmsy=1)
#' brebuild(params)
setGeneric("brebuild", function(object, ...) standardGeneric("brebuild"))

#' @rdname brebuild
#' @export
setMethod("brebuild", signature(object="numeric"),
  function(object, r, p, bmsy=1, ...) {
    if (missing(r) || missing(p)) 
      stop("r and p must be provided when object is numeric")
    
    nyrs <- object
    virgin <- (p+1)^(1/p) * bmsy
    
    objFn <- function(initial) rTime(initial, r, p, bmsy) - nyrs
    res <- uniroot(objFn, lower=1e-6, upper=bmsy-1e-6)
    
    return(res$root)
  })

#' @rdname brebuild
#' @export
setMethod("brebuild", signature(object="FLPar"),
  function(object, ...) {
    if (!all(c("nyrs", "r", "p") %in% names(object)))
      stop("FLPar must contain 'nyrs', 'r', and 'p' parameters")
    
    bmsy <- if ("bmsy" %in% names(object)) object["bmsy"] else FLPar(rep(1, dims(object)$iter))
    
    n_iter <- dims(object)$iter
    result <- numeric(n_iter)
    
    for (i in 1:n_iter) {
      result[i] <- brebuild(as.numeric(object["nyrs", i]), 
                           as.numeric(object["r", i]), 
                           as.numeric(object["p", i]), 
                           as.numeric(bmsy[i]))
    }
    
    return(FLPar(result, dimnames=list(iter=1:n_iter)))
  })

#' Calculate generation time from r and R0
#'
#' @description Calculates generation time using the intrinsic rate of increase and the net reproductive rate.
#'
#' @param R0 Net reproductive rate (mean number of offspring per individual per generation)
#' @param r Intrinsic rate of population increase
#' @return Generation time in years
#' @export
#' @examples
#' gt(R0=2, r=0.2)
gt <- function(R0, r) {
  return(log(R0) / r)
}

#' Calculate generation time from a life table
#'
#' @description Calculates the mean generation time from age-specific survivorship and fecundity.
#'
#' @param ages Vector of ages
#' @param lx Vector of age-specific survivorship (proportion surviving to each age)
#' @param mx Vector of age-specific fecundity (mean number of offspring at each age)
#' @return Generation time in years
#' @export
#' @examples
#' ages <- 0:5
#' lx <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1)
#' mx <- c(0, 0, 2, 3, 1, 0)
#' gtLife(ages, lx, mx)
gtLife <- function(ages, lx, mx) {
  numerator <- sum(ages * lx * mx)
  denominator <- sum(lx * mx)
  return(numerator / denominator)
}

#' Simulate time series of biomass, catch, and harvest rate to a specified biomass target
#'
#' @description Simulates the recovery trajectory of a stock from an initial biomass to a specified target biomass,
#' expressed as a fraction or multiple of BMSY, under the Pella-Tomlinson production model.
#'
#' @param object Initial biomass (relative to BMSY, e.g., 0.3) or FLPar object with parameters
#' @param r Intrinsic rate of population increase (if object is numeric)
#' @param p Shape parameter of the Pella-Tomlinson model (if object is numeric)
#' @param bmsy Biomass at maximum sustainable yield (default=1)
#' @param bmsyTar Target biomass as a fraction or multiple of BMSY (default=1, i.e., BMSY)
#' @param tMax Maximum time (years) to simulate (default=50)
#' @param dt Time step (years) for the simulation (default=0.1)
#' @param hrate Constant fishing mortality rate (default=0)
#' @param ... Additional arguments
#' @return A data.frame with columns: time, biomass, catch, harvest_rate (for numeric) or FLQuant (for FLPar)
#' @export
#' @examples
#' ts <- rTime_series_target(0.3, r=0.2, p=1.5, bmsyTar=0.75)
#' 
#' # With FLPar
#' params <- FLPar(initial=0.3, r=0.2, p=1.5, bmsy=1)
#' q <- rTime_series_target(params, bmsyTar=1.2)
setGeneric("rTime_series_target", function(object, ...) standardGeneric("rTime_series_target"))

#' @rdname rTime_series_target
#' @export
setMethod("rTime_series_target", signature(object="numeric"),
  function(object, r, p, bmsy=1, bmsyTar=1, tMax=50, dt=0.1, hrate=0, ...) {
    if (missing(r) || missing(p)) 
      stop("r and p must be provided when object is numeric")
    
    initial <- object
    virgin <- (p + 1)^(1/p) * bmsy
    target_biomass <- bmsyTar * bmsy
    n_steps <- ceiling(tMax / dt) + 1
    time <- seq(0, tMax, by=dt)
    biomass <- numeric(length(time))
    catch <- numeric(length(time))
    harvest_rate <- numeric(length(time))
    biomass[1] <- initial
    
    for (i in 2:length(time)) {
      catch[i - 1] <- hrate * biomass[i - 1]
      harvest_rate[i - 1] <- ifelse(biomass[i - 1] > 0, catch[i - 1] / biomass[i - 1], 0)
      dB <- r * biomass[i - 1] * (1 - (biomass[i - 1] / virgin)^p) * dt - catch[i - 1] * dt
      biomass[i] <- biomass[i - 1] + dB
      if (biomass[i] < 0) biomass[i] <- 0
      if (biomass[i] >= target_biomass) {
        biomass[i:length(time)] <- target_biomass
        catch[i:length(time)] <- 0
        harvest_rate[i:length(time)] <- 0
        break
      }
    }
    catch[length(time)] <- hrate * biomass[length(time)]
    harvest_rate[length(time)] <- ifelse(biomass[length(time)] > 0, catch[length(time)] / biomass[length(time)], 0)
    
    data.frame(
      time=time,
      biomass=biomass,
      catch=catch,
      harvest_rate=harvest_rate
    )
  })

#' @rdname rTime_series_target
#' @export
setMethod("rTime_series_target", signature(object="FLPar"),
  function(object, bmsyTar=1, tMax=50, dt=0.1, hrate=0, ...) {
    if (!all(c("initial", "r", "p") %in% names(object)))
      stop("FLPar must contain 'initial', 'r', and 'p' parameters")
    
    n_iter <- dims(object)$iter
    r <- object["r"]
    p <- object["p"]
    bmsy <- if ("bmsy" %in% names(object)) object["bmsy"] else FLPar(rep(1, n_iter))
    initial <- object["initial"]
    
    n_steps <- ceiling(tMax / dt) + 1
    nyrs <- seq(0, tMax, by=dt)
    
    out <- FLQuant(dimnames=list(
      quant=c("biomass", "catch", "harvest_rate"),
      year=as.character(nyrs),
      iter=as.character(1:n_iter)
    ))
    
    for(it in 1:n_iter) {
      r_it <- as.numeric(r[it])
      p_it <- as.numeric(p[it])
      bmsy_it <- as.numeric(bmsy[it])
      initial_it <- as.numeric(initial[it])
      virgin_it <- (p_it + 1)^(1/p_it) * bmsy_it
      target_biomass <- bmsyTar * bmsy_it
      
      biomass <- numeric(n_steps)
      catch <- numeric(n_steps)
      harvest_rate <- numeric(n_steps)
      biomass[1] <- initial_it
      
      for (i in 2:n_steps) {
        catch[i-1] <- hrate * biomass[i-1]
        harvest_rate[i-1] <- ifelse(biomass[i-1] > 0, catch[i-1] / biomass[i-1], 0)
        dB <- r_it * biomass[i-1] * (1 - (biomass[i-1] / virgin_it)^p_it) * dt - catch[i-1] * dt
        biomass[i] <- biomass[i-1] + dB
        if (biomass[i] < 0) biomass[i] <- 0
        if (biomass[i] >= target_biomass) {
          biomass[i:n_steps] <- target_biomass
          catch[i:n_steps] <- 0
          harvest_rate[i:n_steps] <- 0
          break
        }
      }
      catch[n_steps] <- hrate * biomass[n_steps]
      harvest_rate[n_steps] <- ifelse(biomass[n_steps] > 0, catch[n_steps] / biomass[n_steps], 0)
      
      out["biomass",,it] <- biomass
      out["catch",,it] <- catch
      out["harvest_rate",,it] <- harvest_rate
    }
    return(out)
  })

# Legacy function name for backward compatibility
#' @rdname rTime_series_target
#' @export
rTime_series_FLR <- function(params, bmsyTar=1, tMax=50, dt=0.1, hrate=0) {
  .Deprecated("rTime_series_target", package="FLRebuild")
  rTime_series_target(params, bmsyTar=bmsyTar, tMax=tMax, dt=dt, hrate=hrate)
}





