#' Calculate recovery time to BMSY under Pella-Tomlinson dynamics
#'
#' @description Numerically integrates the Pella-Tomlinson production function to estimate the time required 
#' for a population to recover from an initial biomass to BMSY, with optional fishing mortality during recovery.
#'
#' @param object Initial biomass (must be less than BMSY) or FLPar object with parameters
#' @param r Intrinsic rate of population increase (if object is numeric)
#' @param p Shape parameter of the Pella-Tomlinson model (if object is numeric)
#' @param bmsy Biomass at maximum sustainable yield (default=1)
#' @param F Fishing mortality rate during recovery (default=0)
#' @param ... Additional arguments
#' @return Recovery time in years
#' @export
#' @examples
#' rTime(0.3, r=0.2, p=1.5)
#' rTime(0.3, r=0.2, p=1.5, F=0.05)
#' 
#' # With FLPar
#' library(FLCore)
#' params=FLPar(initial=0.3, r=0.2, p=1.5, bmsy=1)
#' rTime(params, F=0.05)
setGeneric("rTime", function(object, ...) standardGeneric("rTime"))

#' @rdname rTime
#' @export
setMethod("rTime", signature(object="numeric"),
  function(object, r, p, bmsy=1, F=0, ...) {
    if (missing(r) || missing(p)) 
      stop("r and p must be provided when object is numeric")
    
    initial=object
    
    # Calculate carrying capacity (virgin biomass) from bmsy and p
    if (abs(p) < 1e-6) 
      K=bmsy * exp(1)  # Fox model
    else 
      K=bmsy * (p+1)^(1/p)
    
    # Effective growth rate under fishing mortality
    r_eff=r - F
    
    # Check if recovery is possible under current fishing mortality
    if (r_eff <= 0) {
      warning("Fishing mortality too high: recovery not possible")
      return(NA_real_)
    }
    
    # Initial and target biomass
    B0=initial
    B1=bmsy
    
    if (abs(p) < 1e-6) {
      # Fox (Gompertz) model with fishing mortality
      t=(log(log(K/B0)) - log(log(K/B1)))/r_eff
      return(as.numeric(t))
    }
    
    # General Pella-Tomlinson with fishing mortality
    B0Kp=(B0/K)^p
    B1Kp=(B1/K)^p
    denom0=1 - B0Kp
    denom1=1 - B1Kp
    
    # Check for valid domain
    if (denom0 <= 0 || denom1 <= 0 || r_eff <= 0) return(NA_real_)
    
    t=(1/(r_eff * p)) * log(denom0/denom1)
    
    return(as.numeric(t))
  })

#' @rdname rTime
#' @export
setMethod("rTime", signature(object="FLPar"),
  function(object, F=0, ...) {
    if (!all(c("initial", "r", "p") %in% names(object)))
      stop("FLPar must contain 'initial', 'r', and 'p' parameters")
    
    bmsy=if ("bmsy" %in% names(object)) object["bmsy"] else FLPar(rep(1, dims(object)$iter))
    
    n_iter=dims(object)$iter
    result=numeric(n_iter)
    
    for (i in 1:n_iter) {
      result[i]=rTime(as.numeric(object["initial", i]), 
                        as.numeric(object["r", i]), 
                        as.numeric(object["p", i]), 
                        as.numeric(bmsy[i]), 
                        F=F)
    }
    
    return(FLPar(result, dimnames=list(iter=1:n_iter)))
  })

#' Find initial biomass for a given recovery time to BMSY
#'
#' @description For a specified number of years, returns the biomass from which it takes exactly that time
#' to recover to BMSY under Pella-Tomlinson dynamics with optional fishing mortality.
#'
#' @param object Desired recovery time in years or FLPar object with parameters
#' @param r Intrinsic rate of population increase (if object is numeric)
#' @param p Shape parameter of the Pella-Tomlinson model (if object is numeric)
#' @param bmsy Biomass at maximum sustainable yield (default=1)
#' @param F Fishing mortality rate during recovery (default=0)
#' @param ... Additional arguments
#' @return Initial biomass
#' @export
#' @examples
#' brebuild(3, r=0.2, p=1.5)
#' brebuild(3, r=0.2, p=1.5, F=0.05)
#' 
#' # With FLPar
#' params=FLPar(nyrs=3, r=0.2, p=1.5, bmsy=1)
#' brebuild(params, F=0.05)
setGeneric("brebuild", function(object, ...) standardGeneric("brebuild"))

#' @rdname brebuild
#' @export
setMethod("brebuild", signature(object="numeric"),
  function(object, r, p, bmsy=1, F=0, ...) {
    if (missing(r) || missing(p)) 
      stop("r and p must be provided when object is numeric")
    
    nyrs=object
    
    # Calculate carrying capacity (virgin biomass) from bmsy and p
    if (abs(p) < 1e-6) 
      K=bmsy * exp(1)  # Fox model
    else 
      K=bmsy * (p+1)^(1/p)
    
    # Effective growth rate under fishing mortality
    r_eff=r - F
    
    # Check if recovery is possible under current fishing mortality
    if (r_eff <= 0) {
      warning("Fishing mortality too high: recovery not possible")
      return(NA_real_)
    }
    
    # Grid search approach - more robust than uniroot
    initial_grid <- seq(0.01, 0.99, by=0.001)
    recovery_times <- sapply(initial_grid, function(x) {
      rt <- rTime(x, r, p, bmsy, F)
      if (is.na(rt) || is.infinite(rt)) return(Inf)
      return(rt)
    })
    
    # Find the closest match
    time_diff <- abs(recovery_times - nyrs)
    min_diff_idx <- which.min(time_diff)
    
    if (is.infinite(recovery_times[min_diff_idx])) {
      warning("No valid solution found")
      return(NA_real_)
    }
    
    # Check if the solution is reasonable (within 1% of target time)
    if (time_diff[min_diff_idx] > nyrs * 0.01) {
      warning("Best solution found differs from target by more than 1%")
    }
    
    return(initial_grid[min_diff_idx])
  })

#' @rdname brebuild
#' @export
setMethod("brebuild", signature(object="FLPar"),
  function(object, F=0, ...) {
    if (!all(c("nyrs", "r", "p") %in% names(object)))
      stop("FLPar must contain 'nyrs', 'r', and 'p' parameters")
    
    bmsy=if ("bmsy" %in% names(object)) object["bmsy"] else FLPar(rep(1, dims(object)$iter))
    
    n_iter=dims(object)$iter
    result=numeric(n_iter)
    
    for (i in 1:n_iter) {
      result[i]=brebuild(as.numeric(object["nyrs", i]), 
                           as.numeric(object["r", i]), 
                           as.numeric(object["p", i]), 
                           as.numeric(bmsy[i]), 
                           F=F)
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
gtR0<-function(R0, r) {
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
#' ages<-0:5
#' lx<-c(1, 0.8, 0.6, 0.4, 0.2, 0.1)
#' mx<-c(0, 0, 2, 3, 1, 0)
#' gtLife(ages, lx, mx)
gtLife<-function(ages, lx, mx) {
  numerator<-sum(ages * lx * mx)
  denominator<-sum(lx * mx)
  return(numerator / denominator)
}








