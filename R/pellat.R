#' Calculate Pella-Tomlinson Production
#'
#' @description Calculate production (surplus production) for given Pella-Tomlinson parameters and biomass.
#' This is the main function for calculating production in the Pella-Tomlinson surplus production model.
#'
#' @param params FLPar object with r, p, virgin parameters
#' @param biomass FLQuant object with biomass values
#' @param ... Additional arguments
#' @return FLQuant object with production values
#' @export
#' @examples
#' params = FLPar(r = 0.5, p = 1, virgin = 1000)
#' biomass = FLQuant(500)
#' result = production(params, biomass)
setGeneric("production",
           function(params, biomass, ...) standardGeneric("production"))

#' @rdname production
#' @export
setMethod("production", signature(params="FLPar", biomass="FLQuant"),
          function(params, biomass) {
            # Calculate production using Pella-Tomlinson formula
            production = params["r"] %*% biomass %*% (1 - exp(log((biomass %/% params["virgin"])) %*% params["p"])) %/% params["p"]
            
            # Handle Fox model (p == 0) as special case
            if (any((params["p"] == 0))) {
              prm = iter(params["r"], params["p"] == 0)
              production[,,,,, params["p"] == 0] = prm["r"] %*% biomass[,,,,, params["p"] == 0] %*% log(prm["virgin"] %*/% biomass)
            }
            
            # Ensure non-negative production
            return(qmax(as.FLQuant(production), 0))
          })


#' Pella-Tomlinson Parameter Calculator
#'
#' @description Main function to calculate Pella-Tomlinson parameters from various input combinations.
#' This is the primary function for parameter calculation.
#'
#' @param object FLPar object containing parameters
#' @param biomass Character string specifying biomass type (optional)
#' @param ... Additional arguments
#' @return FLPar object with calculated parameters
#' @export
#' @examples
#' params = FLPar(bmsy = 500, k = 1000, fmsy = 0.2)
#' result = pellat(params)
setGeneric("pellat", function(object, biomass, ...) standardGeneric("pellat"))

#' @rdname pellat
#' @export
setMethod("pellat", signature(object="FLPar"),
          function(object){ 
            
            # Wrap the entire function in error handling for robustness
            tryCatch({
              if ("virgin"%in%dimnames(object)$params)
                dimnames(object)$params["virgin"==dimnames(object)$params]="k"
              
              if ("fmsyMedianC"%in%dimnames(object)$params)
                dimnames(object)$params["fmsyMedianC"==dimnames(object)$params]="fmsy"
              
              # Check if we have bmsy and k (virgin) for bmsyK2p method
              hasBmsyK = all(c("bmsy","k") %in% dimnames(object)$params)
              # Check if we have r and fmsy for rFmsy2p method
              hasRFmsy = all(c("r","fmsy") %in% dimnames(object)$params)
              
              # Always prefer bmsyK2p method when available, as it's more robust
              if (hasBmsyK) {
                # Use bmsyK2p method (preferred)
                hasRFmsy = FALSE  # Don't use rFmsy2p even if available
              } else if (!hasRFmsy) {
                warning("FLPar must contain either (bmsy, k) or (r, fmsy) parameters")
                return(FLPar(array(NA, dim = c(7, 1),
                                  dimnames = list(params = c("r", "p", "k", "virgin", "bmsy", "fmsy", "msy"),
                                                iter = 1))))
              }
            
            # Get dimensions
            dims = dim(object)[2]
            
                         # Initialize output FLPar
             res = FLPar(array(NA, dim = c(7, dims),
                               dimnames = list(params = c("r", "p", "k", "virgin", "bmsy", "fmsy", "msy"),
                                             iter = seq(dims))))
            
            # Loop over iterations
            for(i in seq(dims)) {
              
                             if (hasBmsyK) {
                # Use bmsyK2p method (preferred)
                bmsyK = c(object["bmsy", i]) / c(object["k", i])
                
                # Check if bmsyK ratio is valid for Pella-Tomlinson model
                if (bmsyK <= 0 || bmsyK >= 1) {
                  warning("Invalid bmsyK ratio: ", bmsyK, ". Must be between 0 and 1.")
                  p = NA
                } else {
                  p = bmsyK2p(bmsyK)
                  # Check if bmsyK2p returned NA (calculation failed)
                  if (is.na(p)) {
                    warning("Failed to calculate p for bmsyK ratio: ", bmsyK)
                  }
                }
                
                # Get original MSY and BMSY for consistency
                original_msy = if ("msy" %in% dimnames(object)$params) {
                  c(object["msy", i])
                } else {
                  # Calculate from fmsy and bmsy
                  if ("fmsy" %in% dimnames(object)$params) {
                    c(object["fmsy", i]) * c(object["bmsy", i])
                  } else {
                    NA
                  }
                }
                
                original_bmsy = c(object["bmsy", i])
                k = c(object["k", i])
                
                # Check if p is at a bound (indicating bmsyK2p failed)
                p_at_bound = !is.na(p) && (abs(p - 0.05) < 1e-6)  # Check if p is at lower bound
                
                if (p_at_bound && !is.na(original_msy)) {
                  # Re-estimate parameters to maintain MSY consistency
                  # When p is at bound, we need to adjust r to maintain the original MSY
                  
                  # Calculate the new BMSY that corresponds to p = 0.05
                  if (abs(p) < 1e-10) {
                    # Fox model: bmsy = k * exp(-1)
                    bmsy = k * exp(-1)
                  } else {
                    # Pella-Tomlinson model: bmsy = k * (1/(p+1))^(1/p)
                    bmsy = k * (1/(p+1))^(1/p)
                  }
                  
                  # Calculate r directly from MSY relationship: MSY = r * BMSY * (1 - (BMSY/K)^p) / p
                  # Rearranging: r = MSY * p / (BMSY * (1 - (BMSY/K)^p))
                  if (abs(p) < 1e-10) {
                    # Fox model: MSY = r * BMSY * log(K/BMSY)
                    # Rearranging: r = MSY / (BMSY * log(K/BMSY))
                    r = original_msy / (bmsy * log(k / bmsy))
                  } else {
                    # Pella-Tomlinson model
                    r = original_msy * p / (bmsy * (1 - (bmsy/k)^p))
                  }
                  
                  # Calculate fmsy from the production relationship
                  fmsy = original_msy / bmsy
                  
                  # Use the original MSY
                  msy = original_msy
                  
                  warning("p set to lower bound (0.05). Re-estimated r to maintain original MSY at BMSY.")
                  
                } else {
                  # Standard calculation when p is not at bound
                  if ("fmsy" %in% dimnames(object)$params) {
                    fmsy = c(object["fmsy", i])
                  } else if ("msy" %in% dimnames(object)$params) {
                    fmsy = c(object["msy", i]) / c(object["bmsy", i])
                  } else {
                    stop("Need either fmsy or msy parameter to calculate r")
                  }
                  
                  # Only calculate r if p is valid
                  if (!is.na(p)) {
                    r = fmsy * (p + 1) / p
                  } else {
                    r = NA
                  }
                  
                  # Calculate bmsy from p and k
                  if (is.na(p)) {
                    bmsy = NA
                  } else if (abs(p) < 1e-10) {
                    # Fox model: bmsy = k * exp(-1)
                    bmsy = k * exp(-1)
                  } else {
                    # Pella-Tomlinson model: bmsy = k * (1/(p+1))^(1/p)
                    bmsy = k * (1/(p+1))^(1/p)
                  }
                  
                  # Calculate msy = fmsy * bmsy
                  if (is.na(bmsy) || is.na(fmsy)) {
                    msy = NA
                  } else {
                    msy = fmsy * bmsy
                  }
                }
                
              } else {
                # Use rFmsy2p method (fallback)
                r = c(object["r", i])
                fmsy = c(object["fmsy", i])
                p = rFmsy2p(r, fmsy)
                # Check if rFmsy2p returned NA (calculation failed)
                if (is.na(p)) {
                  warning("Failed to calculate p for r=", r, ", fmsy=", fmsy)
                }
                
                # Calculate k from r and fmsy
                if (is.na(p)) {
                  k = NA
                  bmsy = NA
                  msy = NA
                } else {
                  # Calculate k using the relationship k = r / fmsy * (p + 1) / p
                  k = r / fmsy * (p + 1) / p
                  
                  # Calculate bmsy from p and k
                  if (abs(p) < 1e-10) {
                    # Fox model: bmsy = k * exp(-1)
                    bmsy = k * exp(-1)
                  } else {
                    # Pella-Tomlinson model: bmsy = k * (1/(p+1))^(1/p)
                    bmsy = k * (1/(p+1))^(1/p)
                  }
                  
                  # Calculate msy = fmsy * bmsy
                  msy = fmsy * bmsy
                }
              }
              
              # Store results
              res["r", i] = r
              res["p", i] = p
              res["k", i] = k
              res["virgin", i] = k
              res["bmsy", i] = bmsy
              res["fmsy", i] = fmsy
              res["msy", i] = msy
            }
            
            return(res)
            
          }, error = function(e) {
            warning("Error in pellat calculation: ", e$message)
            return(FLPar(array(NA, dim = c(7, 1),
                              dimnames = list(params = c("r", "p", "k", "virgin", "bmsy", "fmsy", "msy"),
                                            iter = 1))))
          })
          })

#' @rdname pellat
#' @export
setMethod("pellat", signature(object="FLBRP"),
          function(object, biomass = "ssb", ...) {
            # Extract reference points from FLBRP object
            refs = refptsEB(object)[c("msy", "crash", "virgin"),
                                    c("harvest", "yield", "ssb", "eb"), drop = TRUE]
            
            # Create FLPar with extracted parameters
            params = FLPar(
              msy = refs[1, "yield"],
              bmsy = refs[1, biomass],
              k = refs[3, biomass],
              fmsy = refs[1, "harvest"]
            )
            
            # Call pellat with FLPar
            pellat(params)
          })

#' @rdname pellat
#' @export
setMethod("pellat", signature(object="FLBRPs"),
          function(object, biomass = "ssb", ...) {
            # Apply pellat to each FLBRP in the collection
            results = lapply(object, function(x) pellat(x, biomass = biomass, ...))
            
            # Return results as a list
            results
          })

#' Direct calculation of Pella-Tomlinson parameters
#' 
#' @param msy Maximum sustainable yield
#' @param bmsy Biomass at MSY
#' @param virgin Virgin biomass
#' @param ... Additional arguments
#' @return List with estimated parameters
#' @export
#' @examples
#' result = fitPellatDirect(msy=100, bmsy=500, virgin=1000)
fitPellatDirect = function(msy, bmsy, virgin, ...) {
  
  if (missing(msy) || missing(bmsy) || missing(virgin)) {
    stop("msy, bmsy, and virgin must be provided")
  }
  
  if (msy <= 0 || bmsy <= 0 || virgin <= 0) {
    stop("msy, bmsy, and virgin must be positive")
  }
  
  if (bmsy >= virgin) {
    stop("bmsy must be less than virgin")
  }
  
  # Calculate bmsy/virgin ratio
  shape = bmsy / virgin
  
  # Use the standardized bmsyK2p function for consistency
  p = bmsyK2p(shape)
  
  # Check if p calculation was successful
  if (is.na(p)) {
    stop("Failed to calculate p for bmsy/virgin ratio: ", shape)
  }
  
  # Calculate fmsy and r
  fmsy = msy/bmsy
  r = fmsy * (p + 1) / p
  
  # Validate parameters
  if (r <= 0) {
    stop("Calculated r <= 0 (r = ", r, "), this indicates invalid input parameters")
  }
  
  # Return in the same format as pellatParams for consistency
  return(list(
    r = r,
    k = virgin,  # Use 'k' instead of 'virgin' to match pellatParams
    p = p,
    virgin = virgin,
    msy = msy,
    bmsy = bmsy,
    convergence = 0,
    logLik = NA,
    method = "direct"
  ))
}

#' Pella-Tomlinson Maximum Likelihood Estimation
#'
#' Fits Pella-Tomlinson production model using maximum likelihood estimation.
#' Dispatched on different object types: data.frame, FLBRP, FLBRPs, and list.
#'
#' @param object The object to fit the model to:
#'   - data.frame: Must contain biomass and yield columns
#'   - FLBRP: Single FLBRP object
#'   - FLBRPs: Collection of FLBRP objects
#'   - list: List of FLBRP objects
#' @param biomass Character string specifying the biomass column name (default: "ssb")
#' @param ... Additional arguments passed to optimization
#' @return For single objects: FLPar object with estimated parameters
#'         For collections: List of FLPar objects
#' @export
#' @examples
#' # data.frame method
#' df = data.frame(ssb = c(1000, 2000, 3000, 4000, 5000),
#'                 yield = c(100, 200, 250, 200, 150))
#' result = pellatMle(df, biomass = "ssb")
setGeneric("pellatMle", function(object, biomass = "ssb", ...) standardGeneric("pellatMle"))

#' @rdname pellatMle
#' @export
setMethod("pellatMle", "data.frame",
          function(object, biomass = "ssb", ...) {
            # Validate required columns
            if (!biomass %in% names(object)) {
              stop("Biomass column '", biomass, "' not found in data.frame")
            }
            if (!"yield" %in% names(object)) {
              stop("Yield column 'yield' not found in data.frame")
            }
            
            pellaTMleCore(object, biomass = biomass, ...)
          })

#' @rdname pellatMle
#' @export
setMethod("pellatMle", "FLBRP",
          function(object, biomass = "ssb", ...) {
            # Extract data from FLBRP object
            df = model.frame(FLQuants(object, ssb = ssb, eb = FLCandy:::ebiomass, yield = catch), drop = TRUE)
            
            # Filter out NA values and ensure positive values
            df = df[!is.na(df[[biomass]]) & !is.na(df[["yield"]]) & 
                    df[[biomass]] > 0 & df[["yield"]] > 0, ]
            
            pellaTMleCore(df, biomass = biomass, ...)
          })

#' @rdname pellatMle
#' @export
setMethod("pellatMle", "FLBRPs",
          function(object, biomass = "ssb", ...) {
            # Apply pellatMle to each FLBRP in the collection
            results = lapply(object, function(x) pellatMle(x, biomass = biomass, ...))
            
            # Return results as a list
            results
          })

#' @rdname pellatMle
#' @export
setMethod("pellatMle", "list",
          function(object, biomass = "ssb", ...) {
            # Check if all elements are FLBRP objects
            if (!all(sapply(object, function(x) inherits(x, "FLBRP")))) {
              stop("All elements in list must be FLBRP objects")
            }
            
            # Apply pellatMle to each FLBRP in the list
            results = lapply(object, function(x) pellatMle(x, biomass = biomass, ...))
            
            # Return results as a list
            results
          })

#' Core MLE function for Pella-Tomlinson parameter estimation
#' @keywords internal
pellaTMleCore = function(object, biomass = "ssb", ...) {
  
  # Filter out NA values and ensure positive values
  clean = object[!is.na(object[[biomass]]) & !is.na(object[["yield"]]) & 
               object[[biomass]] > 0 & object[["yield"]] > 0, ]
  
  if (nrow(clean) < 3) {
    warning("Insufficient data for Pella-Tomlinson fitting (need at least 3 valid points)")
    return(NULL)
  }
  
  # Extract biomass and yield
  B = clean[[biomass]]
  Y = clean[["yield"]]
  
  # MSY and BMSY
  msyIdx = which.max(Y)
  MSY = Y[msyIdx]
  BMSY = B[msyIdx]
  
  K = max(B)
  
  # Determine if p should be positive or negative based on data
  # If BMSY/K > 0.37 (e^-1), then p should be negative
  # If BMSY/K < 0.37, then p should be positive
  bmsyKRatio = BMSY / K
  pShouldBeNegative = bmsyKRatio > 0.37
  
  # Set appropriate bounds for p based on the data
  if (pShouldBeNegative) {
    pLower = -5
    pUpper = -0.05
    pInitial = -1  # Start with p = -1 (Fox model)
  } else {
    pLower = 0.05
    pUpper = 5
    pInitial = 1   # Start with p = 1 (Schaefer model)
  }
  
  # Estimate r using MSY relationship: MSY=r*K/(p+1)^(1/p)
  # Use the initial p value for better r estimate
  if (abs(pInitial) < 1e-4) {
    rInitial = MSY / (K * exp(-1))  # Fox model
  } else {
    rInitial = MSY * (abs(pInitial) + 1)^(1/pInitial) / K
  }
  
  # Function to calculate production
  calcProd = function(params) {
    r = params[1]
    p = params[2]
    
    # Calculate production
    P = r * B * (1 - (B/K)^p) / p
    
    # Negative log-likelihood (assuming normal errors)
    -sum(dnorm(Y, P, sd(Y) * 0.1, log=TRUE))
  }
  
  # Optimize parameters with correct bounds
  tryCatch({
    opt = optim(c(rInitial, pInitial), calcProd, 
                method = "L-BFGS-B", 
                lower = c(0.01, pLower), 
                upper = c(3, pUpper))
    
    if (opt$convergence == 0) {
      return(FLPar(list(
        r = opt$par[1],
        p = opt$par[2],
        k = K,
        virgin = K,
        bmsy = BMSY,
        fmsy = MSY / BMSY,
        msy = MSY
      )))
    } else {
      warning("MLE optimization failed to converge")
      return(NULL)
    }
    
  }, error = function(e) {
    warning("Error during MLE fitting: ", e$message)
    return(NULL)
  })
} 