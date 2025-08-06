refs<-function(x,biomass="ssb"){
  res=refptsEB(x)[c("msy","crash","virgin"),
                  c("harvest","yield","ssb","eb"),drop=TRUE]
  
  c(msy   =res[1,"yield"],
    bmsy  =res[1,biomass],
    virgin=res[3,biomass],
    fcrash=res[2,"harvest"],
    fmsy  =res[1,"harvest"])}

bmsyK<-function(p) {
  if (p <= -1) stop("Shape parameter p must be > -1.")
  
  (1/(1+p))^(1/p)}



#' pellatParams S4 Class
#'
#' S4 class for estimating and storing Pella-Tomlinson parameters and reference points.
#' Uses FLPar objects for proper unit handling and FLR framework integration.
#'
#' @slot obs FLPar object containing observed reference values (length 3).
#'   Parameter names must be one of: r, virgin, bmsy, msy, fmsy, fcrash.
#' @slot comparison FLPar object comparing observed vs estimated parameters.
setClass("pellatParams",
         slots = c(
           obs = "FLPar",
           comparison = "FLPar"),
         validity = function(object) {
           errors = character()
           
           # Validate obs slot
           if (!inherits(object@obs, "FLPar")) {
             errors = c(errors, "obs slot must be an FLPar object")
           } else {
                           # Check obs dimensions
              if (length(dim(object@obs)) != 2 || dim(object@obs)[1] != 3 || dim(object@obs)[2] != 1) {
                errors = c(errors, "obs must be a 2-dimensional FLPar with exactly 3 parameters and 1 iteration")
              }
             
             # Check parameter names
             validNames = c("r", "virgin", "bmsy", "msy", "fmsy", "fcrash")
             obsNames = dimnames(object@obs)$params
             if (is.null(obsNames) || length(obsNames) != 3) {
               errors = c(errors, "obs must have exactly 3 named parameters")
             } else if (!all(obsNames %in% validNames)) {
               errors = c(errors, paste("Invalid parameter names in obs. Must be one of:", paste(validNames, collapse = ", ")))
             } else if (length(unique(obsNames)) != 3) {
               errors = c(errors, "Parameter names in obs must be unique")
             }
             
             # Check for valid numeric values
             obsValues = c(object@obs)
             if (any(is.na(obsValues)) || any(is.infinite(obsValues))) {
               errors = c(errors, "obs contains NA or infinite values")
             }
             
             # Check for positive values where required
             if ("r" %in% obsNames && obsValues[obsNames == "r"] <= 0) {
               errors = c(errors, "r (intrinsic growth rate) must be positive")
             }
             if ("virgin" %in% obsNames && obsValues[obsNames == "virgin"] <= 0) {
               errors = c(errors, "virgin (carrying capacity) must be positive")
             }
             if ("bmsy" %in% obsNames && obsValues[obsNames == "bmsy"] <= 0) {
               errors = c(errors, "bmsy must be positive")
             }
             if ("msy" %in% obsNames && obsValues[obsNames == "msy"] <= 0) {
               errors = c(errors, "msy must be positive")
             }
             if ("fmsy" %in% obsNames && obsValues[obsNames == "fmsy"] <= 0) {
               errors = c(errors, "fmsy must be positive")
             }
             if ("fcrash" %in% obsNames && obsValues[obsNames == "fcrash"] <= 0) {
               errors = c(errors, "fcrash must be positive")
             }
             
             # Check for reasonable value ranges
             if ("r" %in% obsNames && obsValues[obsNames == "r"] > 5) {
               errors = c(errors, "r (intrinsic growth rate) seems very high (> 5 per year)")
             }
             if ("fmsy" %in% obsNames && obsValues[obsNames == "fmsy"] > 10) {
               errors = c(errors, "fmsy seems very high (> 10 per year)")
             }
             if ("fcrash" %in% obsNames && obsValues[obsNames == "fcrash"] > 10) {
               errors = c(errors, "fcrash seems very high (> 10 per year)")
             }
             
             # Check logical relationships
             if (all(c("bmsy", "virgin") %in% obsNames)) {
               if (obsValues[obsNames == "bmsy"] >= obsValues[obsNames == "virgin"]) {
                 errors = c(errors, "bmsy cannot be greater than or equal to virgin biomass")
               }
             }
           }
           
                       # Validate comparison slot
            if (!inherits(object@comparison, "FLPar")) {
              errors = c(errors, "comparison slot must be an FLPar object")
            } else {
              # Check if comparison is empty (not yet estimated)
              if (all(is.na(object@comparison))) {
                # Empty comparison is valid - parameters not yet estimated
              } else {
                # Check comparison dimensions
                if (length(dim(object@comparison)) != 3 || dim(object@comparison)[1] != 8 || dim(object@comparison)[2] != 2 || dim(object@comparison)[3] != 1) {
                  errors = c(errors, "comparison must be a 3-dimensional FLPar with 8 parameters, 2 methods, and 1 iteration")
                }
                
                # Check parameter names
                expectedParams = c("r", "virgin", "p", "bmsy", "msy", "fmsy", "fcrash", "bmsy/virgin")
                compParams = dimnames(object@comparison)$params
                if (!identical(compParams, expectedParams)) {
                  errors = c(errors, "comparison must have parameters: r, virgin, p, bmsy, msy, fmsy, fcrash, bmsy/virgin")
                }
                
                # Check method names
                expectedMethods = c("age", "biomass")
                compMethods = dimnames(object@comparison)$method
                if (!identical(compMethods, expectedMethods)) {
                  errors = c(errors, "comparison must have methods: age, biomass")
                }
              }
            }
           
           if (length(errors) == 0) TRUE else errors
         }
)

#' Constructor for pellatParams
#'
#' Creates a pellatParams object from a named vector of observed values.
#' Converts the input to an FLPar object for proper unit handling.
#'
#' @param obs Named numeric vector of observed reference values (length 3).
#'   Names must be one of: r, virgin, bmsy, msy, fmsy, fcrash.
#'   Rate parameters (r, fmsy, fcrash) must be instantaneous rates per year.
#'   Biomass and yield parameters can be in any consistent units.
#' @param units Character vector specifying units for each parameter (optional).
#'   If provided, should match the length and order of obs.
#' @return A pellatParams object.
#' @examples
#' pt = pellatParams(c(r = 0.5, fmsy = 0.3, fcrash = 0.25))
#' pt = pellatParams(c(msy = 500, fmsy = 0.3, virgin = 10000),
#'                              units = c("tonnes/year", "per year", "tonnes"))
pellatParams = function(obs, units = NULL) {
  # Validate input
  if (!is.numeric(obs) || length(obs) != 3) {
    stop("obs must be a numeric vector of length 3")
  }
  
  if (is.null(names(obs)) || any(names(obs) == "")) {
    stop("obs must be a named vector")
  }
  
  validNames = c("r", "virgin", "bmsy", "msy", "fmsy", "fcrash")
  if (!all(names(obs) %in% validNames)) {
    stop("Invalid parameter names. Must be one of: ", paste(validNames, collapse = ", "))
  }
  
  if (length(unique(names(obs))) != length(obs)) {
    stop("Parameter names must be unique")
  }
  
  # Create FLPar object for obs
  obsFLPar = FLPar(obs, params = names(obs))
  
  # Set units if provided
  if (!is.null(units)) {
    if (length(units) != length(obs)) {
      stop("units must have the same length as obs")
    }
    units(obsFLPar) = units
  } else {
    # Set default units based on parameter type
    defaultUnits = sapply(names(obs), function(param) {
      if (param %in% c("r", "fmsy", "fcrash")) {
        "per year"
      } else if (param %in% c("virgin", "bmsy")) {
        "tonnes"
      } else if (param == "msy") {
        "tonnes/year"
      } else {
        "NA"
      }
    })
    units(obsFLPar) = defaultUnits
  }
  
  # Create empty FLPar for comparison
  emptyComparison = FLPar(array(NA, dim = c(8, 2), 
                                dimnames = list(params = c("r", "virgin", "p", "bmsy", "msy", "fmsy", "fcrash", "bmsy/virgin"),
                                              method = c("age", "biomass"))))
  
  new("pellatParams", 
      obs = obsFLPar, 
      comparison = emptyComparison)
}

#' Estimate Pella-Tomlinson Parameters
#'
#' Estimate Pella-Tomlinson model (r, virgin, p) parameters from any three reference points.
#' Returns estimated parameters and reference points in the comparison slot, including fcrash and fmsy.
#'
#' @param object An object of class "pellatParams".
#' @return The input object with the comparison slot updated with estimates and a comparison table.
#' @details
#' You can supply any three of r (intrinsic growth rate), virgin (carrying capacity), bmsy (biomass at msy),
#' msy (maximum sustainable yield), fmsy (fishing mortality at msy), or fcrash (crash fishing mortality).
#' 
#' Input requirements:
#' - Rate parameters (r, fmsy, fcrash) must be instantaneous rates per year
#' - Biomass parameters (virgin, bmsy) can be in any consistent biomass units
#' - Yield parameter (msy) can be in any consistent yield units per year
#' 
#' For the Pella-Tomlinson model: fmsy = r * p / (p + 1) and fcrash = r
#' @examples
#' pt = pellatParams(c(r = 0.5, fmsy = 0.3, fcrash = 0.25))
#' pt = estimateParams(pt)
#' pt@comparison
#' pt@comparison["r", "biomass"]
setGeneric("estimateParams", function(object) standardGeneric("estimateParams"))

setMethod("estimateParams", "pellatParams",
          function(object) {
            if (!requireNamespace("rootSolve", quietly = TRUE)) {
              stop("Package 'rootSolve' is required.")
            }
            
            # Extract observed values from FLPar
            obs = c(object@obs)
            names(obs) = dimnames(object@obs)$params
            
            # Validate F values are reasonable rates (per year)
            fParams = c("fmsy", "fcrash")
            for (param in fParams) {
              if (param %in% names(obs)) {
                if (obs[param] <= 0) {
                  stop(paste(param, "must be positive (fishing mortality rate)"))
                }
                if (obs[param] > 10) {
                  warning(paste(param, "seems very high for a fishing mortality rate:", obs[param]))
                }
              }
            }
            
            # Validate r (intrinsic growth rate) is reasonable
            if ("r" %in% names(obs)) {
              if (obs["r"] <= 0) {
                stop("r (intrinsic growth rate) must be positive")
              }
              if (obs["r"] > 5) {
                warning("r (intrinsic growth rate) seems very high:", obs["r"])
              }
            }
            
            objFun = function(x) {
              r = x[1]; virgin = x[2]; p = x[3]
              
              # Check for invalid parameters
              if (r <= 0 || virgin <= 0 || abs(p + 1) < 1e-10) {
                # Return large values to guide optimization away from invalid regions
                return(rep(1e6, length(names(obs))))
              }
              
              # Handle special case where p is very close to 0 (Fox model)
              if (abs(p) < 1e-10) {
                bmsy = virgin * exp(-1)  # Fox model: bmsy = virgin/e
                msy = r * virgin * exp(-1) * (1 - exp(-1))  # Fox model msy
                fmsy = r * exp(-1)  # Fox model fmsy
              } else {
                # Regular Pella-Tomlinson model
                bmsy = virgin * (1/(p+1))^(1/p)
                msy = r * virgin * ((1/(p+1))^(1/p) - (1/(p+1))^((p+1)/p))
                fmsy = r * p / (p+1)
              }
              
              # For Pella-Tomlinson, fcrash ~ r
              fcrash = r
              
              vals = c(r = r, virgin = virgin, bmsy = bmsy,
                        msy = msy, fmsy = fmsy, fcrash = fcrash)
              
              # Check for invalid results
              if (any(is.na(vals)) || any(is.infinite(vals))) {
                return(rep(1e6, length(names(obs))))
              }
              
              sapply(names(obs), function(name) vals[name] - obs[name])
            }
            
            # Provide reasonable start guesses
            startVec = c(
              r = if ("r" %in% names(obs)) obs["r"] else 0.4,
              virgin = if ("virgin" %in% names(obs)) obs["virgin"] else 1000,
              p = 1
            )
            
            # Handle different parameter combinations
            if (all(c("r", "fmsy", "fcrash") %in% names(obs))) {
              # Case 1: r, fmsy, fcrash
              r = obs["r"]
              fmsy = obs["fmsy"]
              fcrash = obs["fcrash"]
              
              # Use rFmsy2p function to calculate p from r and fmsy
              p = rFmsy2p(r, fmsy)
              
              # Check if p calculation was successful
              if (is.na(p)) {
                stop("Failed to calculate p for r = ", r, " and fmsy = ", fmsy)
              }
              
              # From fcrash = r, we can verify
              if (abs(fcrash - r) > 1e-4) {
                warning("fcrash should equal r for Pella-Tomlinson model")
              }
              
              # We need virgin, but we don't have enough information
              virgin = 1000
              
              sol = list(root = c(r, virgin, p))
              
            } else if (all(c("msy", "fmsy", "virgin") %in% names(obs))) {
              # Case 2: msy, fmsy, virgin
              msy = obs["msy"]
              fmsy = obs["fmsy"]
              virgin = obs["virgin"]
              
              # Calculate bmsy from msy and fmsy: bmsy = msy / fmsy
              bmsy = msy / fmsy
              
              # Use bmsyK2p function to calculate p from bmsy/virgin ratio
              bmsyVirginRatio = bmsy / virgin
              p = bmsyK2p(bmsyVirginRatio)
              
              # Check if p calculation was successful
              if (is.na(p)) {
                stop("Failed to calculate p for bmsy/virgin ratio: ", bmsyVirginRatio)
              }
              
              # Now solve for r: r = fmsy * (p+1) / p
              # This is the inverse of fmsy = r * p / (p+1)
              if (abs(p) < 1e-10) {
                stop("p parameter too close to zero")
              }
              r = fmsy * (p + 1) / p
              
              sol = list(root = c(r, virgin, p))
              
            } else if (all(c("bmsy", "fmsy", "virgin") %in% names(obs))) {
              # Case 3: bmsy, fmsy, virgin
              bmsy = obs["bmsy"]
              fmsy = obs["fmsy"]
              virgin = obs["virgin"]
              
              # Use bmsyK2p function to calculate p from bmsy/virgin ratio
              bmsyVirginRatio = bmsy / virgin
              p = bmsyK2p(bmsyVirginRatio)
              
              # Check if p calculation was successful
              if (is.na(p)) {
                stop("Failed to calculate p for bmsy/virgin ratio: ", bmsyVirginRatio)
              }
              
              # Now solve for r: r = fmsy * (p+1) / p
              # This is the inverse of fmsy = r * p / (p+1)
              if (abs(p) < 1e-10) {
                stop("p parameter too close to zero")
              }
              r = fmsy * (p + 1) / p
              
              sol = list(root = c(r, virgin, p))
              
            } else if (all(c("msy", "bmsy", "virgin") %in% names(obs))) {
              # Case 4: msy, bmsy, virgin
              msy = obs["msy"]
              bmsy = obs["bmsy"]
              virgin = obs["virgin"]
              
              # Check if the values are reasonable
              if (bmsy >= virgin) {
                stop("BMSY cannot be greater than or equal to virgin biomass")
              }
              
              # Use bmsyK2p function to calculate p from bmsy/virgin ratio
              bmsyVirginRatio = bmsy / virgin
              p = bmsyK2p(bmsyVirginRatio)
              
              # Check if p calculation was successful
              if (is.na(p)) {
                stop("Failed to calculate p for bmsy/virgin ratio: ", bmsyVirginRatio)
              }
              
              # Check if the solution is reasonable
              if (abs(p + 1) < 1e-10) {
                stop("Solution found has p = -1, which is invalid for Pella-Tomlinson model")
              }
              
              # Now solve for r: r = msy / (virgin * ((1/(p+1))^(1/p) - (1/(p+1))^((p+1)/p)))
              # This is the correct formula for MSY in Pella-Tomlinson model
              if (abs(p) < 1e-10) {
                # Fox model
                r = msy / (virgin * exp(-1) * (1 - exp(-1)))
              } else {
                # For Pella-Tomlinson: MSY = r * virgin * ((1/(p+1))^(1/p) - (1/(p+1))^((p+1)/p))
                # Therefore: r = MSY / (virgin * ((1/(p+1))^(1/p) - (1/(p+1))^((p+1)/p)))
                denominator = virgin * ((1/(p+1))^(1/p) - (1/(p+1))^((p+1)/p))
                r = msy / denominator
              }
              
              # Check if r is positive
              if (r <= 0) {
                stop(paste("Calculated r (intrinsic growth rate) is not positive. p =", p, ", r =", r))
              }
              
              sol = list(root = c(r, virgin, p))
              
            } else {
              # For other combinations, use a simple optimization approach
              warning("Direct solution not available for this parameter combination")
              sol = list(root = startVec)
            }
            
            r = sol$root[1]; virgin = sol$root[2]; p = sol$root[3]
            
            # Validate the solution
            if (r <= 0 || virgin <= 0 || abs(p + 1) < 1e-10) {
              stop("Invalid solution found: r, virgin must be positive and p must not be -1")
            }
            
            # Handle special case where p is very close to 0 (Fox model)
            if (abs(p) < 1e-10) {
              bmsy = virgin * exp(-1)  # Fox model: bmsy = virgin/e
              msy = r * virgin * exp(-1) * (1 - exp(-1))  # Fox model msy
              fmsy = r * exp(-1)  # Fox model fmsy
              fcrash = r  # Fox model: fcrash = r
            } else {
              # Regular Pella-Tomlinson model
              bmsy = virgin * (1/(p+1))^(1/p)
              msy = r*virgin*((1/(p+1))^(1/p) - (1/(p+1))^((p+1)/p))
              fmsy = r * p / (p+1)
              fcrash = r  # Pella-Tomlinson model: fcrash = r (instantaneous rate)
            }
            
            # Create a more robust comparison matrix
            allParams = c("r", "virgin", "p", "bmsy", "msy", "fmsy", "fcrash", "bmsy/virgin")
            paramValues = c(r, virgin, p, bmsy, msy, fmsy, fcrash, bmsy/virgin)
            
            # Create the comparison matrix
            ageBasedValues = sapply(allParams, function(param) {
              if (param %in% names(obs)) {
                obs[param]
              } else if (param == "bmsy/virgin" && all(c("bmsy", "virgin") %in% names(obs))) {
                obs["bmsy"] / obs["virgin"]
              } else {
                NA
              }
            })
            
            # Create FLPar for comparison
            comparisonData = matrix(
              c(ageBasedValues, paramValues),
              nrow = length(allParams),
              ncol = 2,
              dimnames = list(params = allParams, method = c("age", "biomass"))
            )
            
            comparisonFLPar = FLPar(comparisonData)
            
            # Set appropriate units for comparison
            comparisonUnits = sapply(allParams, function(param) {
              if (param %in% c("r", "fmsy", "fcrash")) {
                "per year"
              } else if (param %in% c("virgin", "bmsy")) {
                units(object@obs)[names(obs) == "virgin"][1]  # Use virgin units
              } else if (param == "msy") {
                units(object@obs)[names(obs) == "msy"][1]    # Use msy units
              } else if (param == "p") {
                "dimensionless"
              } else if (param == "bmsy/virgin") {
                "dimensionless"
              } else {
                "NA"
              }
            })
            
            units(comparisonFLPar) = comparisonUnits
            
            # Store results in comparison FLPar
            object@comparison = comparisonFLPar
            
            object
          }
)

#' Print method for pellatParams
#'
#' @param x A pellatParams object
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the object
#' @examples
#' pt = pellatParams(c(msy = 500, fmsy = 0.3, virgin = 10000))
#' print(pt)
setMethod("print", "pellatParams",
          function(x, ...) {
            cat("Pella-Tomlinson Estimator\n")
            cat("========================\n\n")
            
            cat("Observed parameters:\n")
            print(x@obs)
            cat("\n")
            
            if (all(!is.na(x@comparison))) {
              cat("Estimated parameters:\n")
              cat(sprintf("  r (intrinsic growth rate, per year): %.6g\n", x@comparison["r", "biomass"]))
              cat(sprintf("  virgin (carrying capacity): %.6g\n", x@comparison["virgin", "biomass"]))
              cat(sprintf("  p (shape parameter): %.6g\n", x@comparison["p", "biomass"]))
              cat(sprintf("  bmsy: %.6g\n", x@comparison["bmsy", "biomass"]))
              cat(sprintf("  msy: %.6g\n", x@comparison["msy", "biomass"]))
              cat(sprintf("  fmsy (fishing mortality at MSY, per year): %.6g\n", x@comparison["fmsy", "biomass"]))
              cat(sprintf("  fcrash (crash fishing mortality, per year): %.6g\n", x@comparison["fcrash", "biomass"]))
            } else {
              cat("No parameters estimated yet. Use estimateParams() to fit the model.\n")
            }
            
            invisible(x)
          })

#' Summary method for pellatParams
#'
#' @param object A pellatParams object
#' @param ... Additional arguments (not used)
#' @return A summary object
#' @examples
#' pt = pellatParams(c(msy = 500, fmsy = 0.3, virgin = 10000))
#' pt = estimateParams(pt)
#' summary(pt)
setMethod("summary", "pellatParams",
          function(object, ...) {
            if (all(is.na(object@comparison))) {
              stop("No parameters estimated yet. Use estimateParams() to fit the model.")
            }
            
            # Create summary list
            summaryObj = list(
              observed = object@obs,
              estimated = list(
                r = c(object@comparison["r", "biomass"]),
                virgin = c(object@comparison["virgin", "biomass"]),
                p = c(object@comparison["p", "biomass"]),
                bmsy = c(object@comparison["bmsy", "biomass"]),
                msy = c(object@comparison["msy", "biomass"]),
                fmsy = c(object@comparison["fmsy", "biomass"]),
                fcrash = c(object@comparison["fcrash", "biomass"])
              ),
              comparison = object@comparison,
              modelInfo = list(
                bmsyVirginRatio = c(object@comparison["bmsy", "biomass"]) / c(object@comparison["virgin", "biomass"]),
                isFoxModel = abs(c(object@comparison["p", "biomass"])) < 1e-10
              )
            )
            
            # Add class for custom printing
            class(summaryObj) = "summary.pellatParams"
            summaryObj
          })

#' Print method for summary.pellatParams
#'
#' @param x A summary.pellatParams object
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the object
print.summary.pellatParams = function(x, ...) {
  cat("Pella-Tomlinson Model Summary\n")
  cat("=============================\n\n")
  
  cat("Observed Parameters:\n")
  print(x$observed)
  cat("\n")
  
  cat("Estimated Parameters:\n")
  cat(sprintf("  r (intrinsic growth rate, per year): %.6g\n", x$estimated$r))
  cat(sprintf("  virgin (carrying capacity): %.6g\n", x$estimated$virgin))
  cat(sprintf("  p (shape parameter): %.6g\n", x$estimated$p))
  cat(sprintf("  bmsy: %.6g\n", x$estimated$bmsy))
  cat(sprintf("  msy: %.6g\n", x$estimated$msy))
  cat(sprintf("  fmsy (fishing mortality at MSY, per year): %.6g\n", x$estimated$fmsy))
  cat(sprintf("  fcrash (crash fishing mortality, per year): %.6g\n", x$estimated$fcrash))
  cat("\n")
  
  cat("Model Characteristics:\n")
  cat(sprintf("  bmsy/virgin ratio: %.6g\n", x$modelInfo$bmsyVirginRatio))
  if (x$modelInfo$isFoxModel) {
    cat("  Model type: Fox (p ≈ 0)\n")
  } else if (x$estimated$p < 0) {
    cat("  Model type: Depensatory (p < 0)\n")
  } else if (x$estimated$p > 0) {
    cat("  Model type: Compensatory (p > 0)\n")
  }
  cat("\n")
  
  cat("Parameter Comparison:\n")
  print(x$comparison)
  
  invisible(x)
}

#' Get comparison matrix
#'
#' Returns a clean comparison matrix between observed and estimated parameters.
#'
#' @param object A pellatParams object
#' @param digits Number of digits to display (default: 6)
#' @return A formatted matrix
#' @examples
#' pt = pellatParams(c(msy = 500, fmsy = 0.3, virgin = 10000))
#' pt = estimateParams(pt)
#' getComparisonTable(pt)
setGeneric("getComparisonTable", function(object, digits = 6) standardGeneric("getComparisonTable"))

setMethod("getComparisonTable", "pellatParams",
          function(object, digits = 6) {
            if (all(is.na(object@comparison))) {
              stop("No parameters estimated yet. Use estimateParams() to fit the model.")
            }
            
            # Extract data from FLPar and format
            compMatrix = c(object@comparison)
            dim(compMatrix) = dim(object@comparison)
            dimnames(compMatrix) = dimnames(object@comparison)
            
            compMatrix[is.na(compMatrix)] = "NA"
            compMatrix[compMatrix != "NA"] = sprintf(paste0("%.", digits, "g"),
                                                   as.numeric(compMatrix[compMatrix != "NA"]))
            
            compMatrix
          })

#' Plot method for pellatParams
#'
#' Creates a ggplot2 plot showing the production function and reference points.
#'
#' @param x A pellatParams object
#' @param y Not used (required for S4 method signature)
#' @param ... Additional arguments passed to ggplot2 functions
#' @return A ggplot2 object
#' @examples
#' pt = pellatParams(c(msy = 500, fmsy = 0.3, virgin = 10000))
#' pt = estimateParams(pt)
#' plot(pt)
#' 
#' # Save plot to file
#' p = plot(pt)
#' ggplot2::ggsave("production_function.png", p, width = 10, height = 6)
setMethod("plot", "pellatParams",
          function(x, y, ...) {
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
              stop("Package 'ggplot2' is required for plotting.")
            }
            
            if (all(is.na(x@comparison))) {
              stop("No parameters estimated yet. Use estimateParams() to fit the model.")
            }
            
            # Extract parameters
            r = c(x@comparison["r", "biomass"])
            virgin = c(x@comparison["virgin", "biomass"])
            p = c(x@comparison["p", "biomass"])
            
            # Create biomass sequence for plotting
            biomass = seq(0, virgin * 1.1, length.out = 1000)
            
            # Calculate production function
            if (abs(p) < 1e-10) {
              # Fox model
              production = r * biomass * log(virgin / biomass)
            } else {
              # Pella-Tomlinson model
              production = r * biomass * (1 - (biomass / virgin)^p) / p
            }
            
            # Filter out negative production values
            validIndices = production >= 0
            biomass = biomass[validIndices]
            production = production[validIndices]
            
            # Create data frame for plotting
            plotData = data.frame(
              biomass = biomass,
              production = production
            )
            
            # Extract reference points for plotting
            refpts = data.frame(
              biomass = c(
                virgin,  # virgin biomass
                c(x@comparison["bmsy", "biomass"])  # bmsy
              ),
              production = c(
                0,  # virgin biomass has zero production
                c(x@comparison["msy", "biomass"])  # MSY
              ),
              point = c("Virgin", "MSY"),
              type = c("reference", "reference")
            )
            
            # Add observed values if they exist
            obsData = data.frame(
              biomass = numeric(0),
              production = numeric(0),
              point = character(0),
              type = character(0)
            )
            
            # Check for observed biomass parameters
            if ("virgin" %in% dimnames(x@obs)$params) {
              obsVirgin = c(x@obs["virgin"])
              obsData = rbind(obsData, data.frame(
                biomass = obsVirgin,
                production = 0,
                point = "Observed Virgin",
                type = "observed"
              ))
            }
            
            if ("bmsy" %in% dimnames(x@obs)$params) {
              obsbmsy = c(x@obs["bmsy"])
              # Calculate corresponding production for observed bmsy
              if (abs(p) < 1e-10) {
                obsProduction = r * obsbmsy * log(virgin / obsbmsy)
              } else {
                obsProduction = r * obsbmsy * (1 - (obsbmsy / virgin)^p) / p
              }
              obsData = rbind(obsData, data.frame(
                biomass = obsbmsy,
                production = obsProduction,
                point = "Observed bmsy",
                type = "observed"
              ))
            }
            
            # Check for observed MSY
            if ("msy" %in% dimnames(x@obs)$params) {
              obsMsy = c(x@obs["msy"])
              # Find biomass corresponding to observed MSY
              if (abs(p) < 1e-10) {
                # Fox model: MSY occurs at B = virgin/e
                obsbmsy = virgin * exp(-1)
              } else {
                # Pella-Tomlinson: MSY occurs at B = virgin * (1/(p+1))^(1/p)
                obsbmsy = virgin * (1/(p+1))^(1/p)
              }
              obsData = rbind(obsData, data.frame(
                biomass = obsbmsy,
                production = obsMsy,
                point = "Observed MSY",
                type = "observed"
              ))
            }
            
                         # Create the plot
             p = ggplot2::ggplot(plotData, ggplot2::aes(x = biomass, y = production)) +
               ggplot2::geom_line(color = "blue", linewidth = 1) +
                               # Add tangent lines at origin for r and fcrash
                ggplot2::geom_abline(intercept = 0, slope = r, color = "grey", linetype = "solid", alpha = 0.7) +
                ggplot2::geom_abline(intercept = 0, slope = c(x@comparison["fcrash", "biomass"]), color = "red", linetype = "dashed", alpha = 0.7) +
               # Add vertical reference lines
               ggplot2::geom_vline(xintercept = c(x@comparison["bmsy", "biomass"]), color = "red", linetype = "dashed", alpha = 0.7) +
               ggplot2::geom_vline(xintercept = virgin * 0.1, color = "gray", linetype = "dotted", alpha = 0.5) +
               # Plot observed points first (large green dots)
               ggplot2::geom_point(data = obsData, 
                                  ggplot2::aes(color = type), 
                                  size = 4, shape = 16) +
               # Plot reference points second (small red dots) so they overlay on top
               ggplot2::geom_point(data = refpts, 
                                  ggplot2::aes(color = type), 
                                  size = 2, shape = 16) +
               ggplot2::scale_color_manual(values = c("reference" = "red", "observed" = "green")) +
                             ggplot2::labs(
                 title = "Pella-Tomlinson Production Function",
                 subtitle = sprintf("r = %.3f, virgin = %.0f, p = %.3f", r, virgin, p),
                 x = "Biomass",
                 y = "Production",
                 color = "Point Type",
                                   caption = "Grey solid line: slope at origin (r)\nRed dashed line: fcrash slope"
               ) +
              ggplot2::theme_minimal() +
              ggplot2::theme(
                legend.position = "bottom",
                plot.title = ggplot2::element_text(hjust = 0.5),
                plot.subtitle = ggplot2::element_text(hjust = 0.5)
              )
            
            # Add text labels for key points
            p = p + ggplot2::geom_text(
              data = refpts,
              ggplot2::aes(label = point),
              hjust = -0.1,
              vjust = 0.5,
              size = 3
            )
            
            if (nrow(obsData) > 0) {
              p = p + ggplot2::geom_text(
                data = obsData,
                ggplot2::aes(label = point),
                hjust = -0.1,
                vjust = 0.5,
                size = 3,
                color = "green"
              )
            }
            
            p
          })

#' Validate pellatParams object
#'
#' Validates a pellatParams object and returns detailed error messages if invalid.
#'
#' @param object A pellatParams object
#' @return TRUE if valid, character vector of error messages if invalid
#' @examples
#' pt = pellatParams(c(msy = 500, fmsy = 0.3, virgin = 10000))
#' validatePellatParams(pt)
validatePellatParams = function(object) {
  if (!inherits(object, "pellatParams")) {
    return("Object must be of class 'pellatParams'")
  }
  
  # Use the validity function
  validityResult = validObject(object, test = TRUE)
  if (is.logical(validityResult) && validityResult) {
    return(TRUE)
  } else {
    return(validityResult)
  }
}

#' Check if pellatParams object is valid
#'
#' Simple boolean check for object validity.
#'
#' @param object A pellatParams object
#' @return TRUE if valid, FALSE otherwise
#' @examples
#' pt = pellatParams(c(msy = 500, fmsy = 0.3, virgin = 10000))
#' isValidPellatParams(pt)
isValidPellatParams = function(object) {
  result = validatePellatParams(object)
  return(is.logical(result) && result)
}

# # Example usage:
# pt = pellatParams(c(msy=500,fmsy=0.2,virgin=23000))
# pt = estimateParams(pt)
# 
# getComparisonTable(pt) 
# plot(pt)         
