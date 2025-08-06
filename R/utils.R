#' Calculate Tmax using different methods
#' 
#' @description Calculate maximum rebuilding time (Tmax) using various methods
#' as described in NS1 guidelines
#' 
#' @param object An FLR object with life history parameters
#' @param method Method to use: "generation" (Tmin + generation time), 
#'               "fishing" (time to rebuild at 75% MFMT), or "multiply" (Tmin * 2)
#' @param tmin Minimum rebuilding time
#' @param mfmt Maximum fishing mortality threshold
#' @return Tmax value
#' @export
calculateTmax = function(object, method = "generation", tmin = NULL, mfmt = NULL) {
  
  if (method == "generation") {
    # Method 1: Tmin + one generation time
    if (is.null(tmin)) {
      stop("tmin must be provided for generation method")
    }
    gt = gt(object)  # generation time
    return(tmin + gt)
    
  } else if (method == "fishing") {
    # Method 2: Time to rebuild at 75% MFMT
    if (is.null(mfmt)) {
      stop("mfmt must be provided for fishing method")
    }
    target_f = mfmt * 0.75
    # This would need implementation based on specific stock dynamics
    stop("Fishing method not yet implemented")
    
  } else if (method == "multiply") {
    # Method 3: Tmin * 2
    if (is.null(tmin)) {
      stop("tmin must be provided for multiply method")
    }
    return(tmin * 2)
    
  } else {
    stop("Method must be one of: 'generation', 'fishing', 'multiply'")
  }
}

#' Calculate recovery time based on population growth rate
#' 
#' @description Calculate time to recovery using population growth rate
#' 
#' @param initial_biomass Initial biomass relative to target
#' @param target_biomass Target biomass (e.g., BMSY)
#' @param growth_rate Population growth rate
#' @return Time to recovery in years
#' @export
calculateRecoveryTime = function(initial_biomass, target_biomass, growth_rate) {
  if (initial_biomass >= target_biomass) {
    return(0)
  }
  
  # Using exponential growth model: B(t) = B0 * exp(r*t)
  # Solve for t: t = ln(Btarget/B0) / r
  recovery_time = log(target_biomass / initial_biomass) / growth_rate
  
  return(max(0, recovery_time))
}

#' Plot rebuilding trajectories
#' 
#' @description Create a ggplot of rebuilding trajectories
#' 
#' @param rebuild_data Data frame with rebuilding trajectories
#' @param target_line Target biomass level (default = 1 for BMSY)
#' @param title Plot title
#' @return ggplot object
#' @export
plotRebuildTrajectories = function(rebuild_data, target_line = 1, title = "Rebuilding Trajectories") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  ggplot2::ggplot(rebuild_data, ggplot2::aes(x = year, y = data, group = iter)) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::geom_hline(yintercept = target_line, color = "red", linetype = "dashed") +
    ggplot2::labs(title = title, x = "Year", y = "Biomass relative to target") +
    ggplot2::theme_minimal()
} 

#' Find Age with Maximum Value
#' 
#' @description 
#' Identifies which age has the highest value for each year and iteration in an FLQuant object.
#' Excludes the plusgroup (last age) from the calculation.
#' 
#' @param object An FLQuant object with age, year, and iter dimensions
#' @param na.rm Logical. Should NA values be removed? Default is TRUE
#' @param ... Additional arguments (not currently used)
#' 
#' @return An FLQuant object with 1 quant dimension, containing the age indices
#' where the maximum value occurs for each year and iteration
#' 
#' @examples
#' \dontrun{
#' # Example with sample data
#' data(ple4)
#' max_age = ageMax(ssb(ple4))
#' 
#' # Get the actual ages (not indices)
#' ages = an(dimnames(ple4)$age)
#' max_age_actual = ages[c(max_age)]
#' }
#' 
#' @rdname ageMax
#' @export
setGeneric("ageMax", function(object, na.rm = TRUE, ...)
  standardGeneric("ageMax"))

#' @rdname ageMax
#' @export
setMethod("ageMax", signature(object = "FLQuant"),
          function(object, na.rm = TRUE, ...) {
            # Get dimensions
            dmns = dimnames(object)
            nyear = dim(object)[2]
            niter = dim(object)[6]
            
            # Create output FLQuant with 1 quant dimension
            # Remove the age dimension and replace with quant
            result_dmns = dmns
            result_dmns$quant = "all"
            result_dmns$age = NULL
            result = FLQuant(NA, dimnames = result_dmns)
            
            # For each year and iteration, find the age with maximum value (excluding plusgroup)
            for (year in 1:nyear) {
              for (iter in 1:niter) {
                # Extract the age vector for this year and iteration (excluding plusgroup)
                age_vec = c(object[-dim(object)[1], year, , , , iter])
                
                if (na.rm) {
                  # Remove NA values
                  valid_idx = !is.na(age_vec)
                  if (sum(valid_idx) > 0) {
                    age_vec = age_vec[valid_idx]
                    max_idx = which.max(age_vec)
                    # Convert back to original age index (excluding plusgroup)
                    result[1, year, , , , iter] = which(valid_idx)[max_idx]
                  }
                } else {
                  # Include NA values
                  max_idx = which.max(age_vec)
                  result[1, year, , , , iter] = max_idx
                }
              }
            }
            
            return(result)
          })

#' Find Age with Maximum Value (Alternative Implementation)
#' 
#' @description 
#' Alternative implementation using qapply for better performance with large datasets.
#' Excludes the plusgroup (last age) from the calculation.
#' 
#' @param object An FLQuant object with age, year, and iter dimensions
#' @param na.rm Logical. Should NA values be removed? Default is TRUE
#' @param ... Additional arguments (not currently used)
#' 
#' @return An FLQuant object with 1 quant dimension, containing the age indices
#' where the maximum value occurs for each year and iteration
#' 
#' @rdname ageMaxQ
#' @export
setGeneric("ageMaxQ", function(object, na.rm = TRUE, ...)
  standardGeneric("ageMaxQ"))

#' @rdname ageMaxQ
#' @export
setMethod("ageMaxQ", signature(object = "FLQuant"),
          function(object, na.rm = TRUE, ...) {
            # Use qapply to find max age for each iteration
            result = qapply(object, function(x, na.rm) {
              # Get dimensions
              dmns = dimnames(x)
              nyear = dim(x)[2]
              
              # Create output FLQuant with 1 quant dimension
              # Remove the age dimension and replace with quant
              result_dmns = dmns
              result_dmns$quant = "1"
              result_dmns$age = NULL
              result = FLQuant(NA, dimnames = result_dmns)
              
              # For each year, find the age with maximum value (excluding plusgroup)
              for (year in 1:nyear) {
                age_vec = c(x[-dim(x)[1], year])
                
                if (na.rm) {
                  # Remove NA values
                  valid_idx = !is.na(age_vec)
                  if (sum(valid_idx) > 0) {
                    age_vec = age_vec[valid_idx]
                    max_idx = which.max(age_vec)
                    # Convert back to original age index (excluding plusgroup)
                    result[1, year] = which(valid_idx)[max_idx]
                  }
                } else {
                  # Include NA values
                  max_idx = which.max(age_vec)
                  result[1, year] = max_idx
                }
              }
              
              return(result)
            }, na.rm = na.rm)
            
            return(result)
          })




