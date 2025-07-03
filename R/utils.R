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
calculateTmax <- function(object, method = "generation", tmin = NULL, mfmt = NULL) {
  
  if (method == "generation") {
    # Method 1: Tmin + one generation time
    if (is.null(tmin)) {
      stop("tmin must be provided for generation method")
    }
    gt <- gt(object)  # generation time
    return(tmin + gt)
    
  } else if (method == "fishing") {
    # Method 2: Time to rebuild at 75% MFMT
    if (is.null(mfmt)) {
      stop("mfmt must be provided for fishing method")
    }
    target_f <- mfmt * 0.75
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
calculateRecoveryTime <- function(initial_biomass, target_biomass, growth_rate) {
  if (initial_biomass >= target_biomass) {
    return(0)
  }
  
  # Using exponential growth model: B(t) = B0 * exp(r*t)
  # Solve for t: t = ln(Btarget/B0) / r
  recovery_time <- log(target_biomass / initial_biomass) / growth_rate
  
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
plotRebuildTrajectories <- function(rebuild_data, target_line = 1, title = "Rebuilding Trajectories") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }
  
  ggplot2::ggplot(rebuild_data, ggplot2::aes(x = year, y = data, group = iter)) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::geom_hline(yintercept = target_line, color = "red", linetype = "dashed") +
    ggplot2::labs(title = title, x = "Year", y = "Biomass relative to target") +
    ggplot2::theme_minimal()
} 