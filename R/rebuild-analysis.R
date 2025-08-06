#' Rebuild Analysis Functions for FLRebuild Package
#' 
#' This file contains functions for performing rebuild analysis on FLBRP objects
#' with parallel processing capabilities.

#' Calculate Rebuild Times from Multiple Scenarios
#' 
#' Calculate rebuild times from multiple scenarios using time series data,
#' reference point SSB data, and generation time data.
#' 
#' @param rebuild_ts Rebuild time series data
#' @param refpt_ssb Reference point SSB data
#' @param gen_time Generation time data
#' @param scenarios Character vector of scenario names
#' @param srrs Character vector of SRR names
#' @return Rebuild time calculations
#' @export
calculateRebuildTimes = function(rebuild_ts, refpt_ssb, gen_time, scenarios, srrs) {
  results = list()
  
  for (scenario in scenarios) {
    for (srr in srrs) {
      scenario_name = paste0(scenario, "_", srr)
      
      # Filter data for current scenario
      ts_filter = rebuild_ts$scenario == scenario & rebuild_ts$srr == srr
      ssb_filter = refpt_ssb$scenario == scenario & refpt_ssb$srr == srr
      gen_filter = gen_time$scenario == scenario & gen_time$srr == srr
      
      if (any(ts_filter) && any(ssb_filter) && any(gen_filter)) {
        ts_data = rebuild_ts[ts_filter, ]
        ssb_data = refpt_ssb[ssb_filter, ]
        gen_data = gen_time[gen_filter, ]
        
        # Calculate rebuild times
        rebuild_calc = data.frame(
          scenario = scenario,
          srr = srr,
          .id = unique(ts_data$.id),
          stringsAsFactors = FALSE
        )
        
        # Add rebuild time calculations
        rebuild_calc$rebuild_time = tryIt({
          rebuildTime(ts_data, ssb_data, gen_data)
        }, "Error calculating rebuild time")
        
        results[[scenario_name]] = rebuild_calc
      }
    }
  }
  
  # Combine all results
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame(scenario = character(), srr = character(), .id = character(), rebuild_time = numeric()))
  }
}

#' Parallel Rebuild Analysis for Multiple Stocks
#' 
#' Perform parallel rebuild analysis for multiple stocks.
#' 
#' @param flbrp_list List of FLBRP objects
#' @param stock_names Character vector of stock names
#' @param ncores Number of cores to use (default=detectCores()-4)
#' @param rebuild_args List of arguments to pass to rebuild function
#' @return List of rebuild results
#' @export
parallelRebuildAnalysis = function(flbrp_list, stock_names, ncores = NULL, rebuild_args = list()) {
  # Setup parallel processing
  if (is.null(ncores)) {
    ncores = detectCores() - 4
  }
  cl = makeCluster(ncores)
  registerDoParallel(cl)
  
  # Rebuild analysis function for single stock
  rebuildSingleStock = function(stock_name, flbrp_list, rebuild_args) {
    tryCatch({
      flbrp = flbrp_list[[stock_name]]
      if (is.null(flbrp)) return(NULL)
      
      # Set default rebuild arguments
      defaults = list(
        target = 0.3,
        max_year = 50,
        method = "fbar"
      )
      rebuild_args_final = modifyList(defaults, rebuild_args)
      
      # Perform rebuild analysis
      rebuild_result = do.call(rebuild, c(list(flbrp), rebuild_args_final))
      
      return(rebuild_result)
    }, error = function(e) {
      warning(paste("Error in rebuild analysis for", stock_name, ":", e$message))
      return(NULL)
    })
  }
  
  results = foreach(stock_name = stock_names, .packages = c("FLCore", "FLBRP", "FLRebuild")) %dopar% {
    rebuildSingleStock(stock_name, flbrp_list, rebuild_args)
  }
  
  stopCluster(cl)
  names(results) = stock_names
  results = results[!sapply(results, is.null)]
  
  return(results)
}

#' Parallel Reference Point Extraction
#' 
#' Extract reference points from FLBRP objects using parallel processing.
#' 
#' @param flbrp_list List of FLBRP objects
#' @param ncores Number of cores to use (default=detectCores()-4)
#' @param refpts_to_extract Character vector of reference points to extract
#' @return Extracted reference points
#' @export
parallelRefptExtraction = function(flbrp_list, ncores = NULL, refpts_to_extract = c("msy", "f0.1", "spr.30")) {
  # Setup parallel processing
  if (is.null(ncores)) {
    ncores = detectCores() - 4
  }
  cl = makeCluster(ncores)
  registerDoParallel(cl)
  
  # Reference point extraction function for single FLBRP
  extractSingleRefpts = function(id, flbrp_list, refpts_to_extract) {
    tryCatch({
      flbrp = flbrp_list[[id]]
      if (is.null(flbrp)) return(NULL)
      
      # Extract reference points
      refpts_data = data.frame(
        .id = id,
        refpt = refpts_to_extract,
        stringsAsFactors = FALSE
      )
      
      # Add reference point values
      for (refpt in refpts_to_extract) {
        if (refpt %in% dimnames(refpts(flbrp))$refpt) {
          refpts_data$value[refpts_data$refpt == refpt] = 
            as.numeric(refpts(flbrp)[refpt, "ssb", drop = TRUE])
        } else {
          refpts_data$value[refpts_data$refpt == refpt] = NA
        }
      }
      
      return(refpts_data)
    }, error = function(e) {
      warning(paste("Error extracting reference points for", id, ":", e$message))
      return(NULL)
    })
  }
  
  ids = names(flbrp_list)
  results = foreach(id = ids, .packages = c("FLCore", "FLBRP")) %dopar% {
    extractSingleRefpts(id, flbrp_list, refpts_to_extract)
  }
  
  stopCluster(cl)
  names(results) = ids
  results = results[!sapply(results, is.null)]
  
  # Combine all results
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame(.id = character(), refpt = character(), value = numeric()))
  }
}

#' Evaluate Rebuild Potential Across Scenarios
#' 
#' Evaluate rebuild potential across different scenarios.
#' 
#' @param rebuild_time_data Rebuild time data
#' @param gen_time_data Generation time data
#' @param scenarios Character vector of scenario names
#' @param srrs Character vector of SRR names
#' @return Rebuild potential assessment
#' @export
evaluateRebuildPotential = function(rebuild_time_data, gen_time_data, scenarios, srrs) {
  results = list()
  
  for (scenario in scenarios) {
    for (srr in srrs) {
      scenario_name = paste0(scenario, "_", srr)
      
      # Filter data for current scenario
      rebuild_filter = rebuild_time_data$scenario == scenario & rebuild_time_data$srr == srr
      gen_filter = gen_time_data$scenario == scenario & gen_time_data$srr == srr
      
      if (any(rebuild_filter) && any(gen_filter)) {
        rebuild_data = rebuild_time_data[rebuild_filter, ]
        gen_data = gen_time_data[gen_filter, ]
        
        # Calculate rebuild potential metrics
        rebuild_potential = data.frame(
          scenario = scenario,
          srr = srr,
          .id = unique(rebuild_data$.id),
          stringsAsFactors = FALSE
        )
        
        # Add rebuild potential calculations
        rebuild_potential$rebuild_ratio = tryIt({
          rebuild_data$rebuild_time / gen_data$gen_time
        }, NA)
        
        rebuild_potential$rebuild_feasibility = tryIt({
          ifelse(rebuild_potential$rebuild_ratio <= 2, "Feasible", 
                 ifelse(rebuild_potential$rebuild_ratio <= 5, "Challenging", "Difficult"))
        }, "Unknown")
        
        results[[scenario_name]] = rebuild_potential
      }
    }
  }
  
  # Combine all results
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame(scenario = character(), srr = character(), .id = character(), rebuild_ratio = numeric(), rebuild_feasibility = character()))
  }
} 