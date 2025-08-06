#' Reference Point Processing Functions for FLRebuild Package
#' 
#' This file contains functions for extracting and processing reference points
#' from FLBRP objects.

#' Extract Reference Points from FLBRP Objects
#' 
#' Extract reference points from a list of FLBRP objects.
#' 
#' @param flbrp_list List of FLBRP objects
#' @param refpts_to_extract Character vector of reference points to extract (default: c("msy", "f0.1", "spr.30"))
#' @return Data frame of extracted reference points
#' @export
extractRefpts = function(flbrp_list, refpts_to_extract = c("msy", "f0.1", "spr.30")) {
  results = list()
  
  for (id in names(flbrp_list)) {
    flbrp = flbrp_list[[id]]
    if (is.null(flbrp)) next
    
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
    
    results[[id]] = refpts_data
  }
  
  # Combine all results
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame(.id = character(), refpt = character(), value = numeric()))
  }
}

#' Process Breakpoint Data for Multiple Scenarios
#' 
#' Process breakpoint data for multiple scenarios and stock-recruitment relationships.
#' 
#' @param flbrp_list List of FLBRP objects
#' @param sr_params Data frame with stock-recruitment parameters
#' @param scenarios Character vector of scenario names
#' @param srrs Character vector of stock-recruitment relationship names
#' @return Processed breakpoint data
#' @export
processBreakpointData = function(flbrp_list, sr_params, scenarios, srrs) {
  results = list()
  
  for (scenario in scenarios) {
    for (srr in srrs) {
      scenario_name = paste0(scenario, "_", srr)
      
      if (scenario_name %in% names(flbrp_list)) {
        flbrp = flbrp_list[[scenario_name]]
        if (is.null(flbrp)) next
        
        # Extract breakpoint data
        breakpoint_data = data.frame(
          scenario = scenario,
          srr = srr,
          .id = names(flbrp_list)[names(flbrp_list) == scenario_name],
          stringsAsFactors = FALSE
        )
        
        # Add SRR parameters if available
        if (nrow(sr_params) > 0) {
          sr_match = sr_params$scenario == scenario & sr_params$srr == srr
          if (any(sr_match)) {
            breakpoint_data = cbind(breakpoint_data, sr_params[sr_match, ])
          }
        }
        
        results[[scenario_name]] = breakpoint_data
      }
    }
  }
  
  # Combine all results
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame(scenario = character(), srr = character(), .id = character()))
  }
}

#' Standardize Reference Point Data
#' 
#' Clean and standardize reference point data.
#' 
#' @param refpt_data Raw reference point data
#' @param remove_na Logical, whether to remove NA values (default: TRUE)
#' @return Cleaned and standardized data
#' @export
standardizeRefptData = function(refpt_data, remove_na = TRUE) {
  # Remove rows with all NA values
  if (remove_na) {
    refpt_data = refpt_data[!apply(refpt_data, 1, function(x) all(is.na(x))), ]
  }
  
  # Convert character columns to factors where appropriate
  char_cols = sapply(refpt_data, is.character)
  for (col in names(refpt_data)[char_cols]) {
    if (length(unique(refpt_data[[col]])) < nrow(refpt_data) * 0.5) {
      refpt_data[[col]] = as.factor(refpt_data[[col]])
    }
  }
  
  # Ensure numeric columns are properly formatted
  num_cols = sapply(refpt_data, is.numeric)
  for (col in names(refpt_data)[num_cols]) {
    refpt_data[[col]] = as.numeric(refpt_data[[col]])
  }
  
  return(refpt_data)
}

#' Combine Scenario Results
#' 
#' Combine results from multiple scenarios into a single data frame.
#' 
#' @param scenario_results List of scenario results
#' @param scenario_names Character vector of scenario names
#' @param srr_names Character vector of SRR names
#' @return Combined results data frame
#' @export
combineScenarioResults = function(scenario_results, scenario_names, srr_names) {
  results = list()
  
  for (scenario in scenario_names) {
    for (srr in srr_names) {
      scenario_name = paste0(scenario, "_", srr)
      
      if (scenario_name %in% names(scenario_results)) {
        result_data = scenario_results[[scenario_name]]
        if (!is.null(result_data)) {
          result_data$scenario = scenario
          result_data$srr = srr
          results[[scenario_name]] = result_data
        }
      }
    }
  }
  
  # Combine all results
  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame(scenario = character(), srr = character()))
  }
}

#' Merge Time Series and Reference Point Data
#' 
#' Merge time series data with reference point data for operating model analysis.
#' 
#' @param ts_data Time series data
#' @param refpt_data Reference point data
#' @param merge_by Character vector of column names to merge by
#' @return Merged operating model data
#' @export
mergeTimeSeriesRefpts = function(ts_data, refpt_data, merge_by = c(".id", "year")) {
  # Ensure merge columns exist in both datasets
  common_cols = intersect(names(ts_data), names(refpt_data))
  merge_cols = merge_by[merge_by %in% common_cols]
  
  if (length(merge_cols) == 0) {
    stop("No common columns found for merging")
  }
  
  # Merge datasets
  merged_data = merge(ts_data, refpt_data, by = merge_cols, all = TRUE)
  
  return(merged_data)
}

#' Clean Column Names
#' 
#' Clean up column names in merged data.
#' 
#' @param data Data frame with messy column names
#' @return Data frame with clean column names
#' @export
cleanColumnNames = function(data) {
  # Remove .x and .y suffixes from merged columns
  names(data) = gsub("\\.x$", "", names(data))
  names(data) = gsub("\\.y$", "", names(data))
  
  # Remove duplicate column names
  names(data) = make.unique(names(data))
  
  # Convert to lowercase and replace spaces with underscores
  names(data) = tolower(gsub(" ", "_", names(data)))
  
  return(data)
} 