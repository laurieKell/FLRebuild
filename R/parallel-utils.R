#' Parallel Processing Utility Functions for FLRebuild Package
#' 
#' This file contains utility functions for setting up and managing parallel processing
#' in rebuild analysis workflows.

#' Setup Parallel Processing Environment
#' 
#' Setup parallel processing environment for analysis functions.
#' 
#' @param ncores Number of cores to use (default=detectCores()-4)
#' @param packages Character vector of packages to load in parallel environment
#' @return Parallel processing cluster
#' @export
setupParallelProcessing = function(ncores = NULL, packages = c("FLCore", "FLBRP", "FLRebuild")) {
  if (is.null(ncores)) {
    ncores = detectCores() - 4
  }
  
  # Ensure ncores doesn't exceed available cores
  ncores = min(ncores, detectCores())
  
  cl = makeCluster(ncores)
  registerDoParallel(cl)
  
  # Load required packages in parallel environment
  clusterEvalQ(cl, {
    for (pkg in packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        warning(paste("Package", pkg, "not available in parallel environment"))
      }
    }
  })
  
  return(cl)
}

#' Cleanup Parallel Processing Resources
#' 
#' Clean up parallel processing resources.
#' 
#' @param cl Parallel processing cluster
#' @export
cleanupParallelProcessing = function(cl) {
  if (!is.null(cl)) {
    tryCatch({
      stopCluster(cl)
    }, error = function(e) {
      warning("Error stopping cluster:", e$message)
    })
  }
}

#' Validate Input Data for Analysis Functions
#' 
#' Validate input data for analysis functions.
#' 
#' @param data Input data object
#' @param required_cols Character vector of required column names
#' @param data_type Character string describing the type of data
#' @return Validation results
#' @export
validateInputData = function(data, required_cols = NULL, data_type = "data") {
  validation_results = list(
    is_valid = TRUE,
    errors = character(),
    warnings = character()
  )
  
  # Check if data is NULL
  if (is.null(data)) {
    validation_results$is_valid = FALSE
    validation_results$errors = c(validation_results$errors, 
                                 paste(data_type, "is NULL"))
    return(validation_results)
  }
  
  # Check if data is empty
  if (is.data.frame(data) && nrow(data) == 0) {
    validation_results$warnings = c(validation_results$warnings, 
                                   paste(data_type, "is empty"))
  }
  
  # Check required columns if specified
  if (!is.null(required_cols) && is.data.frame(data)) {
    missing_cols = setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      validation_results$is_valid = FALSE
      validation_results$errors = c(validation_results$errors, 
                                   paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Check for NA values in critical columns
  if (is.data.frame(data) && !is.null(required_cols)) {
    for (col in required_cols) {
      if (col %in% names(data)) {
        na_count = sum(is.na(data[[col]]))
        if (na_count > 0) {
          validation_results$warnings = c(validation_results$warnings, 
                                         paste(na_count, "NA values found in column", col))
        }
      }
    }
  }
  
  return(validation_results)
}

#' Safe Execution of Functions with Error Handling
#' 
#' Execute functions safely with comprehensive error handling.
#' 
#' @param expr Expression to evaluate
#' @param error_value Value to return if error occurs (default: NULL)
#' @param warning_msg Custom warning message
#' @return Function result or error_value if error occurs
#' @export
safeExecute = function(expr, error_value = NULL, warning_msg = NULL) {
  tryCatch({
    eval(expr)
  }, error = function(e) {
    if (!is.null(warning_msg)) {
      warning(warning_msg, ": ", e$message)
    } else {
      warning("Error in safe execution: ", e$message)
    }
    return(error_value)
  }, warning = function(w) {
    warning("Warning in safe execution: ", w$message)
    return(eval(expr))
  })
}

#' Check System Resources for Parallel Processing
#' 
#' Check available system resources for parallel processing.
#' 
#' @return System resource information
#' @export
checkSystemResources = function() {
  resources = list(
    total_cores = detectCores(),
    available_cores = detectCores(logical = FALSE),
    memory_gb = as.numeric(system("wmic computersystem get TotalPhysicalMemory", intern = TRUE)[2]) / 1024^3,
    os = Sys.info()["sysname"]
  )
  
  # Estimate recommended cores for parallel processing
  resources$recommended_cores = max(1, min(resources$total_cores - 2, 
                                          floor(resources$total_cores * 0.75)))
  
  return(resources)
}

#' Optimize Parallel Processing Settings
#' 
#' Optimize parallel processing settings based on system resources.
#' 
#' @param task_complexity Character string indicating task complexity ("low", "medium", "high")
#' @param memory_requirement_gb Estimated memory requirement per task in GB
#' @return Optimized parallel processing settings
#' @export
optimizeParallelSettings = function(task_complexity = "medium", memory_requirement_gb = 1) {
  resources = checkSystemResources()
  
  # Base settings
  settings = list(
    ncores = resources$recommended_cores,
    memory_per_core_gb = memory_requirement_gb,
    chunk_size = 1
  )
  
  # Adjust based on task complexity
  if (task_complexity == "low") {
    settings$ncores = min(settings$ncores, 4)
    settings$chunk_size = 10
  } else if (task_complexity == "high") {
    settings$ncores = max(1, settings$ncores - 2)
    settings$chunk_size = 1
  }
  
  # Adjust based on memory requirements
  available_memory_per_core = resources$memory_gb / settings$ncores
  if (memory_requirement_gb > available_memory_per_core * 0.8) {
    settings$ncores = max(1, floor(resources$memory_gb / (memory_requirement_gb * 1.25)))
  }
  
  return(settings)
}

#' Monitor Parallel Processing Progress
#' 
#' Monitor progress of parallel processing tasks.
#' 
#' @param total_tasks Total number of tasks
#' @param completed_tasks Number of completed tasks
#' @param start_time Start time of processing
#' @return Progress information
#' @export
monitorProgress = function(total_tasks, completed_tasks, start_time) {
  progress = list(
    total_tasks = total_tasks,
    completed_tasks = completed_tasks,
    remaining_tasks = total_tasks - completed_tasks,
    completion_rate = completed_tasks / total_tasks * 100,
    elapsed_time = Sys.time() - start_time
  )
  
  # Estimate remaining time
  if (completed_tasks > 0) {
    avg_time_per_task = as.numeric(progress$elapsed_time) / completed_tasks
    progress$estimated_remaining_time = avg_time_per_task * progress$remaining_tasks
  } else {
    progress$estimated_remaining_time = NA
  }
  
  return(progress)
} 