#' Data Management Functions for FLRebuild Package
#' 
#' This file contains functions for organizing, cleaning, and managing data
#' in rebuild analysis workflows.

#' Create Scenario Matrix
#' 
#' Create matrix of different scenarios for analysis.
#' 
#' @param stock_recruitment_relationships Character vector of SRR names
#' @param prior_configurations List of prior configurations
#' @param scenario_names Character vector of scenario names (optional)
#' @return Scenario matrix data frame
#' @export
createScenarioMatrix = function(stock_recruitment_relationships, prior_configurations, scenario_names = NULL) {
  # Create base scenario combinations
  scenarios = expand.grid(
    srr = stock_recruitment_relationships,
    prior_config = names(prior_configurations),
    stringsAsFactors = FALSE
  )
  
  # Add scenario names if provided
  if (!is.null(scenario_names)) {
    scenarios$scenario_name = scenario_names[1:nrow(scenarios)]
  } else {
    scenarios$scenario_name = paste0(scenarios$srr, "_", scenarios$prior_config)
  }
  
  # Add prior configuration details
  for (i in 1:nrow(scenarios)) {
    config = prior_configurations[[scenarios$prior_config[i]]]
    for (param in names(config)) {
      scenarios[i, param] = config[[param]]
    }
  }
  
  return(scenarios)
}

#' Apply Scenario Analysis
#' 
#' Apply analysis across multiple scenarios.
#' 
#' @param base_data Base data for analysis
#' @param scenario_matrix Scenario matrix from createScenarioMatrix
#' @param analysis_function Function to apply to each scenario
#' @param analysis_args Additional arguments for analysis function
#' @return Results for all scenarios
#' @export
applyScenarioAnalysis = function(base_data, scenario_matrix, analysis_function, analysis_args = list()) {
  results = list()
  
  for (i in 1:nrow(scenario_matrix)) {
    scenario = scenario_matrix[i, ]
    scenario_name = scenario$scenario_name
    
    # Prepare scenario-specific arguments
    scenario_args = c(list(data = base_data), as.list(scenario), analysis_args)
    
    # Apply analysis function
    tryCatch({
      result = do.call(analysis_function, scenario_args)
      results[[scenario_name]] = result
    }, error = function(e) {
      warning(paste("Error in scenario", scenario_name, ":", e$message))
      results[[scenario_name]] = NULL
    })
  }
  
  return(results)
}

#' Save Analysis Results
#' 
#' Save analysis results with proper naming and metadata.
#' 
#' @param results Results object to save
#' @param file_path File path for saving
#' @param description Description of the results
#' @param format File format ("RData", "rds", "csv", "txt")
#' @param compress Logical, whether to compress the file
#' @return Saved file path
#' @export
saveAnalysisResults = function(results, file_path, description = "", format = "RData", compress = TRUE) {
  # Create directory if it doesn't exist
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  # Add metadata
  metadata = list(
    description = description,
    created = Sys.time(),
    session_info = sessionInfo(),
    format = format
  )
  
  # Save based on format
  if (format == "RData") {
    save(results, metadata, file = file_path, compress = compress)
  } else if (format == "rds") {
    saveRDS(list(results = results, metadata = metadata), file = file_path, compress = compress)
  } else if (format == "csv" && is.data.frame(results)) {
    write.csv(results, file = file_path, row.names = FALSE)
    # Save metadata separately
    writeLines(paste(names(metadata), metadata, sep = ": "), 
               paste0(file_path, ".meta"))
  } else if (format == "txt") {
    capture.output(results, file = file_path)
    # Save metadata separately
    writeLines(paste(names(metadata), metadata, sep = ": "), 
               paste0(file_path, ".meta"))
  } else {
    stop("Unsupported format: ", format)
  }
  
  return(file_path)
}

#' Export Results Summary
#' 
#' Export summary of analysis results in various formats.
#' 
#' @param results Results object
#' @param format Export format ("csv", "json", "html", "pdf")
#' @param file_path Output file path
#' @param summary_functions List of functions to apply for summarization
#' @return Exported file path
#' @export
exportResultsSummary = function(results, format = "csv", file_path = NULL, summary_functions = NULL) {
  # Generate summary if functions provided
  if (!is.null(summary_functions)) {
    summary_data = list()
    for (func_name in names(summary_functions)) {
      summary_data[[func_name]] = summary_functions[[func_name]](results)
    }
    results = summary_data
  }
  
  # Determine file path if not provided
  if (is.null(file_path)) {
    timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
    file_path = paste0("results_summary_", timestamp, ".", format)
  }
  
  # Export based on format
  if (format == "csv") {
    if (is.data.frame(results)) {
      write.csv(results, file = file_path, row.names = FALSE)
    } else if (is.list(results)) {
      # Combine list elements into data frame
      combined_data = do.call(rbind, lapply(names(results), function(name) {
        data.frame(metric = name, value = results[[name]], stringsAsFactors = FALSE)
      }))
      write.csv(combined_data, file = file_path, row.names = FALSE)
    }
  } else if (format == "json") {
    jsonlite::write_json(results, file_path, pretty = TRUE)
  } else if (format == "html") {
    # Create simple HTML table
    html_content = "<html><body><h1>Results Summary</h1>"
    if (is.data.frame(results)) {
      html_content = paste0(html_content, "<table border='1'>")
      html_content = paste0(html_content, "<tr><th>", paste(names(results), collapse = "</th><th>"), "</th></tr>")
      for (i in 1:nrow(results)) {
        html_content = paste0(html_content, "<tr><td>", paste(results[i, ], collapse = "</td><td>"), "</td></tr>")
      }
      html_content = paste0(html_content, "</table>")
    }
    html_content = paste0(html_content, "</body></html>")
    writeLines(html_content, file_path)
  } else if (format == "pdf") {
    # Create PDF report (requires rmarkdown)
    if (requireNamespace("rmarkdown", quietly = TRUE)) {
      temp_rmd = tempfile(fileext = ".Rmd")
      writeLines(c(
        "---",
        "title: Results Summary",
        "output: pdf_document",
        "---",
        "",
        "```{r results, echo=FALSE}",
        "knitr::kable(results)",
        "```"
      ), temp_rmd)
      rmarkdown::render(temp_rmd, output_file = file_path)
      unlink(temp_rmd)
    } else {
      stop("rmarkdown package required for PDF export")
    }
  } else {
    stop("Unsupported format: ", format)
  }
  
  return(file_path)
}

#' Validate FLBRP Object Structure
#' 
#' Validate FLBRP object structure and contents.
#' 
#' @param flbrp FLBRP object to validate
#' @param required_components Character vector of required components
#' @return Validation status
#' @export
validateFLBRPObject = function(flbrp, required_components = c("refpts", "params")) {
  validation = list(
    is_valid = TRUE,
    errors = character(),
    warnings = character()
  )
  
  # Check if object is FLBRP
  if (!inherits(flbrp, "FLBRP")) {
    validation$is_valid = FALSE
    validation$errors = c(validation$errors, "Object is not an FLBRP")
    return(validation)
  }
  
  # Check required components
  for (component in required_components) {
    if (!exists(component, where = flbrp)) {
      validation$is_valid = FALSE
      validation$errors = c(validation$errors, paste("Missing component:", component))
    }
  }
  
  # Check reference points
  if ("refpts" %in% required_components && exists("refpts", where = flbrp)) {
    refpts_obj = refpts(flbrp)
    if (is.null(refpts_obj) || dim(refpts_obj)[1] == 0) {
      validation$warnings = c(validation$warnings, "Reference points are empty or NULL")
    }
  }
  
  # Check parameters
  if ("params" %in% required_components && exists("params", where = flbrp)) {
    params_obj = params(flbrp)
    if (is.null(params_obj) || dim(params_obj)[1] == 0) {
      validation$warnings = c(validation$warnings, "Parameters are empty or NULL")
    }
  }
  
  return(validation)
}

#' Clean and Prepare Data for Analysis
#' 
#' Clean and prepare data for rebuild analysis.
#' 
#' @param data Input data
#' @param remove_outliers Logical, whether to remove outliers
#' @param outlier_method Method for outlier detection ("iqr", "zscore", "mad")
#' @param outlier_threshold Threshold for outlier detection
#' @param fill_missing Method for filling missing values ("mean", "median", "interpolate")
#' @return Cleaned data
#' @export
cleanAndPrepareData = function(data, remove_outliers = TRUE, outlier_method = "iqr", 
                              outlier_threshold = 1.5, fill_missing = "interpolate") {
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame")
  }
  
  # Identify numeric columns
  numeric_cols = sapply(data, is.numeric)
  
  # Fill missing values
  if (fill_missing != "none") {
    for (col in names(data)[numeric_cols]) {
      if (any(is.na(data[[col]]))) {
        if (fill_missing == "mean") {
          data[[col]][is.na(data[[col]])] = mean(data[[col]], na.rm = TRUE)
        } else if (fill_missing == "median") {
          data[[col]][is.na(data[[col]])] = median(data[[col]], na.rm = TRUE)
        } else if (fill_missing == "interpolate") {
          data[[col]] = zoo::na.approx(data[[col]], na.rm = FALSE)
        }
      }
    }
  }
  
  # Remove outliers
  if (remove_outliers) {
    for (col in names(data)[numeric_cols]) {
      outliers = detectOutliers(data[[col]], method = outlier_method, threshold = outlier_threshold)
      if (length(outliers) > 0) {
        data = data[-outliers, ]
      }
    }
  }
  
  return(data)
}

#' Detect Outliers in Data
#' 
#' Detect outliers using various methods.
#' 
#' @param x Numeric vector
#' @param method Method for outlier detection ("iqr", "zscore", "mad")
#' @param threshold Threshold for outlier detection
#' @return Indices of outliers
#' @export
detectOutliers = function(x, method = "iqr", threshold = 1.5) {
  if (method == "iqr") {
    q1 = quantile(x, 0.25, na.rm = TRUE)
    q3 = quantile(x, 0.75, na.rm = TRUE)
    iqr = q3 - q1
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr
    outliers = which(x < lower_bound | x > upper_bound)
  } else if (method == "zscore") {
    z_scores = abs((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    outliers = which(z_scores > threshold)
  } else if (method == "mad") {
    median_val = median(x, na.rm = TRUE)
    mad_val = median(abs(x - median_val), na.rm = TRUE)
    outliers = which(abs(x - median_val) > threshold * mad_val)
  } else {
    stop("Unsupported outlier detection method: ", method)
  }
  
  return(outliers)
} 