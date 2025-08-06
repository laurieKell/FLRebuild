#' Scenario Management Functions for FLRebuild Package
#' 
#' This file contains functions for managing and organizing different scenarios
#' in rebuild analysis workflows.

#' Create Scenario Configuration
#' 
#' Create configuration for different analysis scenarios.
#' 
#' @param scenario_name Name of the scenario
#' @param stock_recruitment_relationship SRR model to use
#' @param prior_configuration Prior configuration parameters
#' @param bootstrap_settings Bootstrap analysis settings
#' @param rebuild_settings Rebuild analysis settings
#' @return Scenario configuration object
#' @export
createScenarioConfig = function(scenario_name, stock_recruitment_relationship, 
                               prior_configuration = NULL, bootstrap_settings = NULL, 
                               rebuild_settings = NULL) {
  config = list(
    name = scenario_name,
    srr = stock_recruitment_relationship,
    priors = prior_configuration,
    bootstrap = bootstrap_settings,
    rebuild = rebuild_settings,
    created = Sys.time()
  )
  
  # Set default bootstrap settings
  if (is.null(bootstrap_settings)) {
    config$bootstrap = list(
      nits = 100,
      ncores = NULL,
      ftmb_args = list(
        model = "bevholtSV",
        s.est = TRUE,
        s = 0.7,
        s.logitsd = 0.4,
        spr0 = 0.7
      )
    )
  }
  
  # Set default rebuild settings
  if (is.null(rebuild_settings)) {
    config$rebuild = list(
      target = 0.3,
      max_year = 50,
      method = "fbar"
    )
  }
  
  return(config)
}

#' Apply Scenario Configuration
#' 
#' Apply scenario configuration to FLBRP objects.
#' 
#' @param flbrp_list List of FLBRP objects
#' @param scenario_config Scenario configuration object
#' @param analysis_function Function to apply with scenario configuration
#' @return Results of applying scenario configuration
#' @export
applyScenarioConfig = function(flbrp_list, scenario_config, analysis_function) {
  # Validate scenario configuration
  if (!is.list(scenario_config) || is.null(scenario_config$name)) {
    stop("Invalid scenario configuration")
  }
  
  # Apply configuration to each FLBRP object
  results = list()
  
  for (id in names(flbrp_list)) {
    flbrp = flbrp_list[[id]]
    if (is.null(flbrp)) next
    
    tryCatch({
      # Apply scenario-specific settings
      if (!is.null(scenario_config$priors)) {
        # Apply prior configuration
        flbrp = applyPriors(flbrp, scenario_config$priors)
      }
      
      # Apply analysis function with scenario configuration
      result = analysis_function(flbrp, scenario_config)
      results[[id]] = result
      
    }, error = function(e) {
      warning(paste("Error applying scenario", scenario_config$name, "to", id, ":", e$message))
      results[[id]] = NULL
    })
  }
  
  return(results)
}

#' Apply Priors to FLBRP Object
#' 
#' Apply prior configuration to FLBRP object.
#' 
#' @param flbrp FLBRP object
#' @param prior_config Prior configuration
#' @return FLBRP object with applied priors
#' @export
applyPriors = function(flbrp, prior_config) {
  if (is.null(prior_config)) {
    return(flbrp)
  }
  
  # Apply steepness prior if specified
  if (!is.null(prior_config$s) && !is.null(prior_config$cv_s)) {
    params(flbrp)["s", ] = prior_config$s
    # Note: cv_s would be used in the fitting process, not directly in params
  }
  
  # Apply r0 prior if specified
  if (!is.null(prior_config$r0) && !is.null(prior_config$cv_r0)) {
    params(flbrp)["r0", ] = prior_config$r0
    # Note: cv_r0 would be used in the fitting process, not directly in params
  }
  
  return(flbrp)
}

#' Compare Scenarios
#' 
#' Compare results across different scenarios.
#' 
#' @param scenario_results List of results from different scenarios
#' @param comparison_metric Metric to compare
#' @param plot_type Type of comparison plot ("boxplot", "violin", "bar")
#' @return Comparison results and plot
#' @export
compareScenarios = function(scenario_results, comparison_metric, plot_type = "boxplot") {
  # Extract comparison data
  comparison_data = data.frame()
  
  for (scenario_name in names(scenario_results)) {
    scenario_data = scenario_results[[scenario_name]]
    if (!is.null(scenario_data) && is.data.frame(scenario_data)) {
      if (comparison_metric %in% names(scenario_data)) {
        temp_data = data.frame(
          scenario = scenario_name,
          value = scenario_data[[comparison_metric]],
          stringsAsFactors = FALSE
        )
        comparison_data = rbind(comparison_data, temp_data)
      }
    }
  }
  
  if (nrow(comparison_data) == 0) {
    stop("No comparison data found for metric: ", comparison_metric)
  }
  
  # Create comparison plot
  if (plot_type == "boxplot") {
    p = ggplot(comparison_data, aes(x = scenario, y = value)) +
      geom_boxplot(fill = "steelblue", alpha = 0.7) +
      labs(title = paste("Comparison of", comparison_metric, "Across Scenarios"),
           x = "Scenario", y = comparison_metric)
  } else if (plot_type == "violin") {
    p = ggplot(comparison_data, aes(x = scenario, y = value)) +
      geom_violin(fill = "steelblue", alpha = 0.7) +
      labs(title = paste("Comparison of", comparison_metric, "Across Scenarios"),
           x = "Scenario", y = comparison_metric)
  } else if (plot_type == "bar") {
    summary_data = comparison_data %>%
      group_by(scenario) %>%
      summarise(mean_value = mean(value, na.rm = TRUE),
                se_value = sd(value, na.rm = TRUE) / sqrt(n()))
    
    p = ggplot(summary_data, aes(x = scenario, y = mean_value)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                    width = 0.2) +
      labs(title = paste("Comparison of", comparison_metric, "Across Scenarios"),
           x = "Scenario", y = paste("Mean", comparison_metric))
  } else {
    stop("Unsupported plot type: ", plot_type)
  }
  
  # Apply theme
  p = p + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(data = comparison_data, plot = p))
}

#' Generate Scenario Report
#' 
#' Generate a comprehensive report for scenario analysis.
#' 
#' @param scenario_results List of scenario results
#' @param scenario_configs List of scenario configurations
#' @param report_type Type of report ("summary", "detailed", "comparative")
#' @return Report object
#' @export
generateScenarioReport = function(scenario_results, scenario_configs, report_type = "summary") {
  report = list(
    timestamp = Sys.time(),
    type = report_type,
    scenarios = names(scenario_results),
    summary = list()
  )
  
  # Generate summary statistics for each scenario
  for (scenario_name in names(scenario_results)) {
    scenario_data = scenario_results[[scenario_name]]
    scenario_config = scenario_configs[[scenario_name]]
    
    if (!is.null(scenario_data)) {
      # Calculate summary statistics
      if (is.data.frame(scenario_data)) {
        numeric_cols = sapply(scenario_data, is.numeric)
        if (any(numeric_cols)) {
          summary_stats = lapply(scenario_data[, numeric_cols, drop = FALSE], function(x) {
            c(mean = mean(x, na.rm = TRUE),
              median = median(x, na.rm = TRUE),
              sd = sd(x, na.rm = TRUE),
              min = min(x, na.rm = TRUE),
              max = max(x, na.rm = TRUE),
              n = length(x),
              n_missing = sum(is.na(x)))
          })
          
          report$summary[[scenario_name]] = list(
            config = scenario_config,
            statistics = summary_stats
          )
        }
      }
    }
  }
  
  # Add comparative analysis if multiple scenarios
  if (length(scenario_results) > 1 && report_type %in% c("detailed", "comparative")) {
    report$comparison = compareMultipleScenarios(scenario_results)
  }
  
  return(report)
}

#' Compare Multiple Scenarios
#' 
#' Compare multiple scenarios using statistical tests.
#' 
#' @param scenario_results List of scenario results
#' @return Comparison results
#' @export
compareMultipleScenarios = function(scenario_results) {
  comparison = list()
  
  # Extract common metrics across scenarios
  common_metrics = NULL
  for (scenario_name in names(scenario_results)) {
    scenario_data = scenario_results[[scenario_name]]
    if (is.data.frame(scenario_data)) {
      numeric_cols = names(scenario_data)[sapply(scenario_data, is.numeric)]
      if (is.null(common_metrics)) {
        common_metrics = numeric_cols
      } else {
        common_metrics = intersect(common_metrics, numeric_cols)
      }
    }
  }
  
  if (length(common_metrics) == 0) {
    return(NULL)
  }
  
  # Perform statistical comparisons for each metric
  for (metric in common_metrics) {
    metric_data = data.frame()
    
    for (scenario_name in names(scenario_results)) {
      scenario_data = scenario_results[[scenario_name]]
      if (is.data.frame(scenario_data) && metric %in% names(scenario_data)) {
        temp_data = data.frame(
          scenario = scenario_name,
          value = scenario_data[[metric]],
          stringsAsFactors = FALSE
        )
        metric_data = rbind(metric_data, temp_data)
      }
    }
    
    if (nrow(metric_data) > 0) {
      # Perform ANOVA if multiple scenarios
      if (length(unique(metric_data$scenario)) > 1) {
        tryCatch({
          anova_result = aov(value ~ scenario, data = metric_data)
          comparison[[metric]] = list(
            anova = summary(anova_result),
            means = tapply(metric_data$value, metric_data$scenario, mean, na.rm = TRUE),
            sds = tapply(metric_data$value, metric_data$scenario, sd, na.rm = TRUE)
          )
        }, error = function(e) {
          comparison[[metric]] = list(error = e$message)
        })
      }
    }
  }
  
  return(comparison)
}

#' Validate Scenario Configuration
#' 
#' Validate scenario configuration for completeness and correctness.
#' 
#' @param scenario_config Scenario configuration object
#' @return Validation results
#' @export
validateScenarioConfig = function(scenario_config) {
  validation = list(
    is_valid = TRUE,
    errors = character(),
    warnings = character()
  )
  
  # Check required fields
  required_fields = c("name", "srr")
  for (field in required_fields) {
    if (is.null(scenario_config[[field]])) {
      validation$is_valid = FALSE
      validation$errors = c(validation$errors, paste("Missing required field:", field))
    }
  }
  
  # Validate SRR model
  valid_srr_models = c("bevholt", "bevholtSV", "ricker", "segreg", "shepherd")
  if (!is.null(scenario_config$srr) && !scenario_config$srr %in% valid_srr_models) {
    validation$warnings = c(validation$warnings, 
                           paste("Unknown SRR model:", scenario_config$srr))
  }
  
  # Validate prior configuration
  if (!is.null(scenario_config$priors)) {
    prior_validation = validatePriorConfig(scenario_config$priors)
    validation$errors = c(validation$errors, prior_validation$errors)
    validation$warnings = c(validation$warnings, prior_validation$warnings)
  }
  
  # Validate bootstrap settings
  if (!is.null(scenario_config$bootstrap)) {
    bootstrap_validation = validateBootstrapConfig(scenario_config$bootstrap)
    validation$errors = c(validation$errors, bootstrap_validation$errors)
    validation$warnings = c(validation$warnings, bootstrap_validation$warnings)
  }
  
  return(validation)
}

#' Validate Prior Configuration
#' 
#' Validate prior configuration parameters.
#' 
#' @param prior_config Prior configuration
#' @return Validation results
#' @export
validatePriorConfig = function(prior_config) {
  validation = list(
    is_valid = TRUE,
    errors = character(),
    warnings = character()
  )
  
  # Check steepness prior
  if (!is.null(prior_config$s)) {
    if (prior_config$s <= 0 || prior_config$s >= 1) {
      validation$warnings = c(validation$warnings, "Steepness prior should be between 0 and 1")
    }
    if (is.null(prior_config$cv_s)) {
      validation$warnings = c(validation$warnings, "CV for steepness prior not specified")
    }
  }
  
  # Check r0 prior
  if (!is.null(prior_config$r0)) {
    if (prior_config$r0 <= 0) {
      validation$warnings = c(validation$warnings, "r0 prior should be positive")
    }
    if (is.null(prior_config$cv_r0)) {
      validation$warnings = c(validation$warnings, "CV for r0 prior not specified")
    }
  }
  
  return(validation)
}

#' Validate Bootstrap Configuration
#' 
#' Validate bootstrap configuration parameters.
#' 
#' @param bootstrap_config Bootstrap configuration
#' @return Validation results
#' @export
validateBootstrapConfig = function(bootstrap_config) {
  validation = list(
    is_valid = TRUE,
    errors = character(),
    warnings = character()
  )
  
  # Check number of iterations
  if (!is.null(bootstrap_config$nits)) {
    if (bootstrap_config$nits <= 0) {
      validation$errors = c(validation$errors, "Number of bootstrap iterations must be positive")
    }
    if (bootstrap_config$nits < 100) {
      validation$warnings = c(validation$warnings, "Consider using at least 100 bootstrap iterations")
    }
  }
  
  # Check number of cores
  if (!is.null(bootstrap_config$ncores)) {
    if (bootstrap_config$ncores <= 0) {
      validation$errors = c(validation$errors, "Number of cores must be positive")
    }
    if (bootstrap_config$ncores > detectCores()) {
      validation$warnings = c(validation$warnings, "Number of cores exceeds available cores")
    }
  }
  
  return(validation)
} 