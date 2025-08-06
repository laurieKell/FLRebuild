#' Analysis and Visualization Functions for FLRebuild Package
#' 
#' This file contains functions for creating analysis plots and visualizations
#' for rebuild analysis workflows.

#' Plot Reference Point Distribution
#' 
#' Create reference point distribution plots.
#' 
#' @param refpt_data Reference point data
#' @param refpt_type Type of reference point to plot
#' @param col_names Column names for plotting
#' @param plot_type Type of plot ("histogram", "density", "boxplot", "violin")
#' @param facet_by Variable to facet by (optional)
#' @return ggplot object
#' @export
plotRefptDistribution = function(refpt_data, refpt_type, col_names, plot_type = "histogram", facet_by = NULL) {
  # Filter data for specific reference point type
  plot_data = refpt_data[refpt_data$refpt == refpt_type, ]
  
  if (nrow(plot_data) == 0) {
    stop("No data found for reference point type: ", refpt_type)
  }
  
  # Create base plot
  if (plot_type == "histogram") {
    p = ggplot(plot_data, aes_string(x = col_names$value)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      labs(title = paste("Distribution of", refpt_type, "Reference Points"),
           x = "Value", y = "Frequency")
  } else if (plot_type == "density") {
    p = ggplot(plot_data, aes_string(x = col_names$value)) +
      geom_density(fill = "steelblue", alpha = 0.7) +
      labs(title = paste("Density of", refpt_type, "Reference Points"),
           x = "Value", y = "Density")
  } else if (plot_type == "boxplot") {
    p = ggplot(plot_data, aes_string(y = col_names$value)) +
      geom_boxplot(fill = "steelblue", alpha = 0.7) +
      labs(title = paste("Boxplot of", refpt_type, "Reference Points"),
           y = "Value")
  } else if (plot_type == "violin") {
    p = ggplot(plot_data, aes_string(y = col_names$value)) +
      geom_violin(fill = "steelblue", alpha = 0.7) +
      labs(title = paste("Violin Plot of", refpt_type, "Reference Points"),
           y = "Value")
  } else {
    stop("Unsupported plot type: ", plot_type)
  }
  
  # Add faceting if specified
  if (!is.null(facet_by) && facet_by %in% names(plot_data)) {
    p = p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  # Apply theme
  p = p + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  return(p)
}

#' Plot Rebuild Times
#' 
#' Create rebuild time analysis plots.
#' 
#' @param rebuild_data Rebuild time data
#' @param col_names Column names for plotting
#' @param plot_type Type of plot ("histogram", "density", "boxplot", "scatter")
#' @param x_var Variable for x-axis (for scatter plots)
#' @param facet_by Variable to facet by (optional)
#' @return ggplot object
#' @export
plotRebuildTimes = function(rebuild_data, col_names, plot_type = "histogram", x_var = NULL, facet_by = NULL) {
  if (nrow(rebuild_data) == 0) {
    stop("No rebuild time data provided")
  }
  
  # Create base plot
  if (plot_type == "histogram") {
    p = ggplot(rebuild_data, aes_string(x = col_names$rebuild_time)) +
      geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
      labs(title = "Distribution of Rebuild Times",
           x = "Rebuild Time (years)", y = "Frequency")
  } else if (plot_type == "density") {
    p = ggplot(rebuild_data, aes_string(x = col_names$rebuild_time)) +
      geom_density(fill = "darkgreen", alpha = 0.7) +
      labs(title = "Density of Rebuild Times",
           x = "Rebuild Time (years)", y = "Density")
  } else if (plot_type == "boxplot") {
    p = ggplot(rebuild_data, aes_string(y = col_names$rebuild_time)) +
      geom_boxplot(fill = "darkgreen", alpha = 0.7) +
      labs(title = "Boxplot of Rebuild Times",
           y = "Rebuild Time (years)")
  } else if (plot_type == "scatter") {
    if (is.null(x_var) || !x_var %in% names(rebuild_data)) {
      stop("x_var must be specified and present in data for scatter plots")
    }
    p = ggplot(rebuild_data, aes_string(x = x_var, y = col_names$rebuild_time)) +
      geom_point(alpha = 0.7, color = "darkgreen") +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(title = "Rebuild Time vs. Predictor Variable",
           x = x_var, y = "Rebuild Time (years)")
  } else {
    stop("Unsupported plot type: ", plot_type)
  }
  
  # Add faceting if specified
  if (!is.null(facet_by) && facet_by %in% names(rebuild_data)) {
    p = p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  # Apply theme
  p = p + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  return(p)
}

#' Plot Operating Model Performance
#' 
#' Create operating model performance plots.
#' 
#' @param performance_data Performance data
#' @param refpts Reference points to evaluate
#' @param plot_type Type of plot ("scatter", "heatmap", "bar")
#' @param metric Performance metric to plot
#' @param facet_by Variable to facet by (optional)
#' @return ggplot object
#' @export
plotOMPerformance = function(performance_data, refpts, plot_type = "scatter", metric = "bias", facet_by = NULL) {
  if (nrow(performance_data) == 0) {
    stop("No performance data provided")
  }
  
  # Filter data for specified reference points
  if (!is.null(refpts)) {
    performance_data = performance_data[performance_data$refpt %in% refpts, ]
  }
  
  if (nrow(performance_data) == 0) {
    stop("No data found for specified reference points")
  }
  
  # Create base plot
  if (plot_type == "scatter") {
    p = ggplot(performance_data, aes_string(x = "true_value", y = "estimated_value")) +
      geom_point(alpha = 0.7, color = "blue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(title = "Operating Model Performance",
           x = "True Value", y = "Estimated Value")
  } else if (plot_type == "heatmap") {
    # Create correlation matrix for heatmap
    cor_matrix = performance_data %>%
      group_by(refpt) %>%
      summarise(correlation = cor(true_value, estimated_value, use = "complete.obs")) %>%
      spread(refpt, correlation)
    
    p = ggplot(cor_matrix, aes(x = refpt, y = 1, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
      labs(title = "Performance Correlation Heatmap",
           x = "Reference Point", y = "", fill = "Correlation")
  } else if (plot_type == "bar") {
    p = ggplot(performance_data, aes_string(x = "refpt", y = metric)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      labs(title = paste("Performance by Reference Point -", metric),
           x = "Reference Point", y = metric)
  } else {
    stop("Unsupported plot type: ", plot_type)
  }
  
  # Add faceting if specified
  if (!is.null(facet_by) && facet_by %in% names(performance_data)) {
    p = p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  # Apply theme
  p = p + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  return(p)
}

#' Plot Reference Point Correlations
#' 
#' Create correlation plots between reference points.
#' 
#' @param refpt_data Reference point data
#' @param refpts Reference points to include in correlation analysis
#' @param plot_type Type of plot ("correlation", "scatter_matrix", "heatmap")
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @return ggplot object or list of plots
#' @export
plotRefptCorrelations = function(refpt_data, refpts, plot_type = "correlation", method = "pearson") {
  # Filter data for specified reference points
  if (!is.null(refpts)) {
    refpt_data = refpt_data[refpt_data$refpt %in% refpts, ]
  }
  
  if (nrow(refpt_data) == 0) {
    stop("No data found for specified reference points")
  }
  
  # Reshape data for correlation analysis
  wide_data = refpt_data %>%
    select(.id, refpt, value) %>%
    spread(refpt, value)
  
  # Remove .id column for correlation
  cor_data = wide_data[, -which(names(wide_data) == ".id")]
  
  if (plot_type == "correlation") {
    # Create correlation matrix
    cor_matrix = cor(cor_data, use = "complete.obs", method = method)
    
    # Convert to long format for plotting
    cor_long = as.data.frame(cor_matrix) %>%
      rownames_to_column("var1") %>%
      gather(var2, correlation, -var1)
    
    p = ggplot(cor_long, aes(x = var1, y = var2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
      labs(title = paste("Reference Point Correlations (", method, ")"),
           x = "", y = "", fill = "Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  } else if (plot_type == "scatter_matrix") {
    # Create scatter plot matrix
    p = ggpairs(cor_data, 
                title = paste("Reference Point Scatter Matrix (", method, " correlation)"),
                lower = list(continuous = "smooth"))
    
    return(p)
  } else if (plot_type == "heatmap") {
    # Create heatmap with hierarchical clustering
    cor_matrix = cor(cor_data, use = "complete.obs", method = method)
    
    # Convert to long format
    cor_long = as.data.frame(cor_matrix) %>%
      rownames_to_column("var1") %>%
      gather(var2, correlation, -var1)
    
    p = ggplot(cor_long, aes(x = var1, y = var2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
      labs(title = paste("Reference Point Correlation Heatmap (", method, ")"),
           x = "", y = "", fill = "Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  } else {
    stop("Unsupported plot type: ", plot_type)
  }
}

#' Create Summary Statistics Plot
#' 
#' Create summary statistics visualization.
#' 
#' @param data Data for summary statistics
#' @param group_var Variable to group by
#' @param value_var Variable to summarize
#' @param stat_type Type of statistic ("mean", "median", "sd", "range")
#' @param plot_type Type of plot ("bar", "line", "point")
#' @return ggplot object
#' @export
createSummaryStatsPlot = function(data, group_var, value_var, stat_type = "mean", plot_type = "bar") {
  # Calculate summary statistics
  summary_data = data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      mean_val = mean(!!sym(value_var), na.rm = TRUE),
      median_val = median(!!sym(value_var), na.rm = TRUE),
      sd_val = sd(!!sym(value_var), na.rm = TRUE),
      min_val = min(!!sym(value_var), na.rm = TRUE),
      max_val = max(!!sym(value_var), na.rm = TRUE)
    )
  
  # Select appropriate statistic
  if (stat_type == "mean") {
    summary_data$stat_value = summary_data$mean_val
    stat_label = "Mean"
  } else if (stat_type == "median") {
    summary_data$stat_value = summary_data$median_val
    stat_label = "Median"
  } else if (stat_type == "sd") {
    summary_data$stat_value = summary_data$sd_val
    stat_label = "Standard Deviation"
  } else if (stat_type == "range") {
    summary_data$stat_value = summary_data$max_val - summary_data$min_val
    stat_label = "Range"
  } else {
    stop("Unsupported statistic type: ", stat_type)
  }
  
  # Create plot
  if (plot_type == "bar") {
    p = ggplot(summary_data, aes_string(x = group_var, y = "stat_value")) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      labs(title = paste(stat_label, "by", group_var),
           x = group_var, y = stat_label)
  } else if (plot_type == "line") {
    p = ggplot(summary_data, aes_string(x = group_var, y = "stat_value", group = 1)) +
      geom_line(color = "steelblue", size = 1) +
      geom_point(color = "steelblue", size = 3) +
      labs(title = paste(stat_label, "by", group_var),
           x = group_var, y = stat_label)
  } else if (plot_type == "point") {
    p = ggplot(summary_data, aes_string(x = group_var, y = "stat_value")) +
      geom_point(color = "steelblue", size = 3) +
      labs(title = paste(stat_label, "by", group_var),
           x = group_var, y = stat_label)
  } else {
    stop("Unsupported plot type: ", plot_type)
  }
  
  # Apply theme
  p = p + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  return(p)
} 