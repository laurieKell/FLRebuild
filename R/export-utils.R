#' Data Export Utility Functions for FLRebuild Package
#' 
#' This file contains utility functions for exporting data and results
#' from rebuild analysis workflows.

#' Export Results to Multiple Formats
#' 
#' Export analysis results to multiple file formats.
#' 
#' @param results Results object to export
#' @param base_filename Base filename without extension
#' @param output_dir Output directory
#' @param formats Character vector of formats to export ("csv", "rds", "RData", "json", "html")
#' @param include_metadata Logical, whether to include metadata
#' @return List of exported file paths
#' @export
exportResultsToMultipleFormats = function(results, base_filename, output_dir = ".", 
                                        formats = c("csv", "rds"), include_metadata = TRUE) {
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  exported_files = list()
  
  for (format in formats) {
    filename = file.path(output_dir, paste0(base_filename, ".", format))
    
    tryCatch({
      if (format == "csv") {
        if (is.data.frame(results)) {
          write.csv(results, file = filename, row.names = FALSE)
        } else if (is.list(results)) {
          # Export each list element as separate CSV
          for (name in names(results)) {
            if (is.data.frame(results[[name]])) {
              csv_filename = file.path(output_dir, paste0(base_filename, "_", name, ".csv"))
              write.csv(results[[name]], file = csv_filename, row.names = FALSE)
              exported_files[[paste0(name, "_csv")]] = csv_filename
            }
          }
        }
        exported_files[["csv"]] = filename
      } else if (format == "rds") {
        saveRDS(results, file = filename)
        exported_files[["rds"]] = filename
      } else if (format == "RData") {
        save(results, file = filename)
        exported_files[["RData"]] = filename
      } else if (format == "json") {
        jsonlite::write_json(results, filename, pretty = TRUE)
        exported_files[["json"]] = filename
      } else if (format == "html") {
        exportToHTML(results, filename)
        exported_files[["html"]] = filename
      } else {
        warning("Unsupported format: ", format)
      }
    }, error = function(e) {
      warning(paste("Error exporting to", format, "format:", e$message))
    })
  }
  
  # Export metadata if requested
  if (include_metadata) {
    metadata = list(
      export_time = Sys.time(),
      session_info = sessionInfo(),
      formats_exported = formats,
      base_filename = base_filename
    )
    
    metadata_file = file.path(output_dir, paste0(base_filename, "_metadata.json"))
    jsonlite::write_json(metadata, metadata_file, pretty = TRUE)
    exported_files[["metadata"]] = metadata_file
  }
  
  return(exported_files)
}

#' Export to HTML Format
#' 
#' Export results to HTML format with formatting.
#' 
#' @param results Results object
#' @param filename Output filename
#' @param title HTML page title
#' @return HTML file path
#' @export
exportToHTML = function(results, filename, title = "Analysis Results") {
  html_content = c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    paste0("<title>", title, "</title>"),
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 20px; }",
    "table { border-collapse: collapse; width: 100%; margin: 10px 0; }",
    "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
    "th { background-color: #f2f2f2; }",
    "h1, h2 { color: #333; }",
    ".summary { background-color: #f9f9f9; padding: 15px; margin: 10px 0; border-radius: 5px; }",
    "</style>",
    "</head>",
    "<body>",
    paste0("<h1>", title, "</h1>"),
    paste0("<p><strong>Generated:</strong> ", Sys.time(), "</p>")
  )
  
  # Add content based on results type
  if (is.data.frame(results)) {
    html_content = c(html_content, "<h2>Results Table</h2>")
    html_content = c(html_content, "<table>")
    
    # Add header
    html_content = c(html_content, "<tr>")
    for (col in names(results)) {
      html_content = c(html_content, paste0("<th>", col, "</th>"))
    }
    html_content = c(html_content, "</tr>")
    
    # Add data rows
    for (i in 1:nrow(results)) {
      html_content = c(html_content, "<tr>")
      for (j in 1:ncol(results)) {
        html_content = c(html_content, paste0("<td>", results[i, j], "</td>"))
      }
      html_content = c(html_content, "</tr>")
    }
    
    html_content = c(html_content, "</table>")
  } else if (is.list(results)) {
    html_content = c(html_content, "<h2>Results Summary</h2>")
    
    for (name in names(results)) {
      html_content = c(html_content, paste0("<h3>", name, "</h3>"))
      
      if (is.data.frame(results[[name]])) {
        html_content = c(html_content, "<table>")
        
        # Add header
        html_content = c(html_content, "<tr>")
        for (col in names(results[[name]])) {
          html_content = c(html_content, paste0("<th>", col, "</th>"))
        }
        html_content = c(html_content, "</tr>")
        
        # Add data rows (limit to first 10 rows for display)
        n_rows = min(10, nrow(results[[name]]))
        for (i in 1:n_rows) {
          html_content = c(html_content, "<tr>")
          for (j in 1:ncol(results[[name]])) {
            html_content = c(html_content, paste0("<td>", results[[name]][i, j], "</td>"))
          }
          html_content = c(html_content, "</tr>")
        }
        
        html_content = c(html_content, "</table>")
        
        if (nrow(results[[name]]) > 10) {
          html_content = c(html_content, paste0("<p><em>Showing first 10 of ", nrow(results[[name]]), " rows</em></p>"))
        }
      } else {
        html_content = c(html_content, paste0("<div class='summary'>", toString(results[[name]]), "</div>"))
      }
    }
  } else {
    html_content = c(html_content, paste0("<div class='summary'>", toString(results), "</div>"))
  }
  
  html_content = c(html_content, "</body>", "</html>")
  
  writeLines(html_content, filename)
  return(filename)
}

#' Create Analysis Report
#' 
#' Create a comprehensive analysis report in various formats.
#' 
#' @param analysis_results List of analysis results
#' @param report_config Report configuration
#' @param output_dir Output directory
#' @return Report file paths
#' @export
createAnalysisReport = function(analysis_results, report_config, output_dir = ".") {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  report_files = list()
  
  # Generate report content
  report_content = generateReportContent(analysis_results, report_config)
  
  # Export in different formats
  if ("html" %in% report_config$formats) {
    html_file = file.path(output_dir, paste0(report_config$filename, ".html"))
    exportToHTML(report_content, html_file, report_config$title)
    report_files[["html"]] = html_file
  }
  
  if ("pdf" %in% report_config$formats) {
    pdf_file = file.path(output_dir, paste0(report_config$filename, ".pdf"))
    exportToPDF(report_content, pdf_file, report_config$title)
    report_files[["pdf"]] = pdf_file
  }
  
  if ("rds" %in% report_config$formats) {
    rds_file = file.path(output_dir, paste0(report_config$filename, ".rds"))
    saveRDS(report_content, rds_file)
    report_files[["rds"]] = rds_file
  }
  
  return(report_files)
}

#' Generate Report Content
#' 
#' Generate content for analysis report.
#' 
#' @param analysis_results Analysis results
#' @param report_config Report configuration
#' @return Report content
#' @export
generateReportContent = function(analysis_results, report_config) {
  content = list(
    title = report_config$title,
    timestamp = Sys.time(),
    summary = list(),
    details = list()
  )
  
  # Generate summary statistics
  for (result_name in names(analysis_results)) {
    result = analysis_results[[result_name]]
    
    if (is.data.frame(result)) {
      # Calculate summary statistics for data frames
      numeric_cols = sapply(result, is.numeric)
      if (any(numeric_cols)) {
        summary_stats = lapply(result[, numeric_cols, drop = FALSE], function(x) {
          c(mean = mean(x, na.rm = TRUE),
            median = median(x, na.rm = TRUE),
            sd = sd(x, na.rm = TRUE),
            min = min(x, na.rm = TRUE),
            max = max(x, na.rm = TRUE),
            n = length(x),
            n_missing = sum(is.na(x)))
        })
        
        content$summary[[result_name]] = summary_stats
      }
    }
    
    content$details[[result_name]] = result
  }
  
  return(content)
}

#' Export to PDF Format
#' 
#' Export results to PDF format using rmarkdown.
#' 
#' @param results Results object
#' @param filename Output filename
#' @param title PDF title
#' @return PDF file path
#' @export
exportToPDF = function(results, filename, title = "Analysis Results") {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown package required for PDF export")
  }
  
  # Create temporary Rmd file
  temp_rmd = tempfile(fileext = ".Rmd")
  
  rmd_content = c(
    "---",
    paste0("title: ", title),
    "output: pdf_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "```",
    "",
    paste0("**Generated:** ", Sys.time()),
    "",
    "```{r results, echo=FALSE}"
  )
  
  # Add results processing
  if (is.data.frame(results)) {
    rmd_content = c(rmd_content, 
                   "knitr::kable(results, caption = 'Analysis Results')")
  } else if (is.list(results)) {
    for (name in names(results)) {
      rmd_content = c(rmd_content,
                     paste0("cat('### ", name, "\\n')"),
                     paste0("if(is.data.frame(results[['", name, "']])) {"),
                     paste0("  knitr::kable(results[['", name, "']], caption = '", name, "')"),
                     "} else {",
                     paste0("  cat(toString(results[['", name, "']]))"),
                     "}",
                     "cat('\\n\\n')")
    }
  } else {
    rmd_content = c(rmd_content, "cat(toString(results))")
  }
  
  rmd_content = c(rmd_content, "```")
  
  writeLines(rmd_content, temp_rmd)
  
  # Render PDF
  rmarkdown::render(temp_rmd, output_file = filename)
  
  # Clean up temporary file
  unlink(temp_rmd)
  
  return(filename)
}

#' Save Session Information
#' 
#' Save session information for reproducibility.
#' 
#' @param filename Output filename
#' @param include_packages Logical, whether to include package information
#' @param include_platform Logical, whether to include platform information
#' @return Session info file path
#' @export
saveSessionInfo = function(filename, include_packages = TRUE, include_platform = TRUE) {
  session_data = list(
    timestamp = Sys.time(),
    r_version = R.version.string,
    platform = Sys.info()
  )
  
  if (include_packages) {
    session_data$packages = sessionInfo()
  }
  
  if (include_platform) {
    session_data$platform_info = list(
      os = Sys.info()["sysname"],
      release = Sys.info()["release"],
      version = Sys.info()["version"],
      machine = Sys.info()["machine"]
    )
  }
  
  # Save as JSON
  jsonlite::write_json(session_data, filename, pretty = TRUE)
  
  return(filename)
}

#' Export Summary Statistics
#' 
#' Export summary statistics for analysis results.
#' 
#' @param results Results object
#' @param filename Output filename
#' @param statistics Character vector of statistics to calculate
#' @return Summary statistics file path
#' @export
exportSummaryStatistics = function(results, filename, statistics = c("mean", "median", "sd", "min", "max")) {
  summary_data = list()
  
  if (is.data.frame(results)) {
    # Calculate statistics for each numeric column
    numeric_cols = sapply(results, is.numeric)
    if (any(numeric_cols)) {
      for (col in names(results)[numeric_cols]) {
        summary_data[[col]] = calculateStatistics(results[[col]], statistics)
      }
    }
  } else if (is.list(results)) {
    # Process each list element
    for (name in names(results)) {
      if (is.data.frame(results[[name]])) {
        summary_data[[name]] = exportSummaryStatistics(results[[name]], NULL, statistics)
      }
    }
  }
  
  # Save summary data
  if (!is.null(filename)) {
    jsonlite::write_json(summary_data, filename, pretty = TRUE)
  }
  
  return(summary_data)
}

#' Calculate Statistics
#' 
#' Calculate specified statistics for a numeric vector.
#' 
#' @param x Numeric vector
#' @param statistics Character vector of statistics to calculate
#' @return Named vector of statistics
#' @export
calculateStatistics = function(x, statistics = c("mean", "median", "sd", "min", "max")) {
  stats = list()
  
  for (stat in statistics) {
    if (stat == "mean") {
      stats$mean = mean(x, na.rm = TRUE)
    } else if (stat == "median") {
      stats$median = median(x, na.rm = TRUE)
    } else if (stat == "sd") {
      stats$sd = sd(x, na.rm = TRUE)
    } else if (stat == "min") {
      stats$min = min(x, na.rm = TRUE)
    } else if (stat == "max") {
      stats$max = max(x, na.rm = TRUE)
    } else if (stat == "n") {
      stats$n = length(x)
    } else if (stat == "n_missing") {
      stats$n_missing = sum(is.na(x))
    } else if (stat == "quantiles") {
      stats$quantiles = quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    }
  }
  
  return(stats)
} 