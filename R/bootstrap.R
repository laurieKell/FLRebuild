#' Bootstrap Analysis Functions for FLRebuild Package
#' 
#' This file contains functions for performing bootstrap analysis on FLBRP objects
#' with parallel processing capabilities.

#' Bootstrap FLBRP Objects
#' 
#' Perform bootstrap analysis on FLBRP objects with parallel processing.
#' 
#' @param flbrp_list List of FLBRP objects to bootstrap
#' @param nits Number of bootstrap iterations (default=100)
#' @param ftmb_args List of arguments to pass to ftmb3 function
#' @param prior_data Data frame with stock-specific priors (optional)
#' @param ncores Number of cores to use (default=detectCores()-4)
#' @return List of bootstrapped FLBRP objects
#' @export
bootstrapFLBRP = function(flbrp_list, nits = 100, ftmb_args = list(), prior_data = NULL, ncores = NULL) {
  # Setup parallel processing
  if (is.null(ncores)) {
    ncores = detectCores() - 4
  }
  cl = makeCluster(ncores)
  registerDoParallel(cl)
  ids = names(flbrp_list)
  
  # Bootstrap function for single FLBRP
  bootstrapSingleFLBRP = function(id, flbrp_list, nits, ftmb_args, prior_data) {
    tryCatch({
      sr = attributes(flbrp_list[[id]])$sr
      ssbBoot = propagate(ssb(sr), nits)
      recBoot = propagate(rec(sr), nits)
      yrs = dimnames(ssbBoot)$year
      nYrs = dims(ssbBoot)$year
      smpl = sample(yrs, nYrs * nits, TRUE)
      ssbBoot[] = ssb(sr)[, smpl]
      ssb(sr) = ssbBoot
      recBoot[] = rec(sr)[, ac(an(smpl) + dims(rec(sr))$minyear - dims(ssb(sr))$minyear)]
      rec(sr) = recBoot

      # Handle stock-specific priors
      ftmb_args_final = ftmb_args
      if (!is.null(prior_data)) {
        stock_prior = prior_data[prior_data$.id == id, ]
        if (nrow(stock_prior) > 0) {
          if ("s" %in% names(stock_prior) && "cv_s" %in% names(ftmb_args)) {
            ftmb_args_final$prior_s = stock_prior$s
          }
          if ("r0" %in% names(stock_prior) && "cv_r0" %in% names(ftmb_args)) {
            ftmb_args_final$prior_r0 = stock_prior$r0
          }
        }
      }
      
      # Set default arguments if not provided
      defaults = list(
        model = "bevholtSV",
        s.est = TRUE,
        s = 0.7,
        s.logitsd = 0.4,
        spr0 = 0.7,
        prior_s = NULL,
        cv_s = NULL,
        prior_r0 = NULL,
        cv_r0 = NULL
      )
      # Merge defaults with provided arguments
      ftmb_args_final = modifyList(defaults, ftmb_args_final)
      
      sr_fitted = do.call(ftmb3, c(list(object = sr), ftmb_args_final))
      rtn = flbrp_list[[id]]
      refpts(rtn) = rbind(
        FLRebuild:::refCreate(c("virgin", "msy", "crash", "spr.30", "spr.20")),
        FLRebuild:::rmax(rtn, 1.0),
        FLRebuild:::rmax(rtn, 0.3),
        FLRebuild:::rmsy(rtn, 0.5),
        FLRebuild:::rvirgin(rtn, 0.3),
        FLRebuild:::refCreate(attributes(rtn)$benchmark)
      )
      params(rtn) = FLPar(sr_fitted[, 1, drop = TRUE])
      refpts(rtn) = computeRefpts(rtn)
      return(refptsEB(rtn))
    }, error = function(e) {
      warning(paste("Error in bootstrap for", id, ":", e$message))
      return(NULL)
    })
  }
  
  results = foreach(id = ids, .packages = c("FLCore", "FLBRP", "haf", "FLCandy")) %dopar% {
    bootstrapSingleFLBRP(id, flbrp_list, nits, ftmb_args, prior_data)
  }
  stopCluster(cl)
  names(results) = ids
  results = results[!sapply(results, is.null)]
  return(results)
}

#' Bootstrap Single FLBRP Object
#' 
#' Perform bootstrap analysis for a single FLBRP object.
#' 
#' @param flbrp Single FLBRP object to bootstrap
#' @param nits Number of bootstrap iterations
#' @param ftmb_args List of arguments to pass to ftmb3 function
#' @param prior_data Data frame with stock-specific priors (optional)
#' @return Bootstrapped FLBRP object
#' @export
bootstrapSingleFLBRP = function(flbrp, nits = 100, ftmb_args = list(), prior_data = NULL) {
  tryCatch({
    sr = attributes(flbrp)$sr
    ssbBoot = propagate(ssb(sr), nits)
    recBoot = propagate(rec(sr), nits)
    yrs = dimnames(ssbBoot)$year
    nYrs = dims(ssbBoot)$year
    smpl = sample(yrs, nYrs * nits, TRUE)
    ssbBoot[] = ssb(sr)[, smpl]
    ssb(sr) = ssbBoot
    recBoot[] = rec(sr)[, ac(an(smpl) + dims(rec(sr))$minyear - dims(ssb(sr))$minyear)]
    rec(sr) = recBoot

    # Handle stock-specific priors
    ftmb_args_final = ftmb_args
    if (!is.null(prior_data)) {
      stock_prior = prior_data[prior_data$.id == names(flbrp), ]
      if (nrow(stock_prior) > 0) {
        if ("s" %in% names(stock_prior) && "cv_s" %in% names(ftmb_args)) {
          ftmb_args_final$prior_s = stock_prior$s
        }
        if ("r0" %in% names(stock_prior) && "cv_r0" %in% names(ftmb_args)) {
          ftmb_args_final$prior_r0 = stock_prior$r0
        }
      }
    }
    
    # Set default arguments if not provided
    defaults = list(
      model = "bevholtSV",
      s.est = TRUE,
      s = 0.7,
      s.logitsd = 0.4,
      spr0 = 0.7,
      prior_s = NULL,
      cv_s = NULL,
      prior_r0 = NULL,
      cv_r0 = NULL
    )
    # Merge defaults with provided arguments
    ftmb_args_final = modifyList(defaults, ftmb_args_final)
    
    sr_fitted = do.call(ftmb3, c(list(object = sr), ftmb_args_final))
    rtn = flbrp
    refpts(rtn) = rbind(
      FLRebuild:::refCreate(c("virgin", "msy", "crash", "spr.30", "spr.20")),
      FLRebuild:::rmax(rtn, 1.0),
      FLRebuild:::rmax(rtn, 0.3),
      FLRebuild:::rmsy(rtn, 0.5),
      FLRebuild:::rvirgin(rtn, 0.3),
      FLRebuild:::refCreate(attributes(rtn)$benchmark)
    )
    params(rtn) = FLPar(sr_fitted[, 1, drop = TRUE])
    refpts(rtn) = computeRefpts(rtn)
    return(refptsEB(rtn))
  }, error = function(e) {
    warning(paste("Error in bootstrap:", e$message))
    return(NULL)
  })
}

#' Setup Parallel Processing
#' 
#' Setup parallel processing environment for bootstrap analysis.
#' 
#' @param ncores Number of cores to use (default=detectCores()-4)
#' @return Parallel processing cluster
#' @export
setupParallelProcessing = function(ncores = NULL) {
  if (is.null(ncores)) {
    ncores = detectCores() - 4
  }
  cl = makeCluster(ncores)
  registerDoParallel(cl)
  return(cl)
}

#' Cleanup Parallel Processing
#' 
#' Clean up parallel processing resources.
#' 
#' @param cl Parallel processing cluster
#' @export
cleanupParallelProcessing = function(cl) {
  if (!is.null(cl)) {
    stopCluster(cl)
  }
} 