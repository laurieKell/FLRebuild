# Stock Synthesis MLE / MCMC forecast runs
#
#   dirs <- ssSetup(scenario, dir, forecasts)   # dir/mle and dir/mcmc, each with forecasts/<label>
#   ssMcmc(dirs["mcmc"])                         # adnuts NUTS; saves the draws (.psv)
#   ssForecast(dirs["mle"])                      # ss -sdonly per forecast (delta method)
#   ssForecast(dirs["mcmc"])                     # ss -mceval per forecast (posterior)

.ssFirst <- function(dir, files) {
  f <- file.path(dir, files)
  f <- f[file.exists(f)]
  if (length(f)) f[1] else NULL
}

.ssExe <- function(dir) {
  exe <- .ssFirst(dir, c("ss3.exe", "ss.exe", "ss3", "ss"))
  if (is.null(exe)) {
    exe <- list.files(dir, pattern = "\\.exe$", full.names = TRUE)
    if (!length(exe)) stop("no SS3 executable in ", dir)
    exe <- exe[1]
  }
  normalizePath(exe)
}

.ssModel <- function(exe) tools::file_path_sans_ext(basename(exe))

.ssModelFiles <- function(dir) {
  st <- readLines(file.path(dir, "starter.ss"), warn = FALSE)
  st <- trimws(sub("#.*", "", st))
  st <- st[nzchar(st)]
  f  <- c("starter.ss", st[1:2], "forecast.ss", "wtatage.ss")
  f[file.exists(file.path(dir, f))]
}

.ssStarter <- function(file, init = NULL, phase = NULL, burn = NULL, thin = NULL, maxyr = NULL) {
  st  <- readLines(file, warn = FALSE)
  set <- function(st, key, val) {
    if (is.null(val)) return(st)
    i <- grep(key, st, ignore.case = TRUE)
    if (!length(i)) stop("'", key, "' not found in ", file)
    st[i[1]] <- sub("^\\s*-?[0-9.]+", as.character(val), st[i[1]])
    st
  }
  st <- set(st, "#\\s*0=use init values", init)
  st <- set(st, "#\\s*Turn off estimation", phase)
  st <- set(st, "#\\s*MCeval burn", burn)
  st <- set(st, "#\\s*MCeval thin", thin)
  st <- set(st, "#\\s*max year for sdreport", maxyr)
  writeLines(st, file)
}

# The executable reads <name>.par / <name>.psv under its compiled name, which
# can differ from its file name (e.g. ss3 renamed to ss.exe), so write both.
.ssCopyAs <- function(from, to, model, ext) {
  from <- normalizePath(from)
  for (n in unique(c(model, "ss", "ss3"))) {
    dest <- file.path(normalizePath(to), paste0(n, ".", ext))
    if (!identical(normalizePath(dest, mustWork = FALSE), from))
      file.copy(from, dest, overwrite = TRUE)
  }
}

.ssLogDone <- function(dir) {
  log <- file.path(dir, "ss_console.log")
  file.exists(log) &&
    any(grepl("Finished running model|Run has completed", readLines(log, warn = FALSE)))
}

.ssDone <- function(dir, method) {
  out <- if (method == "mcmc") file.exists(file.path(dir, "derived_posteriors.sso"))
         else length(list.files(dir, pattern = "\\.cor$")) > 0
  out && .ssLogDone(dir)
}

.ssRun <- function(dir, exe, args, retry = 5) {
  t0  <- Sys.time()
  owd <- setwd(dir)
  on.exit(setwd(owd))
  for (k in seq_len(retry)) {
    st <- tryCatch(system2(file.path(dir, exe), args, stdout = "ss_console.log",
                           stderr = "ss_console.err"),
                   error = function(e) -1L, warning = function(w) -1L)
    if (identical(as.integer(st), 0L)) break
    Sys.sleep(30)
  }
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
}

#' @rdname ssSetup
#' @param dir Output folder; \code{dir/mle} and \code{dir/mcmc} are created.
#' @param forecasts Named character vector of forecast files (names are the
#'   labels, e.g. catch levels), or a folder of \code{forecast_<label>.ss}
#'   files.
#' @param exe SS3 executable; default the one in \code{object}.
#' @examples
#' \dontrun{
#' dirs <- ssSetup("data/ss3/HW2e_HW2d_InitFprior_0t_recDev1", "runs/recDev1",
#'                 forecasts = "data/ss3/HW2e_HW2d_InitFprior_0t_recDev1/forecasts")
#' }
#' @export
setMethod("ssSetup", signature(object = "character"),
  function(object, dir, forecasts, exe = NULL) {

  if (!dir.exists(object)) stop("scenario folder not found: ", object)
  exe <- if (is.null(exe)) .ssExe(object) else normalizePath(exe)
  model <- .ssModel(exe)

  if (length(forecasts) == 1L && dir.exists(forecasts)) {
    f <- list.files(forecasts, pattern = "^forecast_(.+)\\.ss$", full.names = TRUE)
    forecasts <- stats::setNames(f, sub("^forecast_(.+)\\.ss$", "\\1", basename(f)))
  }
  if (!length(forecasts) || is.null(names(forecasts)) || any(!nzchar(names(forecasts))))
    stop("'forecasts' must be a named vector of files or a folder of forecast_<label>.ss")
  if (!all(file.exists(forecasts)))
    stop("forecast files not found: ", paste(forecasts[!file.exists(forecasts)], collapse = ", "))

  par <- .ssFirst(object, c(paste0(model, ".par"), "ss3.par", "ss.par"))
  if (is.null(par)) stop("no .par file in ", object, "; fit the model first")
  hes <- file.path(object, c("admodel.hes", "admodel.cov"))
  hes <- hes[file.exists(hes)]
  files <- file.path(object, .ssModelFiles(object))

  fill <- function(to, forecast = NULL) {
    dir.create(to, recursive = TRUE, showWarnings = FALSE)
    file.copy(c(files, exe, hes), to, overwrite = TRUE)
    .ssCopyAs(par, to, model, "par")
    if (!is.null(forecast)) file.copy(forecast, file.path(to, "forecast.ss"), overwrite = TRUE)
    .ssStarter(file.path(to, "starter.ss"), init = 1)
  }

  dirs <- c(mle = file.path(dir, "mle"), mcmc = file.path(dir, "mcmc"))
  for (d in dirs) {
    fill(d)
    for (l in names(forecasts)) fill(file.path(d, "forecasts", l), forecasts[[l]])
  }
  stats::setNames(normalizePath(dirs), names(dirs))
})

#' @rdname ssMcmc
#' @param chains,iter,warmup,thin Main NUTS run: \code{iter} includes
#'   \code{warmup}; draws saved per chain are \code{(iter - warmup) / thin}.
#' @param pilot \code{list(iter, warmup, thin)} for a pilot run whose mass
#'   matrix is used by the main run; \code{NULL} uses the MLE mass matrix.
#' @param control Passed to \code{adnuts::sample_nuts} with the metric.
#' @param seed Seed for the chain seeds.
#' @param optimise Re-optimise with \code{-hbf 1 -mcmc 15} first: NUTS needs
#'   a Hessian with the hybrid bounding functions, and \code{-mcmc} turns off
#'   the recruitment bias adjustment so the Hessian matches the sampled model.
#' @examples
#' \dontrun{
#' ssMcmc(dirs["mcmc"])
#' }
#' @export
setMethod("ssMcmc", signature(object = "character"),
  function(object, chains = 3, thin = 100, iter = 11000 * thin, warmup = iter / 10,
           pilot = list(iter = 40000, warmup = 8000, thin = 20),
           control = list(adapt_delta = 0.8, max_treedepth = 2),
           seed = 632, optimise = TRUE) {

  if (!requireNamespace("adnuts", quietly = TRUE)) stop("package 'adnuts' is required")
  object <- normalizePath(object)
  exe    <- .ssExe(object)
  model  <- .ssModel(exe)

  if (optimise) {
    owd <- setwd(object)
    st  <- system2(exe, c("-hbf", "1", "-nox", "-iprint", "200", "-mcmc", "15"),
                   stdout = "ss_optimise.log", stderr = "ss_optimise.err")
    setwd(owd)
    if (!file.exists(file.path(object, "admodel.hes")))
      stop("optimisation did not produce admodel.hes; see ss_optimise.log")
  }

  set.seed(seed)
  seeds  <- sample.int(1e4, chains)
  metric <- "mle"
  if (!is.null(pilot)) {
    p <- adnuts::sample_nuts(model, path = object, iter = pilot$iter, warmup = pilot$warmup,
                             thin = pilot$thin, chains = chains, cores = chains, seeds = seeds,
                             control = c(list(metric = "mle"), control))
    saveRDS(p, file.path(object, "pilot.rds"))
    metric <- p$covar.est
  }
  fit <- adnuts::sample_nuts(model, path = object, iter = iter, warmup = warmup, thin = thin,
                             chains = chains, cores = chains, seeds = seeds, mceval = FALSE,
                             control = c(list(metric = metric), control))
  saveRDS(fit, file.path(object, "fit.rds"))

  psv <- .ssFirst(object, c(paste0(model, ".psv"), "ss3.psv", "ss.psv"))
  if (!is.null(psv)) .ssCopyAs(psv, object, model, "psv")
  invisible(fit)
})

#' @rdname ssForecast
#' @param labels Forecast labels to run; default all in \code{object/forecasts}.
#' @param method \code{"auto"} (MCMC if \code{object} has a \code{.psv}),
#'   \code{"mle"} or \code{"mcmc"}.
#' @param jobs Number of runs at once.
#' @param burn,thin Burn-in and thinning applied by \code{-mceval} to the saved
#'   draws.
#' @param phase Last estimation phase for \code{-sdonly}; must make every
#'   parameter in the Hessian active.
#' @param retry Attempts per run, 30 s apart (e.g. while antivirus scans a new
#'   executable).
#' @param rerun Rerun forecasts that have finished.
#' @examples
#' \dontrun{
#' ssForecast(dirs["mle"])
#' ssForecast(dirs["mcmc"], thin = 20)
#' }
#' @export
setMethod("ssForecast", signature(object = "character"),
  function(object, labels = NULL, method = c("auto", "mle", "mcmc"), jobs = 10,
           burn = 10, thin = 20, phase = 10, retry = 5, rerun = FALSE) {

  object <- normalizePath(object)
  model  <- .ssModel(.ssExe(object))
  psv    <- .ssFirst(object, c(paste0(model, ".psv"), "ss3.psv", "ss.psv"))
  method <- match.arg(method)
  if (method == "auto") method <- if (is.null(psv)) "mle" else "mcmc"
  if (method == "mcmc" && is.null(psv)) stop("no .psv in ", object, "; run ssMcmc first")

  fd <- list.dirs(file.path(object, "forecasts"), recursive = FALSE)
  if (!is.null(labels)) fd <- fd[basename(fd) %in% labels]
  todo <- if (rerun) fd else fd[!vapply(fd, .ssDone, logical(1), method = method)]

  for (d in todo) {
    unlink(file.path(d, c("derived_posteriors.sso", "posteriors.sso", "ss_console.log")))
    if (method == "mcmc") {
      .ssCopyAs(psv, d, model, "psv")
      par <- .ssFirst(object, c(paste0(model, ".par"), "ss3.par", "ss.par"))
      .ssCopyAs(par, d, model, "par")
      .ssStarter(file.path(d, "starter.ss"), init = 1, burn = burn, thin = thin, maxyr = -2)
    } else {
      unlink(list.files(d, pattern = "\\.(cor|std)$", full.names = TRUE))
      .ssStarter(file.path(d, "starter.ss"), init = 1, phase = phase)
    }
  }

  mins <- stats::setNames(rep(NA_real_, length(fd)), fd)
  if (length(todo)) {
    exe  <- basename(.ssExe(todo[1]))
    args <- if (method == "mcmc") c("-mceval", "-nox") else c("-sdonly", "-nox")
    cl   <- parallel::makeCluster(min(jobs, length(todo)))
    on.exit(parallel::stopCluster(cl))
    mins[todo] <- unlist(parallel::parLapplyLB(cl, todo, .ssRun, exe = exe, args = args,
                                               retry = retry))
  }
  data.frame(label = basename(fd), dir = fd,
             ok = vapply(fd, .ssDone, logical(1), method = method),
             minutes = round(mins, 1), row.names = NULL)
})
