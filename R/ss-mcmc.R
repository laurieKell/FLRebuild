#' Read Stock Synthesis MCMC derived quantities as Kobe draws
#'
#' Converts one file of Stock Synthesis MCMC derived quantities into
#' per-iteration Kobe draws \eqn{(SSB/SSB_{MSY},\, F/F_{MSY})}. Two layouts are
#' accepted:
#' \describe{
#'   \item{wide}{one row per iteration, columns \code{SSB_<year>},
#'     \code{Bratio_<year>}, \code{F_<year>}, \code{SSB_MSY}, \code{annF_MSY}
#'     (the \code{derived_posteriors.sso} layout, also as written by
#'     \code{r4ss::SSgetMCMC});}
#'   \item{long}{columns \code{Label} and \code{Value} (the \code{derived_quants}
#'     layout) plus an iteration column (\code{Iter}, \code{iter} or
#'     \code{iteration}).}
#' }
#'
#' \code{stock = "ssbmsy"} divides \code{SSB_<year>} by \code{SSB_MSY} within
#' each iteration; \code{stock = "bratio"} uses \code{Bratio_<year>} as reported,
#' which is \eqn{SSB/SSB_{MSY}} only when the starter depletion basis is
#' \eqn{B_{MSY}}. \code{harvest = "F"} uses \code{F_<year>} as reported (already
#' \eqn{F/F_{MSY}} when the starter F report basis is 2); \code{harvest = "annF"}
#' divides \code{F_<year>} by \code{annF_MSY}.
#'
#' @param file Path to a \code{.csv} (or whitespace-delimited \code{.sso}) file.
#' @param stock \code{"ssbmsy"} or \code{"bratio"}.
#' @param harvest \code{"F"} or \code{"annF"}.
#' @param ssbLag Years subtracted from the SSB year, e.g. \code{1} to align
#'   start-of-year SSB with end-of-year biomass as in joint SS3/JABBA Kobe plots.
#' @param years Optional years to keep.
#' @return A \code{data.frame} with \code{iter}, \code{year}, \code{stock},
#'   \code{harvest}, \code{ssb}, \code{f}.
#' @seealso \code{\link{getMcmcRuns}}, \code{\link{kobeEmpirical}}
#' @export
readMcmcDq <- function(file, stock = c("ssbmsy", "bratio"),
                       harvest = c("F", "annF"), ssbLag = 0, years = NULL) {

  stock   <- match.arg(stock)
  harvest <- match.arg(harvest)

  dat  <- mcmcTable(file)
  iter <- dat$Iter

  yrCols <- function(prefix) {
    nm <- grep(paste0("^", prefix, "_[0-9]{4}$"), names(dat), value = TRUE)
    stats::setNames(nm, sub(paste0("^", prefix, "_"), "", nm))
  }
  longOf <- function(prefix) {
    nm <- yrCols(prefix)
    if (!length(nm)) stop("no ", prefix, "_<year> columns in ", file)
    data.frame(iter = rep(iter, length(nm)),
               year = rep(as.numeric(names(nm)), each = length(iter)),
               value = unlist(lapply(nm, function(n) as.numeric(dat[[n]])),
                              use.names = FALSE))
  }
  refpt <- function(nm) {
    if (!nm %in% names(dat)) stop("column ", nm, " not found in ", file)
    stats::setNames(as.numeric(dat[[nm]]), iter)
  }

  ssb <- longOf("SSB")
  ssb$year <- ssb$year - ssbLag
  names(ssb)[3] <- "ssb"
  ssb$stock <- if (stock == "ssbmsy")
    ssb$ssb / refpt("SSB_MSY")[as.character(ssb$iter)]
  else {
    br <- longOf("Bratio")
    br$year <- br$year - ssbLag
    br$value[match(paste(ssb$iter, ssb$year), paste(br$iter, br$year))]
  }

  f <- longOf("F")
  names(f)[3] <- "f"
  f$harvest <- if (harvest == "annF")
    f$f / refpt("annF_MSY")[as.character(f$iter)] else f$f

  kb <- merge(ssb, f, by = c("iter", "year"))
  if (!is.null(years)) kb <- kb[kb$year %in% years, , drop = FALSE]
  kb <- kb[order(kb$year, kb$iter), c("iter", "year", "stock", "harvest",
                                      "ssb", "f")]
  rownames(kb) <- NULL
  kb
}

mcmcTable <- function(file) {

  dat <- if (grepl("\\.csv$", file, ignore.case = TRUE))
    utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  else {
    txt <- readLines(file, warn = FALSE)
    txt <- txt[nzchar(trimws(txt))]
    hdr <- grep("^\\s*Iter\\b", txt)
    if (!length(hdr)) hdr <- 1L
    ends <- c(hdr[-1] - 1L, length(txt))
    nrw  <- ends - hdr
    if (!any(nrw > 0)) stop("no data rows in ", file)
    blk  <- max(which(nrw == max(nrw)))
    if (length(hdr) > 1L)
      message(basename(file), ": ", length(hdr), " header blocks with ",
              paste(nrw, collapse = "/"), " rows; using block ", blk)
    txt <- txt[hdr[blk]:ends[blk]]
    nf  <- lengths(strsplit(trimws(txt), "[[:space:]]+"))
    bad <- nf != nf[1]
    if (any(bad)) {
      message(basename(file), ": dropped ", sum(bad), " incomplete row(s)")
      txt <- txt[!bad]
    }
    utils::read.table(text = txt, header = TRUE,
                      check.names = FALSE, stringsAsFactors = FALSE)
  }

  if (all(c("Label", "Value") %in% names(dat))) {
    itc <- intersect(c("Iter", "iter", "iteration", "Iteration"), names(dat))
    if (!length(itc)) dat$Iter <- 1L else dat$Iter <- dat[[itc[1]]]
    dat <- stats::reshape(dat[, c("Iter", "Label", "Value")],
                          idvar = "Iter", timevar = "Label", direction = "wide")
    names(dat) <- sub("^Value\\.", "", names(dat))
  }
  itc <- intersect(c("Iter", "iter", "iteration", "Iteration"), names(dat))
  dat$Iter <- if (length(itc)) dat[[itc[1]]] else seq_len(nrow(dat))
  if (nrow(dat) > 1L && any(dat$Iter == 0))
    dat <- dat[dat$Iter != 0, , drop = FALSE]
  rownames(dat) <- NULL
  dat
}

#' Read Stock Synthesis MCMC forecasts as FLQuants
#'
#' Reads a Stock Synthesis MCMC derived-quantities file
#' (\code{derived_posteriors.sso} from \code{-mceval}, or the same layouts as
#' \code{\link{readMcmcDq}}) into an \code{FLQuants} with one iteration per
#' posterior draw.
#'
#' \code{stock} and \code{harvest} are as in \code{\link{readMcmcDq}}:
#' \eqn{SSB/SSB_{MSY}} and \eqn{F/F_{MSY}} when \code{relative = TRUE}, or
#' \eqn{SSB} and \eqn{F/F_{MSY} \times F_{MSY}} (\code{annF_MSY}) when
#' \code{relative = FALSE}. \code{catch} is \code{ForeCatch_<year>}, the
#' forecast catch in the basis set in \code{forecast.ss} (e.g. dead biomass),
#' so it covers forecast years only; earlier years are \code{NA}. When
#' \code{relative = TRUE} it is divided by MSY within each draw
#' (\code{Dead_Catch_MSY}, else \code{Ret_Catch_MSY} or \code{TotYield_MSY}).
#'
#' @param file Path to the derived-quantities file.
#' @param relative Logical; ratios to MSY reference points (default:
#'   \eqn{SSB/SSB_{MSY}}, \eqn{F/F_{MSY}}, \eqn{C/MSY}) or absolute values.
#' @param stock,harvest Passed to \code{\link{readMcmcDq}}.
#' @param ssbLag Years subtracted from the SSB year; \code{1} gives end-of-year
#'   SSB.
#' @param years Optional years to keep.
#' @return An \code{FLQuants} with \code{stock}, \code{harvest} and
#'   \code{catch}, each with dimensions \code{quant = "all"}, \code{year} and
#'   \code{iter} (one per draw, named by the SS3 iteration number).
#' @examples
#' \dontrun{
#' fq <- readMcmcFLQuants("derived_posteriors.sso", ssbLag = 1)
#' plot(fq)
#' iterMedians(fq[["stock"]])
#' }
#' @seealso \code{\link{readMcmcDq}}, \code{\link{getMcmcRuns}}
#' @export
readMcmcFLQuants <- function(file, relative = TRUE,
                             stock = c("ssbmsy", "bratio"),
                             harvest = c("F", "annF"), ssbLag = 0,
                             years = NULL) {

  stock   <- match.arg(stock)
  harvest <- match.arg(harvest)

  kb  <- readMcmcDq(file, stock = stock, harvest = harvest, ssbLag = ssbLag)
  dat <- mcmcTable(file)

  fc <- grep("^ForeCatch_[0-9]{4}$", names(dat), value = TRUE)
  ct <- if (length(fc))
    data.frame(iter = rep(dat$Iter, length(fc)),
               year = rep(as.numeric(sub("^ForeCatch_", "", fc)),
                          each = nrow(dat)),
               catch = unlist(lapply(fc, function(n) as.numeric(dat[[n]])),
                              use.names = FALSE))
  else data.frame(iter = numeric(0), year = numeric(0), catch = numeric(0))

  if (relative && nrow(ct)) {
    msyCol <- intersect(c("Dead_Catch_MSY", "Ret_Catch_MSY", "TotYield_MSY"),
                        names(dat))
    if (!length(msyCol))
      stop("no MSY catch column (Dead_Catch_MSY, Ret_Catch_MSY, ",
           "TotYield_MSY) in ", file, "; use relative = FALSE")
    msy <- stats::setNames(as.numeric(dat[[msyCol[1]]]), dat$Iter)
    ct$catch <- ct$catch / msy[as.character(ct$iter)]
  }

  if (!relative) {
    if (!"annF_MSY" %in% names(dat))
      stop("annF_MSY not found in ", file, "; use relative = TRUE")
    fmsy <- stats::setNames(as.numeric(dat$annF_MSY), dat$Iter)
    kb$stock   <- kb$ssb
    kb$harvest <- if (harvest == "annF") kb$f
                  else kb$f * fmsy[as.character(kb$iter)]
  }

  yrs <- sort(unique(c(kb$year, ct$year)))
  if (!is.null(years)) yrs <- yrs[yrs %in% years]
  its <- sort(unique(dat$Iter))

  toFLQ <- function(d, col, units) {
    a <- array(NA_real_, dim = c(1, length(yrs), 1, 1, 1, length(its)),
               dimnames = list(quant = "all", year = yrs, unit = "unique",
                               season = "all", area = "unique",
                               iter = its))
    d <- d[d$year %in% yrs & d$iter %in% its, , drop = FALSE]
    a[cbind(1, match(d$year, yrs), 1, 1, 1, match(d$iter, its))] <- d[[col]]
    FLCore::FLQuant(a, units = units)
  }

  FLCore::FLQuants(
    stock   = toFLQ(kb, "stock",   if (relative) "SSB/SSBmsy" else "SSB"),
    harvest = toFLQ(kb, "harvest", if (relative) "F/Fmsy" else "f"),
    catch   = toFLQ(ct, "catch",   if (relative) "C/MSY" else "t"))
}

#' Load Stock Synthesis MCMC derived quantities by TAC
#'
#' Reads one derived-quantities file per catch scenario with
#' \code{\link{readMcmcDq}} and binds them with a numeric \code{catch} column,
#' giving a \code{kbtrj}-style table that feeds
#' \code{\link{kobeMomentGridJabba}}, \code{\link{makeK2SM}} and
#' \code{\link{kobeEmpirical}}.
#'
#' @param object Directory of files, or a named character vector of paths
#'   (names = catch levels).
#' @param pattern Regular expression with one capture group for the catch,
#'   matched against paths relative to \code{object} (subfolders included).
#'   Default matches e.g. \code{derived_quant_1200.csv}, \code{dq_1200t.csv};
#'   use \code{"TAC([0-9]+)/derived_posteriors\\\\.sso$"} for one SS run folder
#'   per TAC.
#' @param ... Passed to \code{\link{readMcmcDq}}.
#' @return A \code{data.frame} with \code{catch}, \code{iter}, \code{year},
#'   \code{stock}, \code{harvest}, \code{ssb}, \code{f}.
#' @examples
#' \dontrun{
#' kb <- getMcmcRuns("P:/.../MCMC/derived_quants", ssbLag = 0)
#' k2 <- kobeEmpirical(kb, years = c(2030, 2040, 2050))
#' }
#' @export
getMcmcRuns <- function(object, pattern = "_([0-9]+)t?\\.(csv|sso)$", ...) {

  files <- if (length(object) == 1L && dir.exists(object)) {
    f <- list.files(object, recursive = TRUE)
    f <- f[grepl(pattern, f, ignore.case = TRUE)]
    if (!length(f)) stop("no files matching '", pattern, "' in ", object)
    stats::setNames(file.path(object, f),
                    sub(paste0(".*?", pattern, ".*"), "\\1", f,
                        ignore.case = TRUE, perl = TRUE))
  } else object

  catch <- suppressWarnings(as.numeric(names(files)))
  if (is.null(names(files)) || any(is.na(catch)))
    stop("names of 'object' must be catch levels coercible to numeric")

  kb <- do.call(rbind, lapply(seq_along(files), function(i)
    cbind(catch = catch[i], readMcmcDq(files[[i]], ...))))
  kb[order(kb$catch, kb$year, kb$iter), , drop = FALSE]
}

#' Empirical Kobe probabilities from MCMC draws
#'
#' Counts draws in each Kobe quadrant by scenario and year, without any
#' distributional assumption. Quadrants: green (\eqn{B \ge 1, F < 1}), orange
#' (\eqn{B \ge 1, F \ge 1}), yellow (\eqn{B < 1, F < 1}), red
#' (\eqn{B < 1, F \ge 1}), with \eqn{B = SSB/SSB_{MSY}} and
#' \eqn{F = F/F_{MSY}}.
#'
#' @param kb Draws with \code{year}, \code{stock}, \code{harvest} and a scenario
#'   column \code{catch} or \code{run}.
#' @param years Years to tabulate; default all.
#' @param crit Biomass thresholds for \eqn{P(B < crit)}.
#' @param probs Quantiles reported for \code{stock} and \code{harvest}.
#' @return A \code{data.frame} by scenario and year with quantiles,
#'   \code{green}, \code{orange}, \code{yellow}, \code{red},
#'   \code{pNoOverfishing} (\eqn{P(F \le F_{MSY})}), \code{pNotOverfished}
#'   (\eqn{P(B \ge B_{MSY})}) and one \code{pCrit<100*crit>} column per
#'   threshold.
#' @seealso \code{\link{kobeEmpiricalCI}}, \code{\link{makeK2SM}}
#' @export
kobeEmpirical <- function(kb, years = NULL, crit = c(0.1, 0.2),
                          probs = c(0.025, 0.5, 0.975)) {

  scen <- if ("catch" %in% names(kb)) "catch" else if ("run" %in% names(kb))
    "run" else { kb$run <- "run"; "run" }
  if (!is.null(years)) kb <- kb[kb$year %in% years, , drop = FALSE]
  kb <- kb[is.finite(kb$stock) & is.finite(kb$harvest), , drop = FALSE]

  grp <- split(kb, list(kb[[scen]], kb$year), drop = TRUE)
  out <- do.call(rbind, lapply(grp, function(d) {
    b <- d$stock; h <- d$harvest
    qs <- stats::quantile(b, probs, names = FALSE)
    qh <- stats::quantile(h, probs, names = FALSE)
    res <- data.frame(year = d$year[1], n = nrow(d),
                      t(stats::setNames(qs, paste0("stock", 100 * probs))),
                      t(stats::setNames(qh, paste0("harvest", 100 * probs))),
                      green  = mean(b >= 1 & h < 1),
                      orange = mean(b >= 1 & h >= 1),
                      yellow = mean(b < 1 & h < 1),
                      red    = mean(b < 1 & h >= 1),
                      pNoOverfishing = mean(h <= 1),
                      pNotOverfished = mean(b >= 1),
                      check.names = FALSE)
    for (cr in crit) res[[paste0("pCrit", 100 * cr)]] <- mean(b < cr)
    res[[scen]] <- d[[scen]][1]
    res[c(scen, setdiff(names(res), scen))]
  }))
  out <- out[order(out[[scen]], out$year), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Catch achieving a target probability from tabulated probabilities
#'
#' Empirical counterpart of \code{\link{kobeTargetCatch}}: linear interpolation
#' over catch of a probability column (default \code{green}) in a given year,
#' returning the highest catch at which the target is still met. Works on
#' \code{\link{kobeEmpirical}} output or the \code{long} element of
#' \code{\link{makeK2SM}}.
#'
#' @param prob A \code{data.frame} with \code{catch}, \code{year} and \code{col}.
#' @param year Year by which the target must be met.
#' @param target Target probability.
#' @param col Probability column.
#' @return Numeric catch. \code{NA_real_} if the target is not met at the lowest
#'   catch; the highest catch if it is met throughout.
#' @seealso \code{\link{kobeTargetCatch}}, \code{\link{kobeEmpirical}}
#' @export
kobeTargetCatchEmpirical <- function(prob, year, target = 0.6, col = "green") {

  d <- prob[prob$year == year & is.finite(prob[[col]]), , drop = FALSE]
  d <- d[order(d$catch), , drop = FALSE]
  if (nrow(d) < 2L) return(NA_real_)
  p <- d[[col]] - target
  if (p[1] < 0) return(NA_real_)
  if (all(p >= 0)) return(max(d$catch))
  i <- which(p < 0)[1]
  d$catch[i - 1] + (d$catch[i] - d$catch[i - 1]) * p[i - 1] / (p[i - 1] - p[i])
}

#' Empirical nested intervals from MCMC draws
#'
#' Equal-tailed posterior intervals for \code{stock} or \code{harvest} by
#' scenario and year, in the same long form as \code{\link{lognormalCI}}
#' (\code{level}, \code{ratio}, \code{ratioLo}, \code{ratioHi}).
#'
#' @param kb Draws as for \code{\link{kobeEmpirical}}.
#' @param quant \code{"stock"} or \code{"harvest"}.
#' @param levels Coverage probabilities.
#' @return A \code{data.frame} with the scenario column, \code{year},
#'   \code{level}, \code{ratio} (median), \code{ratioLo}, \code{ratioHi}.
#' @export
kobeEmpiricalCI <- function(kb, quant = c("stock", "harvest"),
                            levels = c(0.5, 0.8, 0.95)) {

  quant <- match.arg(quant)
  scen <- if ("catch" %in% names(kb)) "catch" else if ("run" %in% names(kb))
    "run" else { kb$run <- "run"; "run" }
  kb <- kb[is.finite(kb[[quant]]), , drop = FALSE]

  grp <- split(kb, list(kb[[scen]], kb$year), drop = TRUE)
  out <- do.call(rbind, lapply(grp, function(d) {
    x <- d[[quant]]
    res <- do.call(rbind, lapply(levels, function(l)
      data.frame(year = d$year[1], level = l,
                 ratio = stats::median(x),
                 ratioLo = stats::quantile(x, (1 - l) / 2, names = FALSE),
                 ratioHi = stats::quantile(x, 1 - (1 - l) / 2, names = FALSE))))
    res[[scen]] <- d[[scen]][1]
    res
  }))
  out <- out[order(out[[scen]], out$year, out$level), , drop = FALSE]
  rownames(out) <- NULL
  out
}
