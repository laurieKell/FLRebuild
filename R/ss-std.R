#' Read Stock Synthesis delta-method output from ADMB std and cor files
#'
#' Rebuilds the \code{derived_quants} and \code{CoVar} elements of an
#' \code{\link[r4ss]{SS_output}} object from \code{ss3.std} (or \code{ss.std})
#' and \code{ss3.cor}, for runs whose \code{Report.sso} and \code{covar.sso}
#' were later overwritten, e.g. by \code{-mceval} in an MCMC folder. Labels are
#' taken from the \code{DERIVED_QUANTITIES} section of \code{Report.sso} and
#' matched to the \code{SSB_std}, \code{F_std} and \code{depletion} sdreport
#' vectors in order. The result feeds \code{\link{kobeMomentGrid}}.
#'
#' @param dir Stock Synthesis run directory.
#' @param std,cor,report File names within \code{dir}; \code{std} and
#'   \code{cor} default to \code{ss3.*}, falling back to \code{ss.*}.
#' @return A list with \code{derived_quants} (\code{Label}, \code{Value},
#'   \code{StdDev}) for \code{SSB_}, \code{Recr_}, \code{F_} and \code{Bratio_}, plus the
#'   management quantities from \code{Report.sso} without standard errors, and
#'   \code{CoVar} (\code{label.i}, \code{label.j}, \code{corr}) for all
#'   \code{F_} and \code{Bratio_} pairs.
#' @examples
#' \dontrun{
#' mle  <- readSsStd("SMA North_MCMC_0")
#' grid <- kobeMomentGrid(list("0" = mle), years = 1950:2070, bOffset = 1)
#' }
#' @seealso \code{\link{kobeMomentGrid}}, \code{\link{readMcmcDq}}
#' @export
readSsStd <- function(dir, std = NULL, cor = NULL, report = "Report.sso") {

  pick <- function(f, ext) {
    if (!is.null(f)) return(file.path(dir, f))
    f <- file.path(dir, paste0(c("ss3.", "ss."), ext))
    f <- f[file.exists(f)]
    if (!length(f)) stop("no ss3.", ext, " or ss.", ext, " in ", dir)
    f[1]
  }
  sd <- utils::read.table(pick(std, "std"), header = TRUE,
                          stringsAsFactors = FALSE)

  rep <- readLines(file.path(dir, report), warn = FALSE)
  i0  <- grep("^DERIVED_QUANTITIES", rep)
  if (!length(i0)) stop("no DERIVED_QUANTITIES section in ", report)
  i0  <- i0[1] + grep("^Label", rep[(i0[1] + 1):length(rep)])[1] + 1
  i1  <- i0 - 2 + which(!nzchar(trimws(rep[i0:length(rep)])))[1]
  dqR <- utils::read.table(text = rep[i0:i1], fill = TRUE,
                           stringsAsFactors = FALSE)[, 1:2]
  names(dqR) <- c("Label", "Value")

  block <- function(sdName, pattern) {
    s  <- sd[sd$name == sdName, ]
    lb <- grep(pattern, dqR$Label, value = TRUE)
    if (length(lb) != nrow(s))
      stop(sdName, ": ", nrow(s), " values in std but ", length(lb),
           " labels in ", report)
    data.frame(Label = lb, index = s$index, Value = s$value,
               StdDev = s$std.dev, stringsAsFactors = FALSE)
  }
  dq <- rbind(block("SSB_std",   "^SSB_(Virgin|Initial|[0-9]{4})$"),
              block("recr_std",  "^Recr_(Virgin|Initial|[0-9]{4})$"),
              block("F_std",     "^F_[0-9]{4}$"),
              block("depletion", "^Bratio_[0-9]{4}$"))

  txt  <- readLines(pick(cor, "cor"), warn = FALSE)
  rows <- strsplit(trimws(txt[grepl("^\\s*[0-9]+\\s", txt)]), "\\s+")
  idx  <- as.integer(vapply(rows, `[`, "", 1))
  want <- dq$index[grepl("^(F|Bratio)_", dq$Label)]
  lab  <- dq$Label[match(want, dq$index)]
  R <- matrix(NA_real_, length(want), length(want))
  for (k in seq_along(want)) {
    v <- as.numeric(rows[[match(want[k], idx)]][-(1:4)])
    lo <- want <= want[k]
    R[k, lo] <- v[want[lo]]
  }
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  pr <- which(lower.tri(R), arr.ind = TRUE)

  mgmt <- dqR[!dqR$Label %in% dq$Label, ]
  list(derived_quants = rbind(dq[, c("Label", "Value", "StdDev")],
                              data.frame(Label = mgmt$Label,
                                         Value = suppressWarnings(as.numeric(mgmt$Value)),
                                         StdDev = NA_real_)),
       CoVar = data.frame(label.i = lab[pr[, 1]], label.j = lab[pr[, 2]],
                          corr = R[pr], stringsAsFactors = FALSE))
}
