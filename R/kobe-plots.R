#' Theme and scales for Kobe advice plots
#'
#' \code{theme_kobe} is a plain \code{theme_bw} with no grid and the legend at
#' the bottom. \code{scale_fill_kobe} colours the four Kobe quadrants,
#' \code{scale_fill_k2sm} is the red-to-green probability ramp used for Kobe II
#' strategy matrices, and \code{scale_alpha_ci} shades nested intervals.
#'
#' @param base_size Base font size.
#' @param name Legend title.
#' @param levels Interval coverage levels, as in \code{\link{kobeCI}}.
#' @param alpha Alpha for each level, narrowest first.
#' @param ... Passed to the underlying \code{ggplot2} function.
#' @return A \code{ggplot2} theme or scale.
#' @examples
#' library(ggplot2)
#' d <- data.frame(x = 1:4, quadrant = c("green", "orange", "yellow", "red"))
#' ggplot(d, aes(x, 1, fill = quadrant)) + geom_col() +
#'   scale_fill_kobe() + theme_kobe()
#' @name kobePlotStyle
NULL

#' @rdname kobePlotStyle
#' @export
theme_kobe <- function(base_size = 10, ...) {
  ggplot2::theme_bw(base_size = base_size, ...) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "bottom",
                   strip.background = ggplot2::element_rect(fill = "grey95"))
}

#' @rdname kobePlotStyle
#' @export
scale_fill_kobe <- function(name = NULL, ...) {
  ggplot2::scale_fill_manual(
    values = c(green = "#4DAF4A", orange = "#FF7F00",
               yellow = "#FFD92F", red = "#E41A1C"),
    name = name, ...)
}

#' @rdname kobePlotStyle
#' @export
scale_fill_k2sm <- function(name = "P", ...) {
  ggplot2::scale_fill_gradientn(
    colours = c(grDevices::colorRampPalette(c("red4", "red"))(12),
                grDevices::colorRampPalette(c("yellowgreen", "darkgreen"))(8)),
    limits = c(0, 1), name = name, ...)
}

#' @rdname kobePlotStyle
#' @export
scale_alpha_ci <- function(levels = c(0.5, 0.8, 0.95),
                           alpha = c(0.35, 0.22, 0.12)[seq_along(levels)],
                           name = "Level", ...) {
  ggplot2::scale_alpha_manual(values = stats::setNames(alpha, levels),
                              name = name, ...)
}

#' Intervals of B/BMSY and F/FMSY from moments or draws
#'
#' Stacks the intervals of relative biomass and relative fishing mortality in
#' one long table, so the same plotting code serves every uncertainty method.
#' A moment grid (\code{\link{kobeMomentGrid}}, \code{\link{kobeMomentGridJabba}})
#' gives log-normal intervals via \code{\link{kobeMomentCI}}; draws (columns
#' \code{stock} and \code{harvest}, e.g. from \code{\link{readMcmcDq}},
#' \code{\link{projectJabba}} or \code{\link{jabbaKbtrj}}) give empirical
#' quantiles via \code{\link{kobeEmpiricalCI}}.
#'
#' @param object A moment grid or a table of draws.
#' @param levels Coverage probabilities.
#' @param method Optional label added as a \code{method} column.
#' @return A \code{data.frame} with \code{method} (if given), the scenario
#'   column (\code{catch} or \code{run}), \code{quant} (\code{"B/BMSY"} or
#'   \code{"F/FMSY"}), \code{year}, \code{level}, \code{ratio} (median),
#'   \code{ratioLo} and \code{ratioHi}.
#' @examples
#' \dontrun{
#' ci <- rbind(kobeCI(gridSS, 0.95, "SS3 MVLN"), kobeCI(kbMc, 0.95, "SS3 MCMC"))
#' plotKobeCI(ci, colour = "method")
#' }
#' @seealso \code{\link{plotKobeCI}}
#' @export
kobeCI <- function(object, levels = c(0.5, 0.8, 0.95), method = NULL) {

  moment <- all(c("muB", "muF", "sigmaB", "sigmaF") %in% names(object))
  if (!moment && !all(c("stock", "harvest") %in% names(object)))
    stop("'object' must be a moment grid (muB, muF, sigmaB, sigmaF) ",
         "or draws (stock, harvest)")

  scen <- intersect(c("catch", "run"), names(object))[1]
  if (is.na(scen)) { object$run <- "run"; scen <- "run" }
  keep <- c(scen, "year", "level", "ratio", "ratioLo", "ratioHi")

  one <- function(q, lab) {
    ci <- if (moment) kobeMomentCI(object, q, levels)
          else kobeEmpiricalCI(object, if (q == "B") "stock" else "harvest",
                               levels)
    cbind(quant = lab, ci[, keep])
  }
  out <- rbind(one("B", "B/BMSY"), one("F", "F/FMSY"))
  if (!is.null(method)) out <- cbind(method = method, out)
  rownames(out) <- NULL
  out
}

#' Plot Kobe intervals over time
#'
#' Ribbon-and-median plot of \code{\link{kobeCI}} output. With one method,
#' nested levels are shaded by \code{alpha}; with several, each method has its
#' own colour (use one level).
#'
#' @param ci Output of \code{\link{kobeCI}}, possibly bound across methods.
#' @param colour Column mapped to colour and fill, e.g. \code{"method"};
#'   \code{NULL} for a single series shaded by level.
#' @param facets Facet formula; default \code{quant ~ catch} when \code{ci}
#'   has \code{catch}, otherwise \code{quant ~ .}.
#' @param ylim Upper limits of the y axis by quantity; medians and upper
#'   bounds above them are clipped.
#' @param vline Optional year marking the start of the projection.
#' @param fill Ribbon colour when \code{colour} is \code{NULL}.
#' @return A \code{ggplot}.
#' @export
plotKobeCI <- function(ci, colour = NULL, facets = NULL,
                       ylim = c("B/BMSY" = 3, "F/FMSY" = 5), vline = NULL,
                       fill = "#1b9e77") {

  if (is.null(facets))
    facets <- if ("catch" %in% names(ci)) quant ~ catch else quant ~ .
  ci$quant <- factor(ci$quant, c("B/BMSY", "F/FMSY"))
  cap <- unname(ylim[as.character(ci$quant)])
  cap[is.na(cap)] <- Inf
  ci$ratioHi <- pmin(ci$ratioHi, cap)
  ci$ratio   <- pmin(ci$ratio, cap)
  med <- ci[ci$level == min(ci$level), , drop = FALSE]

  p <- ggplot2::ggplot(ci, ggplot2::aes(.data$year, .data$ratio))
  p <- if (is.null(colour))
    p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ratioLo,
                                          ymax = .data$ratioHi,
                                          alpha = factor(.data$level)),
                             fill = fill) +
      ggplot2::geom_line(data = med, linewidth = 0.6) +
      scale_alpha_ci(sort(unique(ci$level)))
  else
    p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ratioLo,
                                          ymax = .data$ratioHi,
                                          fill = .data[[colour]],
                                          group = .data[[colour]]),
                             alpha = 0.15) +
      ggplot2::geom_line(data = med, ggplot2::aes(colour = .data[[colour]]),
                         linewidth = 0.6)
  p <- p +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    ggplot2::facet_grid(facets, scales = "free_y") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL, fill = NULL) +
    theme_kobe()
  if (!is.null(vline))
    p <- p + ggplot2::geom_vline(xintercept = vline, linetype = "dotted",
                                 colour = "grey50")
  p
}

#' Long form of a Kobe II strategy matrix
#'
#' Stacks the probability columns of \code{\link{makeK2SM}}\code{$long} or
#' \code{\link{kobeEmpirical}} output into one column for faceted plotting.
#'
#' @param x A \code{data.frame} with \code{catch}, \code{year} and the columns
#'   named in \code{panels}.
#' @param panels Named character vector: names are probability columns, values
#'   are panel labels.
#' @param years Optional years to keep.
#' @return A \code{data.frame} with \code{catch}, \code{year}, \code{panel} and
#'   \code{p}, with \code{catch}, \code{year} and \code{panel} as factors.
#' @export
k2smLong <- function(x, panels = c(pNoOverfishing = "P(F <= Fmsy)",
                                   pNotOverfished = "P(B >= Bmsy)",
                                   green = "P(green)"),
                     years = NULL) {
  if (!is.null(years)) x <- x[x$year %in% years, , drop = FALSE]
  extra <- setdiff(names(x), c("catch", "year", names(panels)))
  out <- do.call(rbind, lapply(names(panels), function(nm)
    data.frame(x[, c("catch", "year", intersect(extra, "method")), drop = FALSE],
               panel = panels[[nm]], p = x[[nm]])))
  out$panel <- factor(out$panel, unname(panels))
  out$catch <- factor(out$catch, sort(unique(as.numeric(as.character(out$catch)))))
  out$year  <- factor(out$year, sort(unique(out$year)))
  out
}

#' Plot a Kobe II strategy matrix
#'
#' Tile plot of probabilities by catch and year, labelled in percent.
#'
#' @param x Output of \code{\link{k2smLong}}.
#' @param facet Column to facet by, e.g. \code{"panel"} or \code{"method"}.
#' @param ncol Number of facet columns; default one column for \code{"panel"}
#'   and one row otherwise.
#' @param size Label size.
#' @return A \code{ggplot}.
#' @export
plotK2SM <- function(x, facet = "panel", ncol = NULL, size = 2.2) {
  if (is.null(ncol))
    ncol <- if (facet == "panel") 1 else length(unique(x[[facet]]))
  ggplot2::ggplot(x, ggplot2::aes(.data$year, .data$catch, fill = .data$p)) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f", 100 * .data$p)),
                       size = size, colour = "#00045b") +
    scale_fill_k2sm() +
    ggplot2::scale_y_discrete(limits = rev(levels(x$catch))) +
    ggplot2::facet_wrap(facet, ncol = ncol) +
    ggplot2::labs(x = NULL, y = "Catch (t)") +
    theme_kobe() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
}

#' Log-normal densities from a Kobe moment grid
#'
#' Evaluates the marginal log-normal densities of \eqn{B/B_{MSY}} and
#' \eqn{F/F_{MSY}} implied by each row of a moment grid, without sampling.
#' Rows with zero \eqn{F} have no \eqn{F/F_{MSY}} density.
#'
#' @param grid Output of \code{\link{kobeMomentGrid}} or
#'   \code{\link{kobeMomentGridJabba}}.
#' @param xB,xF Ranges of \eqn{B/B_{MSY}} and \eqn{F/F_{MSY}}.
#' @param n Number of points per curve (log-spaced).
#' @return A \code{data.frame} with the scenario column, \code{year},
#'   \code{quant}, \code{x} and \code{density}.
#' @export
kobeMomentDensity <- function(grid, xB = c(0.2, 3), xF = c(0.05, 4), n = 200) {
  scen <- intersect(c("catch", "run"), names(grid))[1]
  one <- function(g, mu, sig, rng, lab) {
    if (!is.finite(mu) || !is.finite(sig) || sig <= 0) return(NULL)
    x <- exp(seq(log(rng[1]), log(rng[2]), length.out = n))
    data.frame(g[[scen]], year = g$year, quant = lab, x = x,
               density = stats::dlnorm(x, mu, sig))
  }
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    rbind(one(g, g$muB, g$sigmaB, xB, "B/BMSY"),
          one(g, g$muF, g$sigmaF, xF, "F/FMSY"))
  }))
  names(out)[1] <- scen
  out
}

#' Recovery year from tabulated Kobe probabilities
#'
#' Applies \code{\link{kobeRebuildTime}} to each catch in a table of annual
#' probabilities (\code{\link{kobeEmpirical}} output or
#' \code{\link{makeK2SM}}\code{$long}), with no interpolation between catches.
#' The moment-based counterpart, which interpolates, is
#' \code{\link{kobeRebuildCurve}}.
#'
#' @param prob A \code{data.frame} with \code{catch}, \code{year} and \code{col}.
#' @param target Target probabilities.
#' @param from Optional first year to consider (e.g. the first projection year).
#' @param col Probability column.
#' @return A \code{data.frame} with \code{catch} and the columns of
#'   \code{\link{kobeRebuildTime}}, one row per catch and target.
#' @export
kobeRebuildCurveEmpirical <- function(prob, target = c(0.5, 0.6), from = NULL,
                                      col = "green") {
  if (!is.null(from)) prob <- prob[prob$year >= from, , drop = FALSE]
  out <- do.call(rbind, lapply(sort(unique(prob$catch)), function(c0) {
    p0 <- prob[prob$catch == c0, c("year", col)]
    names(p0)[2] <- "green"
    cbind(catch = c0, do.call(rbind, lapply(target, function(t)
      suppressWarnings(kobeRebuildTime(p0, target = t)))))
  }))
  rownames(out) <- NULL
  out
}

#' Subsample draws by scenario
#'
#' Keeps at most \code{n} draws per group, e.g. to give each catch or method
#' equal weight in a Kobe phase plot.
#'
#' @param kb Table of draws.
#' @param n Maximum draws per group.
#' @param by Grouping column(s).
#' @param seed Optional random seed.
#' @return \code{kb} restricted to the sampled rows.
#' @export
kobeSample <- function(kb, n = 800, by = "catch", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  grp <- split(seq_len(nrow(kb)), kb[by], drop = TRUE)
  idx <- unlist(lapply(grp, function(i)
    if (length(i) > n) i[sample.int(length(i), n)] else i), use.names = FALSE)
  kb[sort(idx), , drop = FALSE]
}

#' Local copy of a remote file
#'
#' Copies a file (e.g. a large fit on a network drive) to a local directory
#' once, and again only if its size changes, and returns the local path.
#'
#' @param file Path to the remote file.
#' @param dir Local directory.
#' @return Local path.
#' @export
localCopy <- function(file, dir = tempdir()) {
  if (!file.exists(file)) stop("file not found: ", file)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  loc <- file.path(dir, basename(file))
  if (!file.exists(loc) || file.size(loc) != file.size(file))
    if (!file.copy(file, loc, overwrite = TRUE)) stop("could not copy ", file)
  loc
}
