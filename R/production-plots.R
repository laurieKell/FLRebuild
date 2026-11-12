#' Kobe regions on a production-function plot
#'
#' Polygons for catch/MSY against B/BMSY: red is B < BMSY and catch > MSY;
#' green is B > BMSY and catch < MSY; orange is B > BMSY above the FMSY line
#' (catch/MSY = B/BMSY); yellow otherwise.
#'
#' @param max Upper limit of the polygons; clip with \code{coord_cartesian}.
#' @return A \code{data.frame} with \code{quadrant}, \code{x} and \code{y};
#'   draw with \code{geom_polygon(aes(x, y, group = quadrant, fill = quadrant))}.
#'   The right-hand yellow region is \code{"yellow2"}.
#' @export
kobePhasePolygons <- function(max = 100) {
  data.frame(
    quadrant = rep(c("red", "yellow", "green", "yellow2", "orange"), c(4, 4, 4, 3, 3)),
    x = c(0, 1, 1, 0,   0, 1, 1, 0,   1, max, max, 1,   1, max, max,   1, max, 1),
    y = c(1, 1, max, max,   0, 0, 1, 1,   0, 0, 1, 1,   1, 1, max,   1, max, max))
}

#' Trend with terminal-year marginal density
#'
#' Median, 50% and 90% ribbons by scenario with the distribution in the last
#' year in a right-hand margin, shaded below and above 1.
#'
#' @param q \code{data.frame} with \code{scenario}, \code{year}, \code{q05},
#'   \code{q25}, \code{q50}, \code{q75} and \code{q95}.
#' @param term \code{data.frame} of terminal-year draws with \code{scenario}
#'   and the column \code{var}.
#' @param var Column of \code{term} to draw in the margin.
#' @param ylab Y-axis label.
#' @param below,above Fill below and above 1.
#' @param col Named colours by scenario.
#' @param cap Upper y limit.
#' @param worms Optional individual trajectories: \code{scenario}, \code{year},
#'   \code{iter} and \code{value}. Drawn over the ribbons and under the median.
#'   \code{NULL} leaves the ribbons alone.
#' @return A patchwork object.
#' @export
plotTrendMar <- function(q, term, var, ylab, below, above, col, cap = 5,
                         worms = NULL) {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
  fill <- ggplot2::scale_fill_manual(values = col)
  base <- list(ggplot2::facet_grid(scenario ~ .),
               ggplot2::coord_cartesian(ylim = c(0, cap), expand = FALSE),
               fill, theme_kobe(), ggplot2::theme(legend.position = "none"))
  trend <- ggplot2::ggplot(q, ggplot2::aes(year)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, fill = below, alpha = 0.3) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = cap, fill = above, alpha = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q95, fill = scenario), alpha = 0.4) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75, fill = scenario), alpha = 0.8)
  if (!is.null(worms) && nrow(worms))
    trend <- trend + ggplot2::geom_line(
      ggplot2::aes(year, value, group = iter),
      data = worms, inherit.aes = FALSE, colour = "grey15",
      alpha = 0.45, linewidth = 0.3, na.rm = TRUE)
  trend <- trend +
    ggplot2::geom_line(ggplot2::aes(y = q50), linewidth = 0.6) +
    base + ggplot2::labs(x = NULL, y = ylab) +
    ggplot2::theme(strip.text.y = ggplot2::element_blank())
  marg <- ggplot2::ggplot(term, ggplot2::aes(y = .data[[var]], fill = scenario)) +
    ggplot2::geom_density(ggplot2::aes(x = ggplot2::after_stat(scaled)), colour = "grey40") +
    ggplot2::geom_hline(yintercept = 1, colour = "red") +
    base + ggplot2::labs(x = max(q$year), y = NULL) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  patchwork::wrap_plots(trend, marg, widths = c(5, 1))
}

#' Four-panel process-error plot
#'
#' Rows: process error through time with \code{\link{rod}} regimes; against
#' stock/BMSY with a loess smooth; autocorrelation with 95\% bounds;
#' histogram. One column per scenario.
#'
#' @param pe Output of \code{\link{surplusError}} for a named list of fits.
#' @param col Named colours by scenario.
#' @param lag.max Maximum ACF lag.
#' @return A patchwork object.
#' @export
plotSurplusError <- function(pe, col, lag.max = 20) {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
  pe  <- pe[is.finite(pe$pe), ]
  rg  <- plyr::ddply(pe, "scenario", function(x) rod(x$pe, year = x$year))
  acs <- plyr::ddply(pe, "scenario", function(x) {
    a <- stats::acf(x$pe, lag.max = lag.max, plot = FALSE)
    data.frame(lag = a$lag[-1], acf = a$acf[-1], ci = stats::qnorm(0.975) / sqrt(a$n.used))
  })
  base <- list(ggplot2::facet_wrap(~scenario, nrow = 1, scales = "free"),
               ggplot2::scale_fill_manual(values = col), ggplot2::scale_colour_manual(values = col),
               theme_kobe(), ggplot2::theme(legend.position = "none"))
  zero <- ggplot2::geom_hline(yintercept = 0, colour = "grey60")
  patchwork::wrap_plots(ncol = 1,
    ggplot2::ggplot(pe, ggplot2::aes(year, pe)) +
      ggplot2::geom_polygon(ggplot2::aes(year, data, group = regime, fill = scenario), data = rg, alpha = 0.4) +
      zero + ggplot2::geom_line(linewidth = 0.4) + ggplot2::geom_point(size = 0.7) +
      ggplot2::labs(x = NULL, y = "Process error / MSY") + base,
    ggplot2::ggplot(pe, ggplot2::aes(x, pe)) +
      zero + ggplot2::geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(colour = scenario), size = 0.8) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "grey15", linewidth = 0.6) +
      ggplot2::labs(x = expression(B / B[MSY]), y = "Process error / MSY") + base,
    ggplot2::ggplot(acs, ggplot2::aes(lag, acf)) +
      zero + ggplot2::geom_ribbon(ggplot2::aes(ymin = -ci, ymax = ci), fill = "grey80") +
      ggplot2::geom_col(ggplot2::aes(fill = scenario), width = 0.6) +
      ggplot2::labs(x = "Lag", y = "ACF") + base,
    ggplot2::ggplot(pe, ggplot2::aes(pe)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = scenario), bins = 25, colour = "grey40") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey30") +
      ggplot2::labs(x = "Process error / MSY", y = "Count") + base)
}
