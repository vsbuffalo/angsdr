## diagnostics.R -- diagnostic plots for ANGSD data

#' Create a smooth scatter plot with marginal histograms, with optional lowess
#' curve.
#'
#' This function creates a smooth scatter plot, with optional lowess curve,
#' which is useful for discovering possible technical artifacts or biases in
#' data. This function is a lower-level alternative to \code{\link{windowDiagPlot}},
#' which creates a plot of a particular summary statistic column against the
#' number of sites used by ANGSD using data from a \code{GRanges} object
#' containing ANGSD data.
#'
#' @param x numeric vector of x coordinates.
#' @param y numeric vector of y coordinates.
#' @param lowess logical value indicating whether to include a lowess curve.
#' @param xlab title for the x axis.
#' @param ylab title for the y axis.
#' @param main title for the plot.
#'
#' @export
smoothScatterHist <- function(x, y, lowess=TRUE, xlab=NULL, ylab=NULL, main=NULL) {
  # adapted from r-blogger.
  opar <- par(no.readonly=TRUE)
  zones <- matrix(c(1,1,1,
                    0,5,0,
                    2,6,4,
                    0,3,0), ncol = 3, byrow = TRUE)
  layout(zones, widths=c(0.3, 4, 1), heights = c(1, 3, 10, .75))
  # layout.show(n=6)
  xhist <- hist(x, plot = FALSE)
  yhist <- hist(y, plot = FALSE)
  top <- max(c(xhist$counts, yhist$counts))
  # top plot, the main title
  par(xaxt="n", yaxt="n", bty="n", mar = c(.3, 2, .3, 0) + .05)
  plot(x=1, y=1, type="n", ylim=c(-1, 1), xlim=c(-1, 1))
  if (!is.null(main))
    text(0, 0, paste(main), cex=2)

  # y axis label
  plot(x=1, y=1, type="n", ylim=c(-1, 1), xlim=c(-1, 1))
  if (!is.null(ylab))
    text(0, 0, paste(ylab), cex=1.5, srt=90)

  # x axis label
  plot(x=1, y=1, type="n", ylim=c(-1, 1), xlim=c(-1, 1))
  if (!is.null(xlab))
    text(0, 0, paste(xlab), cex=1.5)

  # y axis histogram, without right margins
  par(mar = c(2, 0, 1, 1))
  barplot(yhist$counts, axes = FALSE, xlim = c(0, top),
          space = 0, horiz = TRUE)

  # x axis histogram, without bottom margins
  par(mar = c(0, 2, 1, 1))
  barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0)

  # smooth scatter with lowess
  par(mar = c(2, 2, .5, .5), xaxt="s", yaxt="s", bty="n")
  smoothScatter(x, y)
  if (lowess)
    lines(lowess(x, y), col="red")
  par(opar)
}
