\name{smoothScatterHist}
\alias{smoothScatterHist}
\title{Create a smooth scatter plot with marginal histograms, with optional lowess
curve.}
\usage{
smoothScatterHist(x, y, lowess = TRUE, xlab = NULL, ylab = NULL,
  main = NULL)
}
\arguments{
  \item{x}{numeric vector of x coordinates.}

  \item{y}{numeric vector of y coordinates.}

  \item{lowess}{logical value indicating whether to include
  a lowess curve.}

  \item{xlab}{title for the x axis.}

  \item{ylab}{title for the y axis.}

  \item{main}{title for the plot.}
}
\description{
This function creates a smooth scatter plot, with optional
lowess curve, which is useful for discovering possible
technical artifacts or biases in data. This function is a
lower-level alternative to \code{\link{windowDiagPlot}},
which creates a plot of a particular summary statistic
column against the number of sites used by ANGSD using data
from a \code{GRanges} object containing ANGSD data.
}

