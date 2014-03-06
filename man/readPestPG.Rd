\name{readPestPG}
\alias{readPestPG}
\title{Read a pestPG file from ANGSD}
\usage{
readPestPG(file, asGRanges = TRUE, physical = FALSE)
}
\arguments{
  \item{file}{the pestPG-formatted file from ANGSD}

  \item{asGRanges}{\code{logical} indicating whether to
  return a \code{GRanges} object. If \code{FALSE}, will
  return a \code{data.frame}.}

  \item{physical}{\code{logical} value indicating whether
  to use the window as the range, or whether to use the
  physical region the analysis included (for \code{GRanges}
  output only).}
}
\description{
This function reads a pestPG file from ANGSD into either a
\code{\link{GRanges}} object or a plain
\code{\link{data.frame}}.
}
\details{
\code{readPestPG} loads a pestPG file into R, as either a
\code{GRanges} object (with \code{asGRanges=TRUE}) or as a
\code{data.frame} (with \code{asGRanges=FALSE}). When a
\code{GRanges} object is returned, the summary statistics
will be stored as \code{elementMetadata} columns.
Additionally, a column named \code{numsites} specifies the
number of bases used by ANGSD in the calculation of the
summary statistics.

ANGSD's pestPG file includes the following estimators of
theta (metadata column names in paranthesis):

\itemize{ \item Watterson's theta (\code{theta_w}) \item
pairwise diversity (\code{theta_pi}) \item Fu and Li
(\code{theta_fuli}) \item Fay's H (\code{theta_fayh}) \item
L (\code{theta_l}) } As well as the follwing neutrality
test statistics: \itemize{ \item Tajima's D (\code{taj_d})
\item Fu and Li F (\code{fuli_f}) \item Fu and Li's D
(\code{fuli_d}) \item Fay's H (\code{fay_h}) \item Zeng's E
(\code{zeng_e}) }
}

