\name{overlapWidths}
\alias{overlapWidths}
\title{Count overlapping bases between two \code{GRanges} objects}
\usage{
overlapWidths(x, anno, type = c("any", "start", "end", "within", "equal"))
}
\arguments{
  \item{anno}{a \code{GenomicRanges} object contaning
  annotation ranges.}

  \item{angsd}{a \code{GenomicRanges} object containing
  windows and summary statistics from ANGSD.}

  \item{type}{the type of overlap to do; this passed
  directly to \code{\link{findOverlaps}}.}
}
\description{
\code{overlapWidths} returns a vector containing the number
of bases each range in a \code{GRanges} object overlap with
a another \code{GRanges} object (e.g. one of annotation
exons). For example, if \code{x} is a \code{GRanges} object
containing ANGSD windows, and \code{anno} is a
\code{GRanges} object containing all exon ranges,
\code{overlapWidths(x, anno)} will return a vector that can
be used as \code{elementMetadata} in \code{x} to indicate
how many bases in each range overlap \code{anno}. Note that
\code{overlapWidths} ignores strand.
}
\examples{
test_file <- system.file("inst", "extdata", "test.pestPG", package="angsdr")
d <- readPestPG(test_file)
anno <- GRanges("chr1", IRanges(c(70000, 120300), c(70300, 120800)))
# create overlaps annotation
overlapWidths(d, anno)
# store in d
d$anno_overlap <- overlapWidths(d, anno)
}

