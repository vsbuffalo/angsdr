## annotate.R -- annotate ANGSD GRanges objects with metadata

#' Count overlapping bases between two \code{GRanges} objects
#'
#' \code{overlapWidths} returns a vector containing the number of bases each
#' range in a \code{GRanges} object overlap with a another \code{GRanges}
#' object (e.g. one of annotation exons). For example, if \code{x} is a
#' \code{GRanges} object containing ANGSD windows, and \code{anno} is a
#' \code{GRanges} object containing all exon ranges, \code{overlapWidths(x,
#' anno)} will return a vector that can be used as \code{elementMetadata} in
#' \code{x} to indicate how many bases in each range overlap \code{anno}. Note
#' that \code{overlapWidths} ignores strand.
#'
#' @param anno a \code{GenomicRanges} object contaning annotation ranges.
#' @param angsd a \code{GenomicRanges} object containing windows and summary statistics from ANGSD.
#' @param type the type of overlap to do; this passed directly to \code{\link{findOverlaps}}.
#' @export
#'
#' @examples
#' test_file <- system.file("inst", "extdata", "test.pestPG", package="angsdr")
#' d <- readPestPG(test_file)
#' anno <- GRanges("chr1", IRanges(c(70000, 120300), c(70300, 120800)))
#' # create overlaps annotation
#' overlapWidths(d, anno)
#' # store in d
#' d$anno_overlap <- overlapWidths(d, anno)
overlapWidths <- function (x, anno,
                           type=c("any", "start", "end", "within", "equal")) {
  if (!is(x, "GRanges") || !is(anno, "GRanges"))
    stop("both anno and x need to be GRanges objects")

  ranno <- reduce(anno, ignore.strand=TRUE)

	# find overlaps between angsd windows and annotation
	olaps <- findOverlaps(x, ranno, type=type, ignore.strand=TRUE)

	# extract overlappging ranges and the widths of these overlapping ranges
	rr <- ranges(olaps, ranges(x), ranges(ranno))
	w <- width(rr)

	# there are n widths, now, where n is the number of overlapping ranges (e.g.
	# length(olaps)). We want to group these by where they fall in the ANGSD
	# window, the findOverlaps query (e.g. length(olaps)). tapply() is the
	# standard solution, but it is slow. As per suggested in http://bit.ly/NCnYPZ
	# we use splitAsList.
	o <- sum(splitAsList(w, factor(queryHits(olaps),
																 seq_len(queryLength(olaps)))))
	return(o)
}
