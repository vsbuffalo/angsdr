## load.R -- load ANGSD data into GRanges objects

#' Read a pestPG file from ANGSD
#'
#' This function reads a pestPG file from ANGSD into either a
#' \code{\link{GRanges}} object or a plain \code{\link{data.frame}}.
#'
#' \code{readPestPG} loads a pestPG file into R, as either a \code{GRanges}
#' object (with \code{asGRanges=TRUE}) or as a \code{data.frame} (with
#' \code{asGRanges=FALSE}). When a \code{GRanges} object is returned, the
#' summary statistics will be stored as \code{elementMetadata} columns.
#' Additionally, a column named \code{numsites} specifies the number of
#' bases used by ANGSD in the calculation of the summary statistics.
#'
#' ANGSD's pestPG file includes the following estimators of theta (metadata column names in paranthesis):
#'
#' \itemize{
#'   \item Watterson's theta (\code{theta_w})
#'    \item pairwise diversity (\code{theta_pi})
#'    \item Fu and Li (\code{theta_fuli})
#'    \item Fay's H (\code{theta_fayh})
#'    \item L (\code{theta_l})
#' }
#' As well as the follwing neutrality test statistics:
#' \itemize{
#'    \item Tajima's D (\code{taj_d})
#'    \item Fu and Li F (\code{fuli_f})
#'    \item Fu and Li's D (\code{fuli_d})
#'    \item Fay's H (\code{fay_h})
#'    \item Zeng's E (\code{zeng_e})
#' }
#'
#' @param file the pestPG-formatted file from ANGSD
#' @param asGRanges \code{logical} indicating whether to return a \code{GRanges} object. If \code{FALSE}, will return a \code{data.frame}.
#' @param physical \code{logical} value indicating whether to use the window as the range, or whether to use the physical region the analysis included (for \code{GRanges} output only).
#' @param asDataFrame \code{logical} value indicating whether to use Bioconductor's \code{DataFrame} rather than a standard \code{data.frame}. Requires \code{asGRanges=FALSE}.
#'
#' @export
readPestPG <- function(file, asGRanges=TRUE, physical=FALSE, asDataFrame=FALSE) {
	# we'll use colClasses in loading, for speed.
  cols <- c(region="character", chrom="factor", window.center="integer",
            theta_w="numeric", theta_pi="numeric", theta_fuli="numeric", theta_fayh="numeric",
            theta_l="numeric", taj_d="numeric", fuli_f="numeric", fuli_d="numeric",
            fay_h="numeric", zeng_e="numeric", numsites="integer")

  angsd_raw <- read.delim(file, header=FALSE, col.names=names(cols),
                          colClasses=cols, sep="\t")

	# extract regions
  region_str <- gsub("\\((\\d+),(\\d+)\\)\\((\\d+),(\\d+)\\)\\((\\d+),(\\d+)\\)",
                     "\\1;;;\\2;;;\\3;;;\\4;;;\\5;;;\\6", angsd_raw$region, perl=TRUE)

	regcols <- c("index_start", "index_stop", "pos_start", "pos_stop", "reg_start", "reg_stop")
  region_df <- setNames(ldply(strsplit(region_str, ";;;"), .fun=as.numeric), regcols)

	# add parsed regions back in
	tmp <- cbind(angsd_raw, region_df)

  if (physical) {
		rcols <- c("pos_start", "pos_stop")
  } else {
		rcols <- c("reg_start", "reg_stop")
  }

  if (asGRanges) {
    out <- GRanges(tmp$chrom, IRanges(tmp[, rcols[1]], tmp[, rcols[2]]))
    elementMetadata(out) <- tmp[, names(cols)[4:14]] # drop first three columns
  } else {
    out <- tmp[, c("chrom", rcols, names(cols)[-c(1,2)])]
		if (asDataFrame)
			out <- as(out, "DataFrame")
  }
  out
}

#' Load another sample another ANGSD sample file from a pestPG file, and append
#' it to an existing \code{GRanges} file as element metadata.
#'
#' ANGSD pestPG data run against the same reference genome with the same window
#' size and step parameters can use the same \code{GRanges} object. Data in
#' this shape facilitates cross-population comparisons.
#' \code{appendPestPGSampleMetaData} loads a pestPG file and attaches it to an
#' existing ANGSD \code{GRanges} object. Note: this can sometimes overload memory,
#' as \code{x} is usually quite large and it must be copied. See
#' \code{readMultiplePestPG}.
#'
#' @param file the pestPG-formatted file from ANGSD.
#' @param x a \code{GRanges} object, representing ANGSD windows.
#' @param prefixes a character vector indicating the prefixes to attach to the metadata columns of \code{file} and \code{x}.
#'
appendPestPGSampleMetaData <- function(file, x, prefixes=c("file_", "x_")) {
	ff <- readPestPG(file, asGRanges=FALSE, asDataFrame=TRUE)

  # do some sanity checking
  stopifnot(nrow(ff) == length(x))
	stopifnot(all(start(x) == ff[, 2]))
	stopifnot(all(end(x) == ff[, 3]))
	stopifnot(all(seqnames(x) == ff[, 1]))

	colnames(mcols(x)) <- paste0(colnames(mcols(x)), prefixes[2])
	colnames(ff) <- paste0(colnames(ff), prefixes[1])
	mcols(x) <- cbind(mcols(x), ff[, -c(1:4)])
	x
}

#' Read multiple PestPG files with the same window coordinates, and combine
#' them into a single metadata \code{DataFrame} on a \code{GRanges} object.
#'
#' \code{readMultiplePestPG} loads multiple PestPG files (all aligned to the
#' same reference and using the same window parameters in ANGSD) into a single
#' \code{GRanges} object by concatenating additional columns onto a single
#' \code{GRanges} object as metadata.
#'
#' @param file_list a list of pestPG files, with names corresponding to sample.
#'
#' @export
readMultiplePestPG <- function(file_list) {
	if (is.null(names(file_list)))
		stop("'file_list' must be a named list, with names corresponding to samples")
	first <- TRUE
	for (i in seq_along(file_list)) {
		file <- file_list[[i]]
		sample <- names(file_list)[i]
		if (first) {
			primary <- readPestPG(file)
			colnames(mcols(primary)) <- paste0(colnames(mcols(primary)), "_", sample)
			first <- FALSE
		} else {
			d <- readPestPG(file, asGRanges=FALSE, asDataFrame=TRUE)
			colnames(d) <- paste0(colnames(d), "_", sample)
			mcols(primary) <- cbind(mcols(primary), d[, -c(1:4)])
			rm(d)
		}
	}
	primary
}
