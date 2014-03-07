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
#' @param file the pestPG-formatted file from ANGSD @param asGRanges
#' \code{logical} indicating whether to return a \code{GRanges} object. If
#' \code{FALSE}, will return a \code{data.frame}. @param physical
#' \code{logical} value indicating whether to use the window as the range, or
#' whether to use the physical region the analysis included (for \code{GRanges}
#' output only).
#'
#' @export
readPestPG <- function(file, asGRanges=TRUE, physical=FALSE) {
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
  }
  out
}


