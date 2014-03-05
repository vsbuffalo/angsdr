## load.R -- load ANGSD data into GRanges objects

library(GenomicRanges)



readPestPG <-
# function to load ANGSD pestPG tables, and parse regions into GRanges object.
function(file, asGRanges=TRUE, useIncluded=FALSE) {
    cols <- c(region="character", chrom="factor", window.center="integer",
              watterson="numeric", pairwise="numeric", fuli="numeric",
              fayH="numeric", L="numeric", tajD="numeric", fuliF="numeric",
              fuliD="numeric", fayH="numeric", zengE="numeric", numsites="integer")

    tmp <- read.delim(x, header=FALSE, col.names=names(cols),
                      colClasses=cols)

    region_str <- gsub("\\((\\d+),(\\d+)\\)\\((\\d+),(\\d+)\\)\\((\\d+),(\\d+)\\)",
                       "\\1;;;\\2;;;\\3;;;\\4;;;\\5;;;\\6", tmp$region, perl=TRUE)

    region_chunks <- do.call(rbind,
                             lapply(region_str,
																		function(x) as.integer(unlist(strsplit(x, ";;;", fixed=TRUE)))))
		# region start
    tmp$start <- region_chunks[, 5]
    tmp$end <- region_chunks[, 6]

		# included positions
		tmp$inc_start <- region_chunks[, 3]
		tmp$inc_start <- region_chunks[, 4]

		if (asGRanges) {
			out <- with(tmp, GRanges(chrom, IRanges(start, end)))
			elementMetadata(out) <- tmp[, names(cols)[4:14]]
		} else {
			out

    out
}
