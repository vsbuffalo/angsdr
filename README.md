# angsdr: load ANGSD data into R

This is a (very unofficial) minimal R package to load data from
[ANGSD](http://www.popgen.dk/angsd/index.php/Main_Page) into R. So far it just
loads ANGSD's `pestPG` file into a
[GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
`GRanges` object (or, optionally a `data.frame`). Eventually, I'll probably add
other convenience functions to work with ANGSD data.

## Installation

Clone this repository:

    $ git clone git@github.com:vsbuffalo/angsdr.git

and install with:

    $ R CMD INSTALL angsdr

## Example

Data from a pestPG file can be loaded with:

    library(angsdr)
    d <- readPestPG("inst/extdata/test.pestPG")

The function `overlapWidths` can be used to annotate a ANGSD `GRanges` object
with the proportion of a window that overlaps a feature, given by another
`GRanges` object. For example:

    library(rtracklayer)
    gtf <- import("your_species.gtf")
    coding <- gtf[gtf$source == "protein_coding"]
    reduced_coding <- reduce(coding, ignore.strand=TRUE)
    d$coding <- overlapWidths(d, reduced_coding)

## Development

Please feel free to fork and participate in this package's development! If
you'd like features added, open an [issue on
Github](https://github.com/vsbuffalo/angsdr/issues).
