# angsdr: load ANGSD data into R

This is a (very unofficial) minimal R package to load data from ANGSD into R.
So far it just loads ANGSD's `pestPG` file into a
[GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
`GRanges` object (or, optionally a `data.frame`). Eventually, I'll probably add
other convenience functions to work with ANGSD data.

## Example

Data from a pestPG file can be loaded with:

    library(angsdr)
		d <- readPestPG("inst/extdata/test.pestPG")

## Development

Please feel free to fork and participate in this package's development!
