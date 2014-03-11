\name{readMultiplePestPG}
\alias{readMultiplePestPG}
\title{Read multiple PestPG files with the same window coordinates, and combine
them into a single metadata \code{DataFrame} on a \code{GRanges} object.}
\usage{
readMultiplePestPG(file_list)
}
\arguments{
  \item{file_list}{a list of pestPG files, with names
  corresponding to sample.}
}
\description{
\code{readMultiplePestPG} loads multiple PestPG files (all
aligned to the same reference and using the same window
parameters in ANGSD) into a single \code{GRanges} object by
concatenating additional columns onto a single
\code{GRanges} object as metadata.
}

