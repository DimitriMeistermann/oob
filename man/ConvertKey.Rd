% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixDataframe.R
\name{ConvertKey}
\alias{ConvertKey}
\title{Convert a list of key from an ID to another.}
\usage{
ConvertKey(keyList, tabKey, colOldKey = 1, colNewKey = 2)
}
\arguments{
\item{keyList}{A vector of character. The list of key from an ID to be
converted (for example gene symbols).}

\item{tabKey}{A dataframe or a matrix of character. The database that
contains correspondence between each key from each ID.}

\item{colOldKey}{Integer. Column of \code{tabKey} that contains the old (same as
keyList) IDs.}

\item{colNewKey}{Integer. Column of \code{tabKey} that contains the new IDs.}
}
\value{
A vector of character.
}
\description{
Convert a list of key from an ID to another.
}
\examples{
data("bulkLogCounts")
data("humanGeneIDtable")
geneSym<-rownames(bulkLogCounts)
geneEntrez<-ConvertKey(geneSym,tabKey = humanGeneIDtable,
    colOldKey = "SYMBOL",colNewKey = "ENTREZID")
}
