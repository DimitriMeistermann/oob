% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixDataframe.R
\name{matrixCoord1D_2D}
\alias{matrixCoord1D_2D}
\title{Return the row and column index (2D coordinate) from a 1D coordinate in a
matrix.}
\usage{
matrixCoord1D_2D(x, Mat)
}
\arguments{
\item{x}{An integer. 1D coordinate to be transposed in 2D.}

\item{Mat}{A matrix.}
}
\value{
A vector of 2 integer. Index of row and column.
}
\description{
Return the row and column index (2D coordinate) from a 1D coordinate in a
matrix.
}
\examples{
test<-matrix(1:25,ncol=5)
test[10]
matrixCoord1D_2D(10,test)
}
