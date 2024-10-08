% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{hierarchicalClustering}
\alias{hierarchicalClustering}
\title{Perform a hierarchical clustering from a matrix/df of observation × features}
\usage{
hierarchicalClustering(
  x,
  transpose = TRUE,
  method.dist = "euclidean",
  method.hclust = "ward.D2",
  bootstrap = FALSE,
  nboot = 10,
  PCAfirst = FALSE,
  nDimPCA = 10
)
}
\arguments{
\item{x}{A matrix or dataframe of numeric.}

\item{transpose}{Logical. If \code{transpose}, samples are columns and features
are rows.}

\item{method.dist}{A method from the "dist" function. Can be also "pearson"
or "bicor" for correlation distance.}

\item{method.hclust}{the agglomeration method to be used. This should be (an
unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
"complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC)
or "centroid" (= UPGMC).}

\item{bootstrap}{Logical. Use bootstrapping for determining the best
clustering.}

\item{nboot}{The number of bootstrap replications.}

\item{PCAfirst}{Compute a PCA before computing teh clustering, if x contains
a lot of features, it can reduce computation time.}

\item{nDimPCA}{Integer. If \code{PCAfirst}, compute a PCA first and take n first
principal components.}
}
\value{
A hclust object.
}
\description{
Perform a hierarchical clustering from a matrix/df of observation × features
}
\examples{
data(iris)
resClust <- hierarchicalClustering(iris[, c(1, 2, 3)], transpose = FALSE)
plot(resClust, hang = -1)
resClust <- hierarchicalClustering(iris[, c(1, 2, 3)], transpose = TRUE,
  bootstrap = TRUE, nboot = 20)
plot(resClust, hang = -1)
}
