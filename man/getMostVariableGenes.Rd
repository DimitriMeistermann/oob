% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markers_diffExpr.R
\name{getMostVariableGenes}
\alias{getMostVariableGenes}
\title{Compute over dispersion values for each gene.}
\usage{
getMostVariableGenes(counts, minCount = 0.01, plot = TRUE, returnPlot = FALSE)
}
\arguments{
\item{counts}{Normalized count table with genes as rows.}

\item{minCount}{Minimum average expression to not be filtered out.}

\item{plot}{Logical. Show the overdispersion plot.}

\item{returnPlot}{Logical, if \code{plot} return it as a ggplot object instead of
printing it.}
}
\value{
A ggplot graph if \code{returnPlot}, otherwise a dataframe with the
following columns:
\itemize{
\item mu: average expression
\item var: variance
\item cv2: squared coefficient of variation. Used as a dispersion value.
\item residuals: y-distance from teh regression. Can be used as an
overdispersion value.
\item residuals2: squared residuals
\item fitted: theoretical dispersion for the gene average (y value of the curve).
}
}
\description{
Compute over dispersion values for each gene.
}
\examples{
data("bulkLogCounts")
normCount<-2^(bulkLogCounts-1)
dispData<-getMostVariableGenes(normCount,minCount=1)
}
