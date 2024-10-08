% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{richUpset}
\alias{richUpset}
\title{Upset plot with additional enrichment values.}
\usage{
richUpset(featurePerGroupList, universe = NULL, pvalThreshold = 0.01)
}
\arguments{
\item{featurePerGroupList}{A list of sets (vector of character)}

\item{universe}{NULL or vector of Character. The entire list of features
(Universal set)}

\item{pvalThreshold}{Numeric. The threshold for the pvalue, represented by a
horizontal red bar. Default is 0.01.}
}
\value{
Plot in the current graphical device. In the pval bar graph, the red
line indicates an adjusted pval of 0.01 (-log10 = 2). U indicates the
universe size (total number of elements).
}
\description{
In addition to Upset plot, this method computes and represents
additional values useful for understanding the relationship between sets. The
main ones are a p-value for each overlap, and the effect size of the
association as the OEdev: \code{(observed - expected) / sqrt(universeSize)}
}
\examples{
lt <- list(set1 = sample(letters, 5),
                    set2 = sample(letters, 10),
                    set3 = sample(letters, 15)
                    )
lt$set4 <- unique(c(lt$set1,lt$set2))

richUpset(lt, universe = letters)

}
\seealso{
\code{\link[=richUpsetStats]{richUpsetStats()}}
}
