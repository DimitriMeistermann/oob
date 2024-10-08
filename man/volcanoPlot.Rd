% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{volcanoPlot}
\alias{volcanoPlot}
\title{General volcano plot}
\usage{
volcanoPlot(
  d = NULL,
  effectSizeCol,
  adjPvalCol,
  labelCol,
  padjThres = 0.05,
  minEffectSize = 0,
  topShownPerSide = 15,
  returnGraph = FALSE,
  neutralVal = 0,
  ...
)
}
\arguments{
\item{d}{A dataframe containing the data needed for the volcano plot. Must
have column names. Can be \code{NULL} if \code{effectSizeCol}, \code{adjPvalCol} and
\code{labelCol} are vectors.}

\item{effectSizeCol}{Column name containing the effect size column
(Log2FoldChange for example). Can also be a vector of numeric containing
the effect size values.}

\item{adjPvalCol}{Column name containing the adjusted pval column. Can also
be a vector of numeric containing the padj values.}

\item{labelCol}{Column name containing the feature name column (gene name for
example).  Can also be a vector of character containing the labels.}

\item{padjThres}{Significativity threshold of adjusted p-value for consider a
feature significant.}

\item{minEffectSize}{Absolute minimum effect size to consider a feature
significant.}

\item{topShownPerSide}{Number of feature shown at the left and right side of
the volcano plot.}

\item{returnGraph}{Logical. Return the graph as a ggplot object instead of
printing it.}

\item{neutralVal}{Value considered as null effect size.}

\item{...}{Parameters passed to geom_repel}
}
\value{
Plot in the current graphical device or a ggplot object if
\code{returnGraph=TRUE}.
}
\description{
General volcano plot
}
\examples{
data("DEgenesPrime_Naive")
data(sampleAnnot)
volcanoPlot(d = DEgenesPrime_Naive,effectSizeCol = "log2FoldChange",
    adjPvalCol = "padj", minEffectSize = 1,
    labelCol = rownames(DEgenesPrime_Naive))

}
