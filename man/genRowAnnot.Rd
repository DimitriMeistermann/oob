% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps.R
\name{genRowAnnot}
\alias{genRowAnnot}
\title{Generate a row annotation for ComplexHeatmap}
\usage{
genRowAnnot(annot, colorScales = NULL, border = TRUE, ...)
}
\arguments{
\item{annot}{A vector of factor, character, numeric or logical. Or, a
dataframe of any of these type of value. The annotation that will be
displayed on the heatmap.}

\item{colorScales}{List or NULL. Precomputed color scales. Color scales will
be only generated for the features not described. Must be in the format of
a list named by columns of \code{annots}. Each element contains the colors at
breaks for continuous values. In the case of factors, the colors are named
to their corresponding level or in the order of the levels.}

\item{border}{Logical. Whether draw border. The value can be logical or a
string of color.}

\item{...}{Other parameters passed to \code{genColorsForAnnots}.}
}
\value{
A HeatmapAnnotation object. Can be used for example in the
\code{top_annotation} argument of \code{Heatmap}.
}
\description{
Generate a row annotation for ComplexHeatmap
}
\examples{
data("bulkLogCounts")
data("sampleAnnot")
data("DEgenesPrime_Naive")

library(ComplexHeatmap)

bestDE <- rownames(DEgenesPrime_Naive)[whichTop(DEgenesPrime_Naive$pvalue,
                                          decreasing = FALSE,
                                          top = 50)]
Heatmap(rowScale(bulkLogCounts[bestDE, ]) |> t(),
        right_annotation  = genRowAnnot(
          sampleAnnot[, c("culture_media", "line")])
        )
}
