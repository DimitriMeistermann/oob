% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimReductions.R
\name{PCA}
\alias{PCA}
\title{Principal Component Analysis}
\usage{
PCA(d, transpose = T, scale = F, center = T)
}
\arguments{
\item{d}{A matrix of numeric (in the RNA-Seq context, log counts).}

\item{transpose}{Logical. If \code{transpose}, samples are columns and features are rows.}

\item{scale}{Logical. Divide features by their standard deviation.}

\item{center}{Logical. Subtract features by their average.}
}
\value{
A list with the following element:
\itemize{
\item sdev: the standard deviations of the principal components (i.e., the square roots of the eigenvalues of the covariance/correlation matrix, though the calculation is actually done with the singular values of the data matrix).
\item rotation: the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors). The function princomp returns this in the element loadings.
\item x: if retx is true the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned. Hence, cov(x) is the diagonal matrix diag(sdev^2). For the formula method, napredict() is applied to handle the treatment of values omitted by the na.action.
\item center, scale: the centering and scaling used, or FALSE.
\item n.obs: Number of samples
\item propExplVar: proportion of explained variance by principals components on total variance.
\item transform: a list containing the scaling (sdeviations) and centering factors (means) for each principal component.
}
}
\description{
Principal Component Analysis
}
\examples{
data(iris)
pca<-PCA(iris[,seq_len(4)],transpose = FALSE,scale = TRUE,center = TRUE)
}
