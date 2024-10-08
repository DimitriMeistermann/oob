% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{mostDistantColor}
\alias{mostDistantColor}
\title{Best theoretical color palette (wrapper for qualpal)}
\usage{
mostDistantColor(
  n,
  colorspace = NULL,
  cvd = c("protan", "deutan", "tritan"),
  cvd_severity = 0
)
}
\arguments{
\item{n}{The number of colors to generate.}

\item{colorspace}{A color space to generate colors from. See
?qualpalr::qualpal. If NULL computed regarding the number of asked color.}

\item{cvd}{Color vision deficiency adaptation. Use cvd_severity to set the
severity of color vision deficiency to adapt to. Permissible values are
"protan", "deutan", and "tritan".}

\item{cvd_severity}{Severity of color vision deficiency to adapt to. Can take
any value from 0, for normal vision (the default), and 1, for dichromatic
vision.}
}
\value{
Colors in hex format.
}
\description{
Best theoretical color palette (wrapper for qualpal)
}
\examples{
mostDistantColor(3)

}
