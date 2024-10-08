% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{formatNumber2Character}
\alias{formatNumber2Character}
\title{Convert numeric to string, add 0 to the number to respect lexicographical
order.}
\usage{
formatNumber2Character(x, digit = max(nchar(as.character(x))))
}
\arguments{
\item{x}{A numeric vector.}

\item{digit}{A single integer value. The maximum number of digits in the
number sequence. It will determine the number of 0 to add.}
}
\value{
A charactervector.
}
\description{
Convert numeric to string, add 0 to the number to respect lexicographical
order.
}
\examples{
formatNumber2Character(seq_len(10))
formatNumber2Character(seq_len(10),digit = 4)
}
