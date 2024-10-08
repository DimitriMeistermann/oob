% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/files.R
\name{write.vectorList}
\alias{write.vectorList}
\title{Write a list in a text file.}
\usage{
write.vectorList(
  list,
  filename,
  sep = "\\t",
  list.names = TRUE,
  vector.names = FALSE
)
}
\arguments{
\item{list}{The list object to save.}

\item{filename}{Path to the file to write.}

\item{sep}{Separator ("\\t"=tab-separated values).}

\item{list.names}{Save list names in the row before the values.}

\item{vector.names}{Save value name vector in the row before the values.}
}
\value{
NULL, write a file.
}
\description{
Write a list in a text file with a specific separator. If
list.names and vector.names are TRUE each element of the list is saved on a
cycle of 3 rows. First row = element name, second = vector names, third =
vector.
}
\examples{
tf <- tempfile(fileext = ".tsv")
write.vectorList(list(1:10,letters[1:10]),tf)
}
\seealso{
read.vectorList
}
