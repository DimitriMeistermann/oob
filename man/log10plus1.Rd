% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{log10plus1}
\alias{log10plus1}
\title{log10(x+1) continuous scale for ggplot2}
\usage{
log10plus1()
}
\value{
A transformation object.
}
\description{
log10(x+1) continuous scale for ggplot2
}
\examples{
ggplot(data.frame(x=c(0,10,100,1000),y=seq_len(4)),
    mapping = aes(x=x,y=y))+geom_point()+
    scale_x_continuous(trans = log10plus1())
}
