% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_variables}
\alias{plot_variables}
\title{Plot bullseye}
\usage{
plot_variables(variables)
}
\arguments{
\item{variables}{A data.frame with the weight for each variable}
}
\value{
A plot with the most important variables
}
\description{
Plot the variables of 16S, RNAseq and other origins. In RNAseq and 16S
plots variables that in both components are above the mean.
}
\seealso{
\code{\link[=variables]{variables()}}
}
