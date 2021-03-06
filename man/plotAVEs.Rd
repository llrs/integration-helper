% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparing-models.R
\name{plotAVEs}
\alias{plotAVEs}
\title{Distribution of inner AVE}
\usage{
plotAVEs(model, loo)
}
\arguments{
\item{model}{Original sgcca object}

\item{loo}{A list with using the same model but using a leave-one-out
strategy.}
}
\value{
A histogram with a vertical line indicating the position of the
original model
}
\description{
Distribution of inner AVE
}
