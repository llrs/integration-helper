% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{dist2d}
\alias{dist2d}
\title{Calculate the distance between a line and a point}
\usage{
dist2d(p, b = c(0, 0), d = c(1, 1))
}
\arguments{
\item{p}{Point c(x, y)}

\item{b, d}{Points c(x, y) defining the line to calculate the distance with.}
}
\value{
The units of distance between the point and the line
}
\description{
The line is defined by points b and d.
}
\note{
Change the d point to change the direction of the diagonal
}
