% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.cv.missoNet.R
\name{plot.cv.missoNet}
\alias{plot.cv.missoNet}
\title{Plot the cross-validation errors produced by \code{cv.missoNet}}
\usage{
\method{plot}{cv.missoNet}(x, type = c("cv.heatmap", "cv.scatter"), detailed.axis = TRUE, ...)
}
\arguments{
\item{x}{Fitted \code{cv.missoNet} object.}

\item{type}{Type of plot, can be either \code{"cv.heatmap"} or \code{"cv.scatter"}. Default is \code{type = "cv.heatmap"}.}

\item{detailed.axis}{Logical: should the detailed axes be plotted? Default is \code{plot.axis = TRUE}.}

\item{...}{Other graphical arguments used by \code{ComplexHeatmap}.}
}
\description{
Plot the cross-validation errors produced by \code{cv.missoNet}
}
