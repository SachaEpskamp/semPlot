\name{semCors}
\alias{semCors}

\title{
Visually inspect implied and observed correlations
}
\description{
This function is still in devellopment.
}
\usage{
semCors(object, include, vertical = TRUE, titles = FALSE, layout, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{semPlotModel} object
}
  \item{include}{
What to include? Can be \code{"observed"}, \code{"implied"} or \code{"difference"}, or a vector containing both. Defaults to showing observed and implied covariances.
}
  \item{vertical}{
Should the layout be vertical or horizontal?
}
  \item{titles}{
Logical, should titles indicating the group and observed/implied correlations be plotted?
}
\item{layout}{
An optional layout matrix send to \code{\link[qgraph]{qgraph}}.
}
\item{maximum}{
The maximum values as used in \code{\link[qgraph]{qgraph}}. Defaults to 1 for observed and implied covariances and 0.1 for difference graph. Important to note: Setting this lower than any of the covariances when comparing observed and implied correlations makes these graphs NOT interpretable.
}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
