\name{semCors}
\alias{semCors}

\title{
Visually inspect implied and observed correlations
}
\description{
This function is still in devellopment.
}
\usage{
semCors(object, include, vertical = TRUE, titles = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{semPlotModel} object
}
  \item{include}{
What to include? Can be \code{"observed"}, \code{"expected"} or a vector containing both.
}
  \item{vertical}{
Should the layout be vertical or horizontal?
}
  \item{titles}{
Logical, should titles indicating the group and observed/implied correlations be plotted?
}
  \item{\dots}{
Arguments sent to \code{\link[qgraph]{qgraph}}
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
