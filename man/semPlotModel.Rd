\name{semPlotModel}
\alias{semPlotModel}
\alias{semPlotModel.default}
\alias{semPlotModel.lm}
\alias{semPlotModel.principal}
\alias{semPlotModel.princomp}
\alias{semPlotModel.loadings}
\alias{semPlotModel.factanal}
% \alias{semPlotModel.lavaan}
\alias{semPlotModel.lisrel}
% \alias{semPlotModel.semspec}
\alias{semPlotModel.mplus.model}
\alias{semPlotModel.sem}
\alias{semPlotModel.msem}
\alias{semPlotModel.msemObjectiveML}
\alias{semPlotModel_Amos}
\alias{semPlotModel_Onyx}
\title{
SEM model representation
}
\description{
Methods to read a SEM object and return a \code{\link{semPlotModel-class}} object.
}
\usage{
\method{semPlotModel}{default}(object)
\method{semPlotModel}{lm}(object)
\method{semPlotModel}{principal}(object)
\method{semPlotModel}{princomp}(object)
\method{semPlotModel}{loadings}(object)
\method{semPlotModel}{factanal}(object)
% \method{semPlotModel}{lavaan}(object)
\method{semPlotModel}{lisrel}(object)
% \method{semPlotModel}{semspec}(object)
\method{semPlotModel}{mplus.model}(object)
\method{semPlotModel}{sem}(object)
\method{semPlotModel}{msem}(object)
\method{semPlotModel}{msemObjectiveML}(object)
semPlotModel_Amos(object)
semPlotModel_Onyx(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object contaning the result of a SEM or GLM analysis, or a string contaning the file path to the output file of a SEM program.
}
}
\details{
A detailed overview of which packages are supported and what is supported for each of them will soon be on my website.
}
\value{
A \code{"semPlotModel"} object. See \code{link{semPlotModel-class}}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{semPaths}}
\code{\link{semCors}}
\code{\link{semPlotModel-class}}
}
