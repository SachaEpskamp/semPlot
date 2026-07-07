\name{semPlotModel}
\alias{semPlotModel}
\alias{semPlotModel.default}
\alias{semPlotModel.lm}
\alias{semPlotModel.principal}
\alias{semPlotModel.psych}
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
\alias{semPlotModel.efaList}
\alias{semPlotModel.psem}
\alias{semPlotModel.wls}
\alias{semPlotModel.meta}
\alias{semPlotModel.tssem1FEM}
\alias{semPlotModel_Amos}
\alias{semPlotModel_Onyx}
\alias{semPlotModel_lavaanModel}
\alias{semPlotModel_psychonetrics}
\alias{semPlotModel_lavaan_mi}
\title{
SEM model representation
}
\description{
Methods to read a SEM object and return a \code{\link{semPlotModel-class}} object.
}
\usage{
\method{semPlotModel}{default}(object, \dots)
\method{semPlotModel}{lm}(object, \dots)
\method{semPlotModel}{principal}(object, \dots)
\method{semPlotModel}{psych}(object, cut, \dots)
\method{semPlotModel}{princomp}(object, \dots)
\method{semPlotModel}{loadings}(object, \dots)
\method{semPlotModel}{factanal}(object, \dots)
% \method{semPlotModel}{lavaan}(object)
\method{semPlotModel}{lisrel}(object, \dots)
% \method{semPlotModel}{semspec}(object)
\method{semPlotModel}{mplus.model}(object, mplusStd = c("std", "stdy", "stdyx"), \dots)
\method{semPlotModel}{sem}(object, \dots)
\method{semPlotModel}{msem}(object, \dots)
\method{semPlotModel}{msemObjectiveML}(object, \dots)
\method{semPlotModel}{efaList}(object, which, \dots)
\method{semPlotModel}{psem}(object, \dots)
\method{semPlotModel}{wls}(object, \dots)
semPlotModel_Amos(object)
semPlotModel_Onyx(object)
semPlotModel_lavaanModel(object, ...)
semPlotModel_psychonetrics(object, delta = c("variance", "ignore", "show"), ...)
semPlotModel_lavaan_mi(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object contaning the result of a SEM or GLM analysis, or a string contaning the file path to the output file of a SEM program. Or a Lavaan model.
}
\item{mplusStd}{
What standardization to use in Mplus models?
}
\item{cut}{
Loadings with absolute value below this cutoff are omitted when importing psych fa/omega objects. Defaults to 0 (keep all).
}
\item{which}{
Which solution to import from a lavaan \code{efa()} result (an \code{efaList}):
a numeric index or a name such as \code{"nf2"}. Defaults to the last (highest
number of factors).
}
\item{model}{
The original sem model (used in cvregsem)
}
\item{delta}{
How to display the Delta scaling parameters of GGM-parameterized blocks in
psychonetrics models. \code{"variance"} (default) replaces them by the implied
variances (the diagonal of \eqn{\Delta (I - \Omega)^{-1} \Delta}), which are
directly interpretable; \code{"ignore"} omits them; \code{"show"} displays the
raw scaling parameters as self-loops.
}
\item{\dots}{Arguments sent to 'lisrelModel', not used in other methods.}
}
\details{
A detailed overview of which packages are supported and what is supported for each of them will soon be on my website.

psychonetrics models (\code{semPlotModel_psychonetrics}, dispatched automatically
by \code{semPlotModel}) are supported for the \code{lvm} and \code{varcov}
frameworks, including all subtypes of their (co)variance structures:
\code{"cov"}, \code{"ggm"} (displayed as undirected network edges),
\code{"chol"} and \code{"prec"} (both displayed as implied covariances).
The development version of psychonetrics
(\code{remotes::install_github("SachaEpskamp/psychonetrics")}) is recommended.

Pooled multiple-imputation models from the lavaan.mi package
(\code{semPlotModel_lavaan_mi}, dispatched automatically) are displayed with
Rubin's-rules pooled estimates and pooled standardized solutions; observed and
implied covariance matrices are not defined for pooled fits, so
\code{\link{semCors}} is unavailable for them.

metaSEM stage-2 models (tssem2()/wls()) are supported through metaSEM's meta2semPlot(); stage-1 pooling objects give an informative error. umx models (umxRAM() returns an MxRAMModel) are supported via the OpenMx importer. psych's fa() and omega() outputs are supported (omega as a bifactor structure; pair it with semPaths(..., bifactor = 'g', layout = 'tree2')); principal() was already supported.

piecewiseSEM's psem() models are supported: directed paths from the component models (with piecewiseSEM's standardized estimates in the std column) and correlated errors (\%~~\%) as bidirectional edges.
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
