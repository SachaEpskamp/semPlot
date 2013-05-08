\name{ramModel}
\alias{ramModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct SEM model using RAM matrix specification.
}
\description{
This function creates a 'semPlotModel' object using matrices of the RAM model (McArdle & McDonald, 1984). 
}
\usage{
ramModel(A, S, F, manNames, latNames, Names, ObsCovs, ImpCovs, modelLabels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Specification of the assymmetric (A) matrix, see details.
}
  \item{S}{
Specification of the symmetric (S) matrix, see details.
}
  \item{F}{
Specification of the filter (F) matrix, see details.
}
  \item{manNames}{
Character vector of the manifest names.
}
  \item{latNames}{
Character vector of the latent names.
}
  \item{Names}{
Character vector containing all names. Defaults to \code{c(manNames,latNames)}.
}
  \item{ObsCovs}{
Observed covariancem matrix.
}
  \item{ImpCovs}{
Implied covariancem matrix.
}
  \item{modelLabels}{
Logical. If \code{TRUE} all latents are named \code{l1, l2, ...} and all manifests \code{m1, m2, ...}
}
}
\details{
The matrices can be assigned in various ways, depending on the amount of information that should be stored in the resulting model.

First, the a single matrix can be used. The values of this matrix correspond to the parameter estimates in the 'semPlotModel'. For multiple groups, a list of such matrices can be used.

to store more information, a named list of multiple matrices of the same dimensions can be used. Included in this list can be the following (but only estimates is nessesary):
\describe{
    \item{\code{est}}{Parameter estimates} 
    \item{\code{std}}{standardized parameter estimates} 
    \item{\code{par}}{Parameter numbers. 0 indicating fixed variables and parameters with the same parameter number are constrained to be equal.} 
    \item{\code{fixed}}{Logical matrix indicating if the parameter is fixed.}
}
    
If \code{std} is missing the function tries to compute standardized solutions (not yet working for intercepts). If \code{fixed} is missing it is computed from the \code{par} matrix. For multiple groups, a list containing such lists can be used.

The number of variables is extracted from the assigned matrices.
}
\value{
A 'semPlotModel' object.
}
\references{
McArdle, J. J., & McDonald, R. P. (1984). Some algebraic properties of the reticular action model for moment structures. British Journal of Mathematical and Statistical Psychology, 37(2), 234-251.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}


\seealso{
\code{\link{semPlotModel}}
 \code{\link{semCors}}
 \code{\link{semPaths}}
 \code{\link{lisrelModel}}
}

