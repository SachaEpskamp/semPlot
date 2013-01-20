\name{lisrelModel}
\alias{lisrelModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct SEM model using LISREL matrix specification.
}
\description{
This function creates a 'semPlotModel' object using matrices of the extended LISREL model (Joreskog & Sorbom, 1996). This function has two main purposes. First, it can be used to easilly create path diagrams of arbitrary SEM models without having to run an actual analysis. And second, it is specifically designed to work with the output of the 'lisrelToR' package (using \code{do.call(lisrelModel,output$matrices)}). Using  \code{\link{semPaths}} or  \code{\link{semPlotModel}} on the file path of a LISREL output file will automatically first run \code{\link[lisrelToR]{readLisrel}} and then this function.
}
\usage{
lisrelModel(LY, PS, BE, TE, TY, AL, manNamesEndo, latNamesEndo, LX, PH, GA, TD, TX, KA, manNamesExo, latNamesExo, ObsCovs, ImpCovs, setExo, modelLabels = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LY}{
Specification of the Lambda-Y matrix. See details.
}
  \item{PS}{
Specification of the Psi matrix. See details.
}
  \item{BE}{
Specification of the Beta matrix. See details.
}
  \item{TE}{
Specification of the Theta-Epsilon matrix. See details.
}
  \item{TY}{
Specification of the Tau-Y matrix. See details.
}
  \item{AL}{
Specification of the Alpha matrix. See details.
}
  \item{manNamesEndo}{
Character vector of names for the endogenous manifests.
}
  \item{latNamesEndo}{
Character vector of names for the endogenous latents.
}
  \item{LX}{
Specification of the Lambda-X matrix. See details.
}
  \item{PH}{
Specification of the Phi matrix. See details.
}
  \item{GA}{
Specification of the Gamma matrix. See details.
}
  \item{TD}{
Specification of the Theta-Delta matrix. See details.
}
  \item{TX}{
Specification of the Tau-X matrix. See details.
}
  \item{KA}{
Kappa
}
  \item{manNamesExo}{
Character vector of names for the exogenous manifests.
}
  \item{latNamesExo}{
Character vector of names for the exogenous latents.
}
  \item{ObsCovs}{
The observed covariance matrix, or a list of such matrices for each group.
}
  \item{ImpCovs}{
The implied covariance matrix, or a list of such matrices for each group.
}
  \item{setExo}{
Logical. If TRUE the 'exogenous' variable in the Variables data frame is specified. This forces \code{\link{sempaths}} to not attempt to identify which variables are endogenous and exogenous.
}
  \item{modelLabels}{
Logical. If TRUE all labels are set to the LISREL model matrix terms, as expressions. When plotted with \code{\link{sempaths}} this requires the argument \code{as.expression=c("nodes","edges")}.
}
}
\details{
The LISREL matrices can be assigned in various ways, depending on the amount of information that should be stored in the resulting model.

First, the a single matrix can be used. The values of this matrix correspond to the parameter estimates in the 'semPlotModel'. For multiple groups, a list of such matrices can be used.

to store more information, a named list of multiple matrices of the same dimensions can be used. Included in this list can be the following (but only estimates is nessesary):
\describe{
    \item{\code{est}}{Parameter estimates} 
    \item{\code{std}}{standardized parameter estimates} 
    \item{\code{par}}{Parameter numbers. 0 indicating fixed variables and parameters with the same parameter number are constrained to be equal.} 
    \item{\code{fixed}}{Logical matrix indicating if the parameter is fixed.}
}
    
If \code{std} is missing the function tries to compute standardized solutions (not yet working for intercepts). If \code{fixed} is missing it is computed from the \code{par} matrix. For multiple groups, a list containing such lists can be used.

The number of variables is extracted from the assigned matrices. Matrices that are not assigned are assumed to be empty matrices of the appropriate dimensions. e.g., Lambda-Y is assumed to be a 0 by 0 matrix if there are no endogenous variables.
}
\value{
A 'semPlotModel' object.
}
\references{
Joreskog, K. G., & Sorbom, D. (1996). LISREL 8 user's reference guide. Scientific Software.

https://github.com/SachaEpskamp/lisrelToR
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{semPlotModel}}
 \code{\link{semCors}}
 \code{\link{semPaths}}
}

\examples{
## Example of a Full LISREL model path diagram with the same number of exgenous and endogenous variables:

# Lambda matrices:
Loadings <- rbind(diag(1,2,2),diag(1,2,2),diag(1,2,2))

# Phi and Psi matrices:
LatVar <- diag(1,2,2)

# Beta matrix:
Beta <- matrix(0,2,2)
Beta[1,2] <- 1

# Theta matrices:
ManVar <- diag(1,nrow(Loadings),nrow(Loadings))

# Gamma matrix:
Gamma <- diag(1,2,2)

# Tau matrices:
ManInts <- rep(1,6)

# Alpha and Kappa matrices:
LatInts <- rep(1,2)

# Combine model:
mod <- lisrelModel(LY=Loadings,PS=LatVar,BE=Beta,TE=ManVar,
                   LX=Loadings,PH=LatVar,GA=Gamma,TD=ManVar,
                   TY=ManInts,TX=ManInts,AL=LatInts,KA=LatInts)

# Plot path diagram:
semPaths(mod, as.expression=c("nodes","edges"), sizeMan = 3, sizeInt = 1, sizeLat = 4)

# Plot path diagram with more graphical options:
semPaths(mod, as.expression=c("nodes","edges"), sizeMan = 3, sizeInt = 1, sizeLat = 4, label.prop=0.5, curve=0.5, bg="black", groups="latents", intercepts=FALSE, borders=FALSE, label.norm="O")

}
