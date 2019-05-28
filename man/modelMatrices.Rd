\name{modelMatrices}
\alias{modelMatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract SEM model matrices
}
\description{
Create a \code{"semMatriModel"} object. Use \code{\link{semMatrixAlgebra}} to extract or compute with these models. The structure of \code{"semMatriModel"} objects is chosen such that they can be used to create a \code{\link{semPlotModel-class}} object using \code{do.call} in combination with \code{\link{ramModel}}, \code{\link{lisrelModel}} or \code{mplusModel} (not yet implemented). See details.
}
\usage{
modelMatrices(object, model = "ram", endoOnly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{"semPlotModel"} object or any of the input types that can be used in \code{\link{semPlotModel}} directly. 
}
  \item{model}{
Model to be used, \code{"mplus"}, \code{"ram"} or \code{"lisrel"}
}
  \item{endoOnly}{
Only needed when the model is \code{"lisrel"}, sets all variables to endogenous.
}
}
\details{
The \code{"lisrel"} model uses the following matrix names: \code{LY}, \code{TE}, \code{PS}, \code{BE}, \code{LX}, \code{TD}, \code{PH}, \code{GA}, \code{TY}, \code{TX}, \code{AL} and \code{KA}. Regressions on manifest variables will cause dummy latents to be included in the model.

The \code{"mplus"} model uses the following matrix names: \code{Lambda}, \code{Nu}, \code{Theta}, \code{Kappa}, \code{Alpha}, \code{Beta}, \code{Gamma} and \code{Psi}.

The \code{"ram"} model uses the following matrix names: \code{F}, \code{A} and \code{S}.
}
\value{
a \code{"semMatriModel"} object
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{semPlotModel}}
 \code{\link{semPlotModel-class}}
 \code{\link{semMatrixAlgebra}}
 \code{\link{lisrelModel}}
 \code{\link{ramModel}}
}

\examples{
## Mplus user guide SEM example:
outfile <- tempfile(fileext=".out")
download.file("http://www.statmodel.com/usersguide/chap5/ex5.11.out",outfile)

# Plot model:
semPaths(outfile, intercepts = FALSE)

# Extract RAM:
RAM <- modelMatrices(outfile, "ram")
semPaths(do.call(ramModel, RAM), as.expression = "edges", intercepts = FALSE)

# Extract LISREL:
LISREL <- modelMatrices(outfile, "lisrel")
semPaths(do.call(lisrelModel, LISREL), as.expression = "edges", intercepts = FALSE)
}

