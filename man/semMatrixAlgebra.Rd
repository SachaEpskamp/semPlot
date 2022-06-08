\name{semMatrixAlgebra}
\alias{semMatrixAlgebra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract or calculate with model matrices
}
\description{
This function can be used to extract or calculate with model matrices given a \code{"semMatriModel"} object (from \code{\link{modelMatrices}}) or a \code{"semPlotModel"} object or any of the input types that can be used in \code{\link{semPlotModel}} directly.

If the model is not specified it is attempted to be identified by the given algebra.
}
\usage{
semMatrixAlgebra(object, algebra, group, simplify = TRUE, model, endoOnly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A \code{"semMatriModel"} object (from \code{\link{modelMatrices}}) or a \code{"semPlotModel"} object or any of the input types that can be used in \code{\link{semPlotModel}} directly.
}
  \item{algebra}{
An R expression to use.
}
  \item{group}{
Groups the algebra should be used on. If more than one a list is returned with the result for each group.
}
  \item{simplify}{
If TRUE and only one group is used, return output as is instead of in a list.
}
  \item{model}{
Model to be used in \code{\link{modelMatrices}}, \code{"mplus"}, \code{"ram"} or \code{"lisrel"}
}
  \item{endoOnly}{
Only needed when the model is \code{"lisrel"}, sets all variables to endogenous.
}
}
\details{
The \code{"lisrel"} model uses the following matrix names: \code{LY}, \code{TE}, \code{PS}, \code{BE}, \code{LX}, \code{TD}, \code{PH}, \code{GA}, \code{TY}, \code{TX}, \code{AL} and \code{KA}.

The \code{"mplus"} model uses the following matrix names: \code{Lambda}, \code{Nu}, \code{Theta}, \code{Kappa}, \code{Alpha}, \code{Beta}, \code{Gamma} and \code{Psi}.

The \code{"ram"} model uses the following matrix names: \code{F}, \code{A} and \code{S}.
}
\value{
A list containing output per group
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{semPlotModel}}
 \code{\link{semPlotModel-class}}
 \code{\link{modelMatrices}}
 \code{\link{lisrelModel}}
 \code{\link{ramModel}}
}

\examples{
## Mplus user guide SEM example:
outfile <- tempfile(fileext=".out")
tryres <- try({
  download.file("http://www.statmodel.com/usersguide/chap5/ex5.11.out",outfile)
})

if (!is(tryres,"try-error")){
# Plot model:
semPaths(outfile,intercepts=FALSE)

# Obtain latent regressions (mplus)
semMatrixAlgebra(outfile, Beta)

# mplus model implied covariance:
mat1 <- semMatrixAlgebra(outfile, 
  Lambda \%*\% Imin(Beta, TRUE) \%*\% Psi \%*\% t(Imin(Beta, TRUE)) \%*\% t(Lambda) + Theta)

# Lisrel model implied covariance:
mat2 <- semMatrixAlgebra(outfile, 
  LY \%*\% Imin(BE, TRUE) \%*\% PS \%*\% t(Imin(BE, TRUE)) \%*\% t(LY) + TE, endoOnly = TRUE)

# RAM model implied covariance:
mat3 <- semMatrixAlgebra(outfile, 
                 F \%*\% Imin(A,TRUE) \%*\% S \%*\% t(Imin(A, TRUE)) \%*\% t(F))

\dontrun{
# Plot:
library("qgraph")

pdf("Models.pdf",width=15,height=5)
layout(t(1:3))
qgraph(round(cov2cor(mat1),5), maximum=1, edge.labels=TRUE, layout = "spring", 
  cut = 0.4, minimum = 0.1)
title("Mplus model")
qgraph(round(cov2cor(mat2),5), maximum=1, edge.labels=TRUE, layout = "spring", 
  cut = 0.4, minimum = 0.1)
title("LISREL model")
qgraph(round(cov2cor(mat3),5), maximum=1, edge.labels=TRUE, layout = "spring", 
  cut = 0.4, minimum = 0.1)
title("RAM model")
dev.off()
}
# They are the same.
}
}
