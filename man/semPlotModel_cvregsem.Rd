\name{semPlotModel_cvregsem}
\alias{cvregsemModel}
\title{
Bridge between cv_regsem output and sempaths
}
\description{
 The package regsem (Jacobucci, 2017) is designed for a specific type of SEM called regularized structural equation modelling (RegSEM). For more information about RegSEM and the implementation in R we refer to the manual written by Jacobucci (2017).This function creates a bridge between the regsem and semplot packages, making it possible to use output from the regsem() and cv_regsem() functions to create models in sempaths.}
\usage{
semPlotModel_cvregsem(object,model)
\method{semPlotModel}{cv_regsem}(object, model, \dots)
}
\arguments{
  \item{object}{
  The regsem output
}
  \item{model}{
  The cfa output used as input for the cv_regsem function
  }
}

\value{
A 'semPlotModel' object.
}
\references{
Jacobucci, R. (2017). regsem: Regularized Structural Equation Modeling. arXiv preprint  arXiv:1703.08489. 
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
Jason Nak <jasonnak@hotmail.com>
Myrthe Veenman <myrthe.veenman@hotmail.com>
}

\seealso{
\code{\link{semPlotModel}}
  \code{\link{semPaths}}
}
\examples{
## Example of fitting and plotting a cv_regsem model in semPaths

library(psych)
library(lavaan)

# use a subset of the BFI
bfi2 <- bfi[1:250,c(1:5,18,22)]
bfi2[,1] <- reverse.code(-1,bfi2[,1])

# specify a SEM model
mod <- "
f1 =~ NA*A1+A2+A3+A4+A5+O2+N3
f1~~1*f1
"

# fit the model
fit <- cfa(mod, bfi2)
out.reg <- cv_regsem(fit, type="lasso", pars_pen=c(1:7), n.lambda=23, jump =.05)

# plot the model
semPaths(semPlotModel.cv_regsem(object = out.reg, model = fit))
}

