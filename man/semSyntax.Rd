\name{semSyntax}
\alias{semSyntax}

\title{
Produce model syntax for various SEM software
}
\description{
This function produces a model object or model syntax for SEM software based on a \code{\link{semPlotModel-class}} object. If the input is not a \code{"semPlotModel"} object the \code{\link{semPlotModel}} function is run on the input. This allows to create model syntax for one program based on the output of another program.

Currently only the R packages 'lavaan'  (Rosseel, 2012) and 'sem' (Fox, Nie & Byrnes, 2012) are supported.
}
\usage{
semSyntax(object, syntax = "lavaan", allFixed = FALSE, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A "semPlotModel" object or any of the input possibilities for \code{\link{semPlotModel}}.
}
  \item{syntax}{
A string indicating which syntax to be used for the output. Currently supported are \code{'lavaan'} and \code{'sem'}.
}
\item{allFixed}{
Logical, should all parameters be fixed to their estimate. Useful for simulating data.
}
  \item{file}{
Path of a file the model should be written to. 
}
}
\value{
A string containing the \code{lavaan} model syntax or a \code{"semmod"} object for the \code{sem} package.
}
\references{
  Yves Rosseel (2012). lavaan: An R Package for Structural
  Equation Modeling. Journal of Statistical Software, 48(2),
  1-36. URL http://www.jstatsoft.org/v48/i02/.
  
  John Fox, Zhenghua Nie and Jarrett Byrnes (2012). sem:
  Structural Equation Models. R package version 3.0-0.
  http://CRAN.R-project.org/package=sem  
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{semPlotModel}}
\code{\link{semPlotModel-class}}
\code{\link{semPaths}}
}

\examples{
# MIMIC model, example 5.8 from mplus user guide:
Data <- read.table("http://www.statmodel.com/usersguide/chap5/ex5.8.dat")
names(Data) <- c(paste("y", 1:6, sep=""),
                 paste("x", 1:3, sep=""))

# Data <- Data[,c(7:9,1:6)]

# Model:
model.Lavaan <- 'f1 =~ y1 + y2 + y3
f2 =~ y4 + y5 + y6
f1 + f2 ~ x1 + x2 + x3 '

# Run Lavaan:
library("lavaan")
fit.Lavaan <- lavaan:::cfa(model.Lavaan, data=Data, std.lv=TRUE)

# Obtain Lavaan syntax:
model.Lavaan2 <- semSyntax(fit.Lavaan, "lavaan")

# Run Lavaan again:
fit.Lavaan2 <- lavaan:::lavaan(model.Lavaan2, data=Data)

# Compare models:
layout(t(1:2))
semPaths(fit.Lavaan,"std",title=FALSE)
title("Lavaan model 1",line=3)
semPaths(fit.Lavaan2, "std",title=FALSE)
title("Lavaan model 2",line=3)

# Convert to sem model:
model.sem <- semSyntax(fit.Lavaan, "sem")

# Run sem:
library("sem")
fit.sem <- sem:::sem(model.sem, data = Data)

# Compare models:
layout(t(1:2))
semPaths(fit.Lavaan,"std",title=FALSE)
title("Lavaan",line=3)
semPaths(fit.sem, "std",title=FALSE)
title("sem",line=3)
}
