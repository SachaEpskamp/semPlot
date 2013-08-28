\name{semPlot-tricks}
\alias{+.semPlotModel}
\alias{semPlotModel.list}
\title{
Tricks that can be used in semPlot.
}
\description{
Use a list contaning several SEM objects (from any source) to plot them as the same model. Also, the '+' operator can be used to combine two models, including in calls in \code{\link{semPaths}} and \code{\link{semPlotModel}}. See examples.
}
\usage{
\method{+}{semPlotModel}(x,y)
\method{semPlotModel}{list}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{"semPlotModel"} object
}
  \item{y}{
A \code{"semPlotModel"} object
}
  \item{object}{
An object contaning the result of a SEM or GLM analysis, or a string contaning the file path to the output file of a sEM program.
}
  \item{\dots}{
Not used.
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{semPlotModel}}
\code{\link{semPaths}}
\code{\link{semCors}}
}


\examples{
# A silly dataset:
A <- rnorm(100)
B <- A + rnorm(100)
C <- B + rnorm(100)
DF <- data.frame(A,B,C)

# Two regressions:
res1 <- lm(B ~ C, data = DF)
res2 <- lm(A ~ B + C, data = DF)

# Plot both in the same path diagram in two ways:
semPaths(res1 + res2, "model", "est", intercepts=FALSE)
semPaths(list(res1,res2), "model", "est", intercepts=FALSE)
}