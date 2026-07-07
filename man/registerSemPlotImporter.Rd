\name{registerSemPlotImporter}
\alias{registerSemPlotImporter}
\alias{validateSemPlotModel}
\title{Register a third-party semPlotModel importer}
\description{
\code{registerSemPlotImporter} lets other packages (or users) plug a converter
for their model class into \code{\link{semPlotModel}} without modifying
semPlot: when \code{semPlotModel} (and hence \code{\link{semPaths}}) receives
an object whose class matches a registered importer, that importer is called
instead of the built-in dispatch.

\code{validateSemPlotModel} checks the invariants that \code{\link{semPaths}}
relies on and stops with a message listing all violations; it is called
automatically by \code{semPaths} and is exported so importer authors can test
their converters.
}
\usage{
registerSemPlotImporter(class, fun)
validateSemPlotModel(object)
}
\arguments{
  \item{class}{A single class name (character).}
  \item{fun}{A function \code{function(object, ...)} returning a valid
  \code{\link{semPlotModel-class}} object.}
  \item{object}{A \code{"semPlotModel"} object to validate.}
}
\details{
The returned model must contain a \code{Pars} data frame with the columns
\code{label}, \code{lhs}, \code{edge}, \code{rhs}, \code{est}, \code{std},
\code{group}, \code{fixed} and \code{par}, and a \code{Vars} data frame with
\code{name}, \code{manifest} and \code{exogenous}. Recognized edge types are
\code{"->"} (factor loading), \code{"~>"} (regression), \code{"<->"}
((co)variance), \code{"--"} (undirected network edge), \code{"int"}
(intercept; \code{lhs} must be \code{""}) and \code{"|"} (threshold).
Registered importers take precedence over the built-in dispatch.
}
\value{
\code{registerSemPlotImporter} returns \code{NULL} invisibly;
\code{validateSemPlotModel} returns the object invisibly or stops.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\examples{
\dontrun{
# In a package providing models of class "myModel":
registerSemPlotImporter("myModel", function(object, ...){
  # ... build and return a semPlotModel object ...
})
# Afterwards this works directly:
# semPaths(myModelObject)
}
}
