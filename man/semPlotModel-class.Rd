\name{semPlotModel-class}
\Rdversion{1.1}
\docType{class}
\alias{semPlotModel-class}

\title{Class \code{"semPlotModel"}}
\description{
Representation of SEM models, can be used by \code{\link{semPaths}}, \code{\link{semCors}} and \code{\link{semSyntax}.} See \code{\link{semPlotModel-edit}} for utility functions on how to edit this model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("semPlotModel", ...)}.
%%  ~~ describe objects here ~~ 
}

    
\section{Slots}{
  \describe{
    \item{\code{Pars}:}{Object of class \code{"data.frame"} indicating the parameters used in the SEM model. this must contain the following elements, in order:
    \describe{
    \item{\code{label}}{The name of the parameter, used as edge label in the graph.}
    \item{\code{lhs}}{Name of the variable on the left hand side of the path.}
    \item{\code{edge}}{String as indicator of the edge. This can be one of the following:
    \describe{
      \item{\code{->}}{Factor loading}
      \item{\code{~>}}{Regression. The same as \code{'->'} in that it results in a directed edge from the left hand side to the right hand side, but \code{'~>'} differs in that if the right hand side is manifest and the left hand side is an exogenous latent the right hand side is interpreted as an endogenous variable rather than an exogenous variable.}
      \item{\code{<->}}{(co)variance}
      \item{\code{int}}{intercept, The left hand side should be "" and the right hand side indicates the variable to which the intercept belongs.}
      \item{\code{--}}{Undirected edge. Only used as dummy encoding and in cases the parameter can not be interpreted (usually this indicates something that is not yet supported)}

    }}
    
    \item{\code{rhs}}{Name of the variable on the left hand side of the path.}
    \item{\code{est}}{Parameter estimate.}
    \item{\code{est}}{Standardized parameter estimate.}
    \item{\code{group}}{Character of the name of the group the parameter belongs to.}
    \item{\code{fixed}}{Logical indicating if the parameter is fixed.}
    \item{\code{par}}{Parameter number. 0 indicates the parameter is fixed and parameters with the same parameter number are constrained to be equal.}
    \item{\code{knot}}{Knot number. 0 indicates the edge is not knotted and edges with the same knot number are knotted together. Only used to indicate interactions in 'lm' models and can be omitted.}
    
    }
    }
    \item{\code{Vars}:}{Object of class \code{"data.frame"} indicating the variables used in the SEM model. Must have the following elements:
    \describe{
      \item{name}{Name of the variable}
      \item{manifest}{Logical indicating if the variable is manifest}
      \item{exogenous}{Logical indicating if the variable is exogenous. If \code{NA}} (the default) \code{\link{semPaths}} will attempt to detect which variables are exogenous.
    
    }}
    \item{\code{Thresholds}:}{Object of class \code{"data.frame"} indicating the thresholds in the SEM model. It is the same as \code{Pars} except it does not have the elements \code{'edge'} and \code{'rhs'}.}
    
    \item{\code{Computed}:}{Object of class \code{"logical"} indicating if the SEM model was computed or if the object only indicates a structure.}
    \item{\code{ObsCovs}:}{Object of class \code{"list"} containing observed covariance matrices for each group. If available.}
    \item{\code{ImpCovs}:}{Object of class \code{"list"} containing implied covariance matrices for each group. If available.}
    \item{\code{Original}:}{Object of class \code{"list"} containing the original object used as input (or multiple objects if the \code{'+'} operator was used to combine objects.) }
  }
}


\section{Methods}{
No methods defined with class "semPlotModel" in the signature.
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{semPlotModel}}
\code{\link{semPaths}}
\code{\link{semCors}}
\code{\link{semSyntax}}
\code{\link{semPlotModel-edit}}
}
\examples{
showClass("semPlotModel")
}
\keyword{classes}
