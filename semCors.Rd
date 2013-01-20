\name{semCors}
\alias{semCors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
semCors(object, include = c("observed", "expected"), vertical = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{include}{
%%     ~~Describe \code{include} here~~
}
  \item{vertical}{
%%     ~~Describe \code{vertical} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, include = c("observed", "expected"), vertical = FALSE, 
    ...) 
{
    if (!"semPlotModel" \%in\% class(object)) 
        object <- semPlotModel(object)
    if (!object@Computed) 
        stop("SEM model has not been evaluated; there are no implied covariances")
    Ng <- length(object@ObsCovs)
    Groups <- unique(object@RAM$group)
    l <- matrix(1:(Ng * length(include)), length(include), )
    if (vertical) 
        layout(t(l))
    else layout(l)
    Res <- list()
    for (g in 1:Ng) {
        Res[[g]] <- list()
        if (any(grepl("obs", include, ignore.case = TRUE))) {
            Res[[g]][[1]] <- qgraph(round(cov2cor(object@ObsCovs[[g]]), 
                5), maximum = 1, ...)
            title(paste("Group", Groups[g], "(observed)"), line = 3)
        }
        if (any(grepl("obs", include, ignore.case = TRUE)) | 
            any(grepl("exp", include, ignore.case = TRUE)) | 
            any(grepl("imp", include, ignore.case = TRUE))) {
            Res[[g]][[1]] <- qgraph(round(cov2cor(object@ImpCovs[[g]]), 
                5), maximum = 1, ...)
            title(paste("Group", Groups[g], "(implied)"), line = 3)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
