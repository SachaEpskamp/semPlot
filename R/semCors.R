semCors <- function(object,include,vertical=FALSE,...){
  if (!"semPlotModel"%in%class(object)) object <- semPlotModel(object) 
  
  if (!object@Computed) stop("SEM model has not been evaluated; there are no implied covariances")
  
  Ng <- max(sapply(list(object@ObsCovs,object@ImpCovs),length))
  if (missing(include))
  {
    include <- c("observed","expected")[c(length(object@ObsCovs)==Ng,length(object@ImpCovs)==Ng)]
  }
  Groups <- unique(object@Pars$group)
  
  l <- matrix(1:(Ng*length(include)),length(include),)
  if (vertical) layout(t(l)) else layout(l)
  
  Res <- list()
  
  for (g in 1:Ng)
  {
    Res[[g]] <- list()
    
    if (any(grepl("obs",include,ignore.case=TRUE)))
    {
      Res[[g]][[1]] <- qgraph(round(cov2cor(object@ObsCovs[[g]]),5),maximum=1,...)
      title(paste("Group",Groups[g],"(observed)"),line=3)
    }

    if (any(grepl("obs",include,ignore.case=TRUE)) | any(grepl("exp",include,ignore.case=TRUE)) | any(grepl("imp",include,ignore.case=TRUE)))
    {
      Res[[g]][[1]] <- qgraph(round(cov2cor(object@ImpCovs[[g]]),5),maximum=1,...)
      title(paste("Group",Groups[g],"(implied)"),line=3)
    }
  }
}