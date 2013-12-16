semCors <- function(object,include,vertical=TRUE,titles=FALSE,layout,maximum,...){
  if (!"semPlotModel"%in%class(object)) object <- semPlotModel(object) 
  
  if (!object@Computed) stop("SEM model has not been evaluated; there are no implied covariances")
  
  if (missing(layout)) layout <- NULL
  
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
      Res[[g]]$Observed <- qgraph(round(cov2cor(object@ObsCovs[[g]]),5),maximum=ifelse(missing(maximum),1,maximum),layout=layout,...)
      layout <- Res[[g]]$Observed$layout
      if (titles) 
      {
        if (Ng > 1)
        {
          text(mean(par('usr')[1:2]),par("usr")[4],paste("Group",Groups[g],"(observed)"), adj = c(0.5,1))
        } else {
          text(mean(par('usr')[1:2]),par("usr")[4],"Observed", adj = c(0.5,1))
        }
      }
    }
    
    if (any(grepl("exp",include,ignore.case=TRUE)) | any(grepl("imp",include,ignore.case=TRUE)))
    {
      Res[[g]]$Implied <- qgraph(round(cov2cor(object@ImpCovs[[g]]),5),maximum=ifelse(missing(maximum),1,maximum),layout=layout,...)
      layout <- Res[[g]]$Implied$layout
      if (titles) 
      {
        if (Ng > 1)
        {
          text(mean(par('usr')[1:2]),par("usr")[4],paste("Group",Groups[g],"(implied)"), adj = c(0.5,1))
        } else {
          text(mean(par('usr')[1:2]),par("usr")[4],"Implied", adj = c(0.5,1))
        }
      }
    }
    
    
    if (any(grepl("dif",include,ignore.case=TRUE)) | any(grepl("res",include,ignore.case=TRUE)))
    {
      Res[[g]]$Difference <- qgraph(round(cov2cor(object@ObsCovs[[g]]) - cov2cor(object@ImpCovs[[g]]),5),maximum=ifelse(missing(maximum),.1,maximum),layout=layout,diag = TRUE, ...)
      if (titles) 
      {
        if (Ng > 1)
        {
          text(mean(par('usr')[1:2]),par("usr")[4],paste("Group",Groups[g],"(observed - implied)"), adj = c(0.5,1))
        } else {
          text(mean(par('usr')[1:2]),par("usr")[4],"Observed - Implied", adj = c(0.5,1))
        }
      }
    }
    
    
  }
  
  invisible(Res)
}