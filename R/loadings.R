# semPaths.loadings <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 

### SINGLE GROUP MODEL ###
semPlotModel.loadings <- function(object, ...)
{
  
  # Check if object is of class "sem":
  if (!"loadings"%in%class(object)) stop("Input must be a 'factanal' object")
  
  
  manNames <- rownames(object)
  latNames <- colnames(object)
  
  # Define Pars:
  Pars <- data.frame(
    label = "", 
    lhs = rep(latNames,each=length(manNames)),
    edge = "--",
    rhs = rep(manNames,times=length(latNames)),
    est = c(object),
    std = c(object),
    group = "",
    fixed = FALSE,
    par = 1:length(object),
    stringsAsFactors=FALSE)
  
  
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = c(manNames[order(apply(abs(object),1,which.max))],latNames),
    manifest = c(rep(TRUE,nrow(object)),rep(FALSE,ncol(object))),
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- FALSE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  
  return(semModel)
}


          