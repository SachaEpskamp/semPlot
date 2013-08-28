### Path diagrams ###
# 
# setMethod("semPaths.S4",signature("lavaan"),function(object,...){
#   invisible(semPaths(semPlotModel(object),...))
# })
# 


## EXTRACT MODEL ###
semPlotModel_lavaanModel <- function(object, ...)
{

  # Check if parTable, otherwise run lavaanify:
  if (!is.data.frame(object) &  !is.list(object))
  {
    object <- lavaanify(object, ...)
  }

  varNames <- lavaanNames(object, type="ov")
  factNames <- lavaanNames(object, type="lv")
#   rm(Lambda)
  
  factNames <- factNames[!factNames%in%varNames]
  
  # Extract number of variables and factors
  n <- length(varNames)
  k <- length(factNames)
  
  # Extract parameter names:
  if (is.null(object$label)) object$label <- rep("",nrow(object))
  
  semModel <- new("semPlotModel")
  
  # Set estimates to 1 or ustart:
  object$est <- ifelse(is.na(object$ustart),1,object$ustart)
  
  if (is.null(object$group)) object$group <- ""
  
  # Create edges dataframe
  semModel@Pars <- data.frame(
     label = object$label,
    lhs = ifelse(object$op=="~"|object$op=="~1",object$rhs,object$lhs),
    edge = "--",
    rhs = ifelse(object$op=="~"|object$op=="~1",object$lhs,object$rhs),
    est = object$est,
    std = NA,
    group = object$group,
    fixed = object$free==0,
    par = object$free,
    stringsAsFactors=FALSE)

  semModel@Pars$edge[object$op=="~~"] <- "<->"  
  semModel@Pars$edge[object$op=="~*~"] <- "<->"  
  semModel@Pars$edge[object$op=="~"] <- "~>"
  semModel@Pars$edge[object$op=="=~"] <- "->"
  semModel@Pars$edge[object$op=="~1"] <- "int"
  semModel@Pars$edge[grepl("\\|",object$op)] <- "|"
  
  # Move thresholds to Thresholds slot:
  semModel@Thresholds <- semModel@Pars[grepl("\\|",semModel@Pars$edge),-(3:4)]
  # Remove thresholds from Pars:
#   semModel@Pars <- semModel@Pars[!grepl("\\|",semModel@Pars$edge),]
  
  # Remove weird edges:
  semModel@Pars <- semModel@Pars[!object$op%in%c(':=','<','>','==','|','<', '>'),]

  
  semModel@Vars <- data.frame(
    name = c(varNames,factNames),
    manifest = c(varNames,factNames)%in%varNames,
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  semModel@ObsCovs <- list()  
  semModel@ImpCovs <- list()
  semModel@Computed <- FALSE
  semModel@Original <- list(object)
  
  return(semModel)
}



