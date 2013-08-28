# semPaths.sem <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 
# semPaths.msem <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 
# semPaths.msemObjectiveML <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 


### SINGLE GROUP MODEL ###
semPlotModel.sem <- function(object, ...)
{
  
  # Check if object is of class "sem":
  if (!any(class(object)%in%c("sem","semmod"))) stop("Input must be a 'sem' object")
  
  # Define Pars:
  Pars <- data.frame(
    label = rownames(object$ram), 
    lhs = object$ram[,3],
    edge = "--",
    rhs = object$ram[,2],
    est = object$ram[,5],
    std = standardizedCoefficients(object)[,2],
    group = 1,
    fixed = object$ram[,4]==0,
    par = object$ram[,4],
    stringsAsFactors=FALSE)
  
  # Extract parameter estimates:
  Pars$est[object$ram[,4]!=0] <- object$coef[object$ram[,4]]
  
  # Fix labels:
  for (i in unique(object$ram[,4][object$ram[,4]!=0]))
  {
    if (any(Pars$label[object$ram[,4]==i]=="") & any(Pars$label[object$ram[,4]==i]!="")) 
    {
      Pars$label[object$ram[,4]==i & Pars$label==""] <- Pars$label[object$ram[,4]==i & Pars$label!=""] 
    }
  }
  
  # Name variables:
  Pars$lhs <- object$var.names[Pars$lhs]
  Pars$rhs <- object$var.names[Pars$rhs]
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = object$var.names,
    manifest = object$var.names %in% colnames(object$S),
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  # Define operators:
  Pars$edge[object$ram[,1]==2] <- "<->"
  Pars$edge[object$ram[,1]==1] <- "~>"
#   Pars$op[object$ram[,1]==1 & !Vars$manifest[match(Pars$lhs,Vars$name)] & Vars$manifest[match(Pars$rhs,Vars$name)]] <- "->"
  
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list(object$S)
  semModel@ImpCovs <- list(object$C)
  
  return(semModel)
}




### MUTLI GROUP MODEL ###
semPlotModel.msem <- semPlotModel.msemObjectiveML <- function(object, ...)
{
  
  nGroup <- length(object$ram)
  GroupNames <- object$groups
  
  ParsS <- list()
  stdobject <- standcoefmsem(object)
  
  for (g in 1:nGroup)
  {
    # Define Pars:
    Pars <- data.frame(
      label = rownames(object$ram[[g]]), 
      lhs = object$ram[[g]][,3],
      edge = "",
      rhs = object$ram[[g]][,2],
      est = object$ram[[g]][,5],
      std = stdobject[[g]][,2],
      group = GroupNames[g],
      fixed = object$ram[[g]][,4]==0,
      par = object$ram[[g]][,4],
      stringsAsFactors=FALSE)
    
    # Extract parameter estimates:
    Pars$est[object$ram[[g]][,4]!=0] <- object$coef[object$ram[[g]][,4]]
    
    # Fix labels:
    for (i in unique(object$ram[[g]][,4][object$ram[[g]][,4]!=0]))
    {
      if (any(Pars$label[object$ram[[g]][,4]==i]=="") & any(Pars$label[object$ram[[g]][,4]==i]!="")) 
      {
        Pars$label[object$ram[[g]][,4]==i & Pars$label==""] <- Pars$label[object$ram[[g]][,4]==i & Pars$label!=""] 
      }
    }
    
    # Name variables:
    Pars$lhs <- object$var.names[[g]][Pars$lhs]
    Pars$rhs <- object$var.names[[g]][Pars$rhs]
    
    
    # Define operators:
    Pars$edge[object$ram[[g]][,1]==2] <- "<->"
    Pars$edge[object$ram[[g]][,1]==1] <- "->"
    
    ParsS[[g]] <- Pars
  }
  
  Pars <- do.call("rbind",ParsS)
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = unique(unlist(object$var.names)),
    manifest = unique(unlist(object$var.names)) %in% unique(c(sapply(object$S,colnames))),
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  
  #   Pars$op[object$ram[,1]==1 & !Vars$manifest[match(Pars$lhs,Vars$name)] & Vars$manifest[match(Pars$rhs,Vars$name)]] <- "->"
  
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- object$S
  semModel@ImpCovs <- object$C
  
  return(semModel)
}
          