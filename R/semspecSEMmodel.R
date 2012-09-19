
SEMpaths.semspec <- function(object,...) 
{
  invisible(SEMpaths(SEMmodel(object),...))
}



### SINGLE GROUP MODEL ###
SEMmodel.semspec <- function(object)
{
  
  # Load 'semspec':
  if (!require("semspec")) stop('semspec is required: install.packages("semspec", repos="http://R-Forge.R-project.org")')
  
  semreprObject <- semrepr(object)
  sumObject <- summary(object)
  
  # Define RAM:
  RAM <- data.frame(
    label = "", 
    lhs = semreprObject$lhs,
    edge = "",
    rhs = semreprObject$rhs,
    est = NA,
    std = NA,
    group = ifelse(is.na(semreprObject$group),"",semreprObject$group),
    fixed = FALSE,
    par = 1:nrow(semreprObject),
    stringsAsFactors=FALSE)
  
  # Extract parameter estimates:
  RAM$est[object$ram[,4]!=0] <- object$coef[object$ram[,4]]
  
  # Switch sides in regression:
  RAM[c("lhs","rhs")][semreprObject$type=="regression",] <- RAM[c("rhs","lhs")][semreprObject$type=="regression",]
  
  # Set edges:
  RAM$edge[semreprObject$type%in%c("regression","latent")] <- "->"
  RAM$edge[semreprObject$type=="covariance"] <- "<->"
  RAM$edge[semreprObject$type=="intercept"] <- "int"
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = object$var.names,
    manifest = object$var.names %in% colnames(object$S),
    stringsAsFactors=FALSE)
  
  # Define operators:
  RAM$edge[object$ram[,1]==2] <- "<->"
  RAM$edge[object$ram[,1]==1] <- "->"
#   RAM$op[object$ram[,1]==1 & !Vars$manifest[match(RAM$lhs,Vars$name)] & Vars$manifest[match(RAM$rhs,Vars$name)]] <- "->"
  
  semModel <- new("SEMmodel")
  semModel@RAM <- RAM
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list(object$S)
  semModel@ImpCovs <- list(object$C)
  
  return(semModel)
}




### MUTLI GROUP MODEL ###
SEMmodel.msem <- SEMmodel.msemObjectiveML <- function(object)
{
  
  nGroup <- length(object$ram)
  GroupNames <- object$groups
  
  RAMS <- list()
  stdobject <- standcoefmsem(object)
  
  for (g in 1:nGroup)
  {
    # Define RAM:
    RAM <- data.frame(
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
    RAM$est[object$ram[[g]][,4]!=0] <- object$coef[object$ram[[g]][,4]]
    
    # Fix labels:
    for (i in unique(object$ram[[g]][,4][object$ram[[g]][,4]!=0]))
    {
      if (any(RAM$label[object$ram[[g]][,4]==i]=="") & any(RAM$label[object$ram[[g]][,4]==i]!="")) 
      {
        RAM$label[object$ram[[g]][,4]==i & RAM$label==""] <- RAM$label[object$ram[[g]][,4]==i & RAM$label!=""] 
      }
    }
    
    # Name variables:
    RAM$lhs <- object$var.names[[g]][RAM$lhs]
    RAM$rhs <- object$var.names[[g]][RAM$rhs]
    
    
    # Define operators:
    RAM$edge[object$ram[[g]][,1]==2] <- "<->"
    RAM$edge[object$ram[[g]][,1]==1] <- "->"
    
    RAMS[[g]] <- RAM
  }
  
  RAM <- do.call("rbind",RAMS)
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = unique(unlist(object$var.names)),
    manifest = unique(unlist(object$var.names)) %in% unique(c(sapply(object$S,colnames))),
    stringsAsFactors=FALSE)
  
  
  #   RAM$op[object$ram[,1]==1 & !Vars$manifest[match(RAM$lhs,Vars$name)] & Vars$manifest[match(RAM$rhs,Vars$name)]] <- "->"
  
  semModel <- new("SEMmodel")
  semModel@RAM <- RAM
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- object$S
  semModel@ImpCovs <- object$S
  
  return(semModel)
}
          