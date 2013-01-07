# 
# semPaths.semspec <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 


### SINGLE GROUP MODEL ###
semPlotModel.semspec <- function(object)
{
  
  # Load 'semspec':
  if (!require("semspec")) stop('semspec is required: install.packages("semspec", repos="http://R-Forge.R-project.org")')
  
  semreprObject <- semrepr(object)
  sumObject <- summary(object)

  # Define RAM:
  RAM <- data.frame(
    label = "", 
    lhs = semreprObject$lhs,
    edge = "--",
    rhs = semreprObject$rhs,
    est = NA,
    std = NA,
    group = ifelse(is.na(semreprObject$group),"",semreprObject$group),
    fixed = FALSE,
    par = 0,
    stringsAsFactors=FALSE)
  
  # Label:
  if (!is.null(semreprObject$param)) RAM$label <- semreprObject$param
  
  
  # Fixed:
#   if (!is.null(semreprObject$free)) RAM$fixed <- !semreprObject$free
  if (length(sumObject$constraints$details$Constraint)>0)
  {
    spl <- strsplit(sumObject$constraints$details$Constraint,split=" == ")[grepl("==",sumObject$constraints$details$Constraint)]
    parNum <- sapply(spl,function(x)sum(x%in%RAM$label))
    parIt <- 1
    for (p in 1:length(spl))
    {
      if (parNum[p]==1)
      {
        RAM$fixed[RAM$label%in%spl[[p]]] <- TRUE
      } else if (parNum[p]==2)
      {
        RAM$par[RAM$label%in%spl[[p]]] <- parIt
        parIt <- parIt + 1
      } else warning("Error in computation of equality constraints.")
    }
  }
  
  if (max(RAM$par) < nrow(RAM))
  {
    RAM$par[RAM$par==0] <- max(RAM$par)+(1:sum(RAM$par==0))
  }
  
  # Extract parameter estimates:
  RAM$est[object$ram[,4]!=0] <- object$coef[object$ram[,4]]
  
  # Switch sides in regression:
  RAM[c("lhs","rhs")][semreprObject$type=="regression",] <- RAM[c("rhs","lhs")][semreprObject$type=="regression",]
  
  # Set edges:
  RAM$edge[semreprObject$type=="regression"] <- "~>"
  RAM$edge[semreprObject$type=="latent"] <- "->"
  RAM$edge[semreprObject$type=="covariance"] <- "<->"
  RAM$edge[semreprObject$type=="intercept"] <- "int"
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = sumObject$variables$details$Variable,
    manifest = sumObject$variables$details$Type == "Manifest",
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  # If all are latent, make guess at which are latent:
  if (all(!Vars$manifest))
  {
    for (i in 1:nrow(Vars))
    {
      Vars$manifest[i] <- !any(semreprObject$type[semreprObject$lhs==Vars$name[i]]=="latent")
    }
  }
  
  semModel <- new("semPlotModel")
  semModel@RAM <- RAM
  semModel@Vars <- Vars
  semModel@Computed <- FALSE
  semModel@Original <- list()
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  
  return(semModel)
}

