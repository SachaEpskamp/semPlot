### SINGLE GROUP MODEL ###
matModel <- function(Lambda,Psi,Beta,Theta,manNames,latNames,ObsCovs,ImpCovs)
{
  # If names missing, set default::
  if (missing(manNames))
  {
    if (!missing(Lambda)) 
    {
      manNames <- paste0("y[",1:nrow(Lambda),"]")
    } else if (!missing(Theta))
    {
      manNames <- paste0("y[",1:nrow(Theta),"]")
    } else manNames <- character(0)
  }
  
  if (missing(latNames))
  {
    if (!missing(Lambda)) 
    {
      latNames <- paste0("eta[",1:ncol(Lambda),"]")
    } else if (!missing(Psi))
    {
        latNames <- paste0("eta[",1:ncol(Psi),"]")
    } else if (!missing(Beta))
    {
      latNames <- paste0("eta[",1:ncol(Beta),"]")
    } else latNames <- character(0)
  }
  
  # Emtpy dummy df:
  dumdf <- data.frame(
    label = character(0), 
    lhs = character(0),
    edge = character(0),
    rhs = character(0),
    est = numeric(0),
    std = numeric(0),
    group = character(0),
    fixed = logical(0),
    par = numeric(0),
    stringsAsFactors=FALSE)
  
  # Define Lambda RAM:
  if (!missing(Lambda))
  {
    LambdaRAM <- data.frame(
      label = c(outer(1:nrow(Lambda),1:ncol(Lambda),function(x,y)paste0("lambda[",x,y,"]"))), 
      lhs = rep(latNames,each=length(manNames)),
      edge = "->",
      rhs = rep(manNames,times=length(latNames)),
      est = c(Lambda),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else LambdaRAM <- dumdf
  
  # Define Theta RAM:
  if (!missing(Theta))
  {
    if (!isSymmetric(Theta)) stop("'Theta' matrix must be symmetrical.")
    Theta[upper.tri(Theta)] <- 0
    
    ThetaRAM <- data.frame(
      label = c(outer(1:nrow(Theta),1:ncol(Theta),function(x,y)paste0("theta[",x,y,"]"))), 
      lhs = rep(manNames,each=length(manNames)),
      edge = "<->",
      rhs = rep(manNames,times=length(manNames)),
      est = c(Theta),
      std = c(Theta),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else ThetaRAM <- dumdf

  # Define Psi RAM:
  if (!missing(Psi))
  {
    if (!isSymmetric(Psi)) stop("'Psi' matrix must be symmetrical.")
    Psi[upper.tri(Psi)] <- 0
    
    PsiRAM <- data.frame(
      label = c(outer(1:nrow(Psi),1:ncol(Psi),function(x,y)paste0("psi[",x,y,"]"))), 
      lhs = rep(latNames,each=length(latNames)),
      edge = "<->",
      rhs = rep(latNames,times=length(latNames)),
      est = c(Psi),
      std = c(Psi),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else PsiRAM <- dumdf
  
  # Define Psi RAM:
  if (!missing(Beta))
  {
    BetaRAM <- data.frame(
      label = c(outer(1:nrow(Beta),1:ncol(Beta),function(x,y)paste0("beta[",x,y,"]"))), 
      lhs = rep(latNames,each=length(latNames)),
      edge = "->",
      rhs = rep(latNames,times=length(latNames)),
      est = c(Beta),
      std = c(Beta),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else BetaRAM <- dumdf
  
  # Combine RAMS:
  RAM <- rbind(LambdaRAM,ThetaRAM,PsiRAM,BetaRAM)
  
  # Remove zeroes:
  RAM <- RAM[RAM$est!=0,]
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = c(manNames,latNames),
    manifest = c(manNames,latNames)%in%manNames,
    stringsAsFactors=FALSE)
  
  semModel <- new("semPlotModel")
  semModel@RAM <- RAM
  semModel@Vars <- Vars
  semModel@Computed <- FALSE
  semModel@Original <- list()
  
  if (!missing(ObsCovs))
  {
    semModel@ObsCovs <- list(ObsCovs)
  } else {
    semModel@ObsCovs <- list()
  }

  if (!missing(ImpCovs))
  {
    semModel@ImpCovs <- list(ImpCovs)
  } else {
    semModel@ImpCovs <- list()
  }
  
  return(semModel)
}


          