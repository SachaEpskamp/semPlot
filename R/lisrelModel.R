### SINGLE GROUP MODEL ###
lisrelModel <- function(LambdaY,Psi,Beta,ThetaY,TauY,Alpha,manNamesEndo,latNamesEndo,LambdaX,Phi,Gamma,ThetaX,TauX,Kappa,manNamesExo,latNamesExo,ObsCovs,ImpCovs,setExo)
{
  
  ### ENDOGENOUS MODEL ###
  # If names missing, set default::
  if (missing(manNamesEndo))
  {
    if (!missing(LambdaY)) 
    {
      manNamesEndo <- paste0("y[",1:nrow(LambdaY),"]")
    } else if (!missing(ThetaY))
    {
      manNamesEndo <- paste0("y[",1:nrow(ThetaY),"]")
    } else if (!missing(TauY))
    {
      manNamesEndo <- paste0("y[",1:length(TauY),"]")
    } else manNamesEndo <- character(0)
  }
  
  if (missing(latNamesEndo))
  {
    if (!missing(LambdaY)) 
    {
      latNamesEndo <- paste0("eta[",1:ncol(LambdaY),"]")
    } else if (!missing(Psi))
    {
      latNamesEndo <- paste0("eta[",1:ncol(Psi),"]")
    } else if (!missing(Beta))
    {
      latNamesEndo <- paste0("eta[",1:ncol(Beta),"]")
    } else if (!missing(Alpha))
    {
      latNamesEndo <- paste0("eta[",1:length(Alpha),"]")
    } else latNamesEndo <- character(0)
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
  
  # Define LambdaY RAM:
  if (!missing(LambdaY))
  {
    LambdaYRAM <- data.frame(
      label = c(outer(1:nrow(LambdaY),1:ncol(LambdaY),function(x,y)paste0("lambda[",x,y,"]^{(y)}"))), 
      lhs = rep(latNamesEndo,each=length(manNamesEndo)),
      edge = "->",
      rhs = rep(manNamesEndo,times=length(latNamesEndo)),
      est = c(LambdaY),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else LambdaYRAM <- dumdf
  
  # Define ThetaY RAM:
  if (!missing(ThetaY))
  {
    if (!isSymmetric(ThetaY)) stop("'ThetaY' matrix must be symmetrical.")
    ThetaY[upper.tri(ThetaY)] <- 0
    
    ThetaYRAM <- data.frame(
      label = c(outer(1:nrow(ThetaY),1:ncol(ThetaY),function(x,y)paste0("theta[",x,y,"]^{(y)}"))), 
      lhs = rep(manNamesEndo,each=length(manNamesEndo)),
      edge = "<->",
      rhs = rep(manNamesEndo,times=length(manNamesEndo)),
      est = c(ThetaY),
      std = c(ThetaY),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else ThetaYRAM <- dumdf

  # Define Psi RAM:
  if (!missing(Psi))
  {
    if (!isSymmetric(Psi)) stop("'Psi' matrix must be symmetrical.")
    Psi[upper.tri(Psi)] <- 0
    
    PsiRAM <- data.frame(
      label = c(outer(1:nrow(Psi),1:ncol(Psi),function(x,y)paste0("psi[",x,y,"]"))), 
      lhs = rep(latNamesEndo,each=length(latNamesEndo)),
      edge = "<->",
      rhs = rep(latNamesEndo,times=length(latNamesEndo)),
      est = c(Psi),
      std = c(Psi),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else PsiRAM <- dumdf
  
  # Define Beta RAM:
  if (!missing(Beta))
  {
    BetaRAM <- data.frame(
      label = c(outer(1:nrow(Beta),1:ncol(Beta),function(x,y)paste0("beta[",x,y,"]"))), 
      lhs = rep(latNamesEndo,each=length(latNamesEndo)),
      edge = "->",
      rhs = rep(latNamesEndo,times=length(latNamesEndo)),
      est = c(Beta),
      std = c(Beta),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else BetaRAM <- dumdf

  # Define TauY RAM:
  if (!missing(TauY))
  {
    TauYRAM <- data.frame(
      label = paste0("tau[",1:length(TauY),"]^{(y)}"), 
      lhs = "",
      edge = "int",
      rhs = manNamesEndo,
      est = TauY,
      std = TauY,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else TauYRAM <- dumdf
  
  # Define Alpha RAM:
  if (!missing(Alpha))
  {
    AlphaRAM <- data.frame(
      label = paste0("alpha[",1:length(Alpha),"]"), 
      lhs = "",
      edge = "int",
      rhs = latNamesEndo,
      est = Alpha,
      std = Alpha,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else AlphaRAM <- dumdf
  
    ### ExoGENOUS MODEL ###
  # If names missing, set default::
  if (missing(manNamesExo))
  {
    if (!missing(LambdaX)) 
    {
      manNamesExo <- paste0("x[",1:nrow(LambdaX),"]")
    } else if (!missing(ThetaX))
    {
      manNamesExo <- paste0("x[",1:nrow(ThetaX),"]")
    } else if (!missing(TauX))
    {
      manNamesExo <- paste0("x[",1:length(TauX),"]")
    } else manNamesExo <- character(0)
  }
  
  if (missing(latNamesExo))
  {
    if (!missing(LambdaX)) 
    {
      latNamesExo <- paste0("xi[",1:ncol(LambdaX),"]")
    } else if (!missing(Phi))
    {
      latNamesExo <- paste0("xi[",1:ncol(Phi),"]")
    } else if (!missing(Gamma))
    {
      latNamesExo <- paste0("xi[",1:ncol(Gamma),"]")
    } else  if (!missing(Kappa))
    {
      latNamesExo <- paste0("xi[",1:length(Kappa),"]")
    } else latNamesExo <- character(0)
  }  
  
  # Define LambdaX RAM:
  if (!missing(LambdaX))
  {
    LambdaXRAM <- data.frame(
      label = c(outer(1:nrow(LambdaX),1:ncol(LambdaX),function(x,y)paste0("lambda[",x,y,"]^{(x)}"))), 
      lhs = rep(latNamesExo,each=length(manNamesExo)),
      edge = "->",
      rhs = rep(manNamesExo,times=length(latNamesExo)),
      est = c(LambdaX),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else LambdaXRAM <- dumdf
  
  # Define ThetaX RAM:
  if (!missing(ThetaX))
  {
    if (!isSymmetric(ThetaX)) stop("'ThetaX' matrix must be symmetrical.")
    ThetaX[upper.tri(ThetaX)] <- 0
    
    ThetaXRAM <- data.frame(
      label = c(outer(1:nrow(ThetaX),1:ncol(ThetaX),function(x,y)paste0("theta[",x,y,"]^{(x)}"))), 
      lhs = rep(manNamesExo,each=length(manNamesExo)),
      edge = "<->",
      rhs = rep(manNamesExo,times=length(manNamesExo)),
      est = c(ThetaX),
      std = c(ThetaX),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else ThetaXRAM <- dumdf

  # Define Phi RAM:
  if (!missing(Phi))
  {
    if (!isSymmetric(Phi)) stop("'Phi' matrix must be symmetrical.")
    Phi[upper.tri(Phi)] <- 0
    
    PhiRAM <- data.frame(
      label = c(outer(1:nrow(Phi),1:ncol(Phi),function(x,y)paste0("phi[",x,y,"]"))), 
      lhs = rep(latNamesExo,each=length(latNamesExo)),
      edge = "<->",
      rhs = rep(latNamesExo,times=length(latNamesExo)),
      est = c(Phi),
      std = c(Phi),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else PhiRAM <- dumdf
  
  # Define Gamma RAM:
  if (!missing(Gamma))
  {
    GammaRAM <- data.frame(
      label = c(outer(1:nrow(Gamma),1:ncol(Gamma),function(x,y)paste0("beta[",x,y,"]"))), 
      lhs = rep(latNamesExo,each=length(latNamesEndo)),
      edge = "->",
      rhs = rep(latNamesEndo,times=length(latNamesExo)),
      est = c(Gamma),
      std = c(Gamma),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else GammaRAM <- dumdf

  # Define TauY RAM:
  if (!missing(TauX))
  {
    TauXRAM <- data.frame(
      label = paste0("tau[",1:length(TauX),"]^{(x)}"), 
      lhs = "",
      edge = "int",
      rhs = manNamesExo,
      est = TauX,
      std = TauX,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else TauXRAM <- dumdf
  
  # Define Kappa RAM:
  if (!missing(Kappa))
  {
    KappaRAM <- data.frame(
      label = paste0("kappa[",1:length(Alpha),"]"), 
      lhs = "",
      edge = "int",
      rhs = latNamesExo,
      est = Kappa,
      std = Kappa,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else KappaRAM <- dumdf
  
  
  #######
  
  # Combine RAMS:
  RAM <- rbind(LambdaYRAM,ThetaYRAM,PsiRAM,BetaRAM,LambdaXRAM,ThetaXRAM,PhiRAM,GammaRAM,TauYRAM,TauXRAM,AlphaRAM,KappaRAM)
  
  # Remove zeroes:
  RAM <- RAM[RAM$est!=0,]
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo),
    manifest = c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo)%in%c(manNamesEndo,manNamesExo),
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  # Set exogenous:
  if (missing(setExo))
  { 
	setExo <- !(missing(ThetaX) & missing(LambdaX) & missing(Phi) & missing(Gamma))
  }
  
  if (setExo)
  {
    Vars$exogenous <- c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo)%in%c(manNamesExo,latNamesExo)
  }
  
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


          