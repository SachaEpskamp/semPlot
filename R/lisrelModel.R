### SINGLE GROUP MODEL ###
lisrelModel <- function(LY,Psi,BE,TE,TY,AL,manNamesEndo,latNamesEndo,LX,Phi,GA,TX,TX,KA,manNamesExo,latNamesExo,ObsCovs,ImpCovs,setExo)
{
  
  ### ENDOGENOUS MODEL ###
  # If names missing, set default::
  if (missing(manNamesEndo))
  {
    if (!missing(LY)) 
    {
      manNamesEndo <- paste0("y[",1:nrow(LY),"]")
    } else if (!missing(TE))
    {
      manNamesEndo <- paste0("y[",1:nrow(TE),"]")
    } else if (!missing(TY))
    {
      manNamesEndo <- paste0("y[",1:length(TY),"]")
    } else manNamesEndo <- character(0)
  }
  
  if (missing(latNamesEndo))
  {
    if (!missing(LY)) 
    {
      latNamesEndo <- paste0("eta[",1:ncol(LY),"]")
    } else if (!missing(Psi))
    {
      latNamesEndo <- paste0("eta[",1:ncol(Psi),"]")
    } else if (!missing(BE))
    {
      latNamesEndo <- paste0("eta[",1:ncol(BE),"]")
    } else if (!missing(AL))
    {
      latNamesEndo <- paste0("eta[",1:length(AL),"]")
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
  
  # Define LY RAM:
  if (!missing(LY))
  {
    LYRAM <- data.frame(
      label = c(outer(1:nrow(LY),1:ncol(LY),function(x,y)paste0("lambda[",x,y,"]^{(y)}"))), 
      lhs = rep(latNamesEndo,each=length(manNamesEndo)),
      edge = "->",
      rhs = rep(manNamesEndo,times=length(latNamesEndo)),
      est = c(LY),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else LYRAM <- dumdf
  
  # Define TE RAM:
  if (!missing(TE))
  {
    if (!isSymmetric(TE)) stop("'TE' matrix must be symmetrical.")
    TE[upper.tri(TE)] <- 0
    
    TERAM <- data.frame(
      label = c(outer(1:nrow(TE),1:ncol(TE),function(x,y)paste0("theta[",x,y,"]^{(y)}"))), 
      lhs = rep(manNamesEndo,each=length(manNamesEndo)),
      edge = "<->",
      rhs = rep(manNamesEndo,times=length(manNamesEndo)),
      est = c(TE),
      std = c(TE),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else TERAM <- dumdf

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
  
  # Define BE RAM:
  if (!missing(BE))
  {
    BERAM <- data.frame(
      label = c(outer(1:nrow(BE),1:ncol(BE),function(x,y)paste0("beta[",x,y,"]"))), 
      lhs = rep(latNamesEndo,each=length(latNamesEndo)),
      edge = "->",
      rhs = rep(latNamesEndo,times=length(latNamesEndo)),
      est = c(BE),
      std = c(BE),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else BERAM <- dumdf

  # Define TY RAM:
  if (!missing(TY))
  {
    TYRAM <- data.frame(
      label = paste0("tau[",1:length(TY),"]^{(y)}"), 
      lhs = "",
      edge = "int",
      rhs = manNamesEndo,
      est = TY,
      std = TY,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else TYRAM <- dumdf
  
  # Define AL RAM:
  if (!missing(AL))
  {
    ALRAM <- data.frame(
      label = paste0("alpha[",1:length(AL),"]"), 
      lhs = "",
      edge = "int",
      rhs = latNamesEndo,
      est = AL,
      std = AL,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else ALRAM <- dumdf
  
    ### ExoGENOUS MODEL ###
  # If names missing, set default::
  if (missing(manNamesExo))
  {
    if (!missing(LX)) 
    {
      manNamesExo <- paste0("x[",1:nrow(LX),"]")
    } else if (!missing(TX))
    {
      manNamesExo <- paste0("x[",1:nrow(TX),"]")
    } else if (!missing(TX))
    {
      manNamesExo <- paste0("x[",1:length(TX),"]")
    } else manNamesExo <- character(0)
  }
  
  if (missing(latNamesExo))
  {
    if (!missing(LX)) 
    {
      latNamesExo <- paste0("xi[",1:ncol(LX),"]")
    } else if (!missing(Phi))
    {
      latNamesExo <- paste0("xi[",1:ncol(Phi),"]")
    } else if (!missing(GA))
    {
      latNamesExo <- paste0("xi[",1:ncol(GA),"]")
    } else  if (!missing(KA))
    {
      latNamesExo <- paste0("xi[",1:length(KA),"]")
    } else latNamesExo <- character(0)
  }  
  
  # Define LX RAM:
  if (!missing(LX))
  {
    LXRAM <- data.frame(
      label = c(outer(1:nrow(LX),1:ncol(LX),function(x,y)paste0("lambda[",x,y,"]^{(x)}"))), 
      lhs = rep(latNamesExo,each=length(manNamesExo)),
      edge = "->",
      rhs = rep(manNamesExo,times=length(latNamesExo)),
      est = c(LX),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else LXRAM <- dumdf
  
  # Define TX RAM:
  if (!missing(TX))
  {
    if (!isSymmetric(TX)) stop("'TX' matrix must be symmetrical.")
    TX[upper.tri(TX)] <- 0
    
    TXRAM <- data.frame(
      label = c(outer(1:nrow(TX),1:ncol(TX),function(x,y)paste0("theta[",x,y,"]^{(x)}"))), 
      lhs = rep(manNamesExo,each=length(manNamesExo)),
      edge = "<->",
      rhs = rep(manNamesExo,times=length(manNamesExo)),
      est = c(TX),
      std = c(TX),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else TXRAM <- dumdf

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
  
  # Define GA RAM:
  if (!missing(GA))
  {
    GARAM <- data.frame(
      label = c(outer(1:nrow(GA),1:ncol(GA),function(x,y)paste0("beta[",x,y,"]"))), 
      lhs = rep(latNamesExo,each=length(latNamesEndo)),
      edge = "->",
      rhs = rep(latNamesEndo,times=length(latNamesExo)),
      est = c(GA),
      std = c(GA),
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else GARAM <- dumdf

  # Define TY RAM:
  if (!missing(TX))
  {
    TXRAM <- data.frame(
      label = paste0("tau[",1:length(TX),"]^{(x)}"), 
      lhs = "",
      edge = "int",
      rhs = manNamesExo,
      est = TX,
      std = TX,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else TXRAM <- dumdf
  
  # Define KA RAM:
  if (!missing(KA))
  {
    KARAM <- data.frame(
      label = paste0("kappa[",1:length(AL),"]"), 
      lhs = "",
      edge = "int",
      rhs = latNamesExo,
      est = KA,
      std = KA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
  } else KARAM <- dumdf
  
  
  #######
  
  # Combine RAMS:
  RAM <- rbind(LYRAM,TERAM,PsiRAM,BERAM,LXRAM,TXRAM,PhiRAM,GARAM,TYRAM,TXRAM,ALRAM,KARAM)
  
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
	setExo <- !(missing(TX) & missing(LX) & missing(Phi) & missing(GA))
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


          