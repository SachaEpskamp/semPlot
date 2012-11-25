fixMatrix <- function(m)
{
  if (!is.list(m))
  {
    m <- list(est=m)
  }
  if (is.null(m[['fixed']]))
  {
    if (!is.null(m[['par']]))
    {
      m[['fixed']] <- m[['par']]==0
    } else if (!is.null(m[['parSpec']]))
    {
      m[['fixed']] <- m[['parSpec']]==0
    }
  }
  if (!is.null(m[['stdComp']]))
  {
    m[['std']] <- m[['stdComp']]
  }
  if (is.null(m[['est']])) m <- list()
  return(m)
}


### SINGLE GROUP MODEL ###
lisrelModel <- function(LY,PS,BE,TE,TY,AL,manNamesEndo,latNamesEndo,LX,PH,GA,TD,TX,KA,manNamesExo,latNamesExo,ObsCovs,ImpCovs,setExo)
{
  # Input matrices either in matrix form or list containing  'est', 'std', ; fixed', and 'par' or 'parSpec' matrices. If 'stdComp' is in the list it overwrites 'std' (compatibility with 'lisrelToR' package):
  # Check input, replace matrices with list: 
  mats <- c("LY","PS","BE","TE","TY","AL","LX","PH","GA","TD","TX","KA")
  for (m in mats)
  {
    if (!do.call(missing,list(m)))
    {
      assign(m,fixMatrix(get(m)))
    } else {
      assign(m,list())
    }
  }
  
  ### ENDOGENOUS MODEL ###
  # If names missing, set default::
  if (missing(manNamesEndo))
  {
    if (length(LY)>0) 
    {
      manNamesEndo <- paste0("y[",1:nrow(LY$est),"]")
    } else if (length(TE)>0)
    {
      manNamesEndo <- paste0("y[",1:nrow(TE$est),"]")
    } else if (length(TY)>0)
    {
      manNamesEndo <- paste0("y[",1:length(TY$est),"]")
    } else manNamesEndo <- character(0)
  }
  
  if (missing(latNamesEndo))
  {
    if (length(LY)>0) 
    {
      latNamesEndo <- paste0("eta[",1:ncol(LY$est),"]")
    } else if (length(PS)>0)
    {
      latNamesEndo <- paste0("eta[",1:ncol(PS$est),"]")
    } else if (length(BE)>0)
    {
      latNamesEndo <- paste0("eta[",1:ncol(BE$est),"]")
    } else if (length(AL)>0)
    {
      latNamesEndo <- paste0("eta[",1:length(AL$est),"]")
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
  if (length(LY)>0)
  {
    LYRAM <- data.frame(
      label = c(outer(1:nrow(LY$est),1:ncol(LY$est),function(x,y)paste0("lambda[",x,y,"]^{(y)}"))), 
      lhs = rep(latNamesEndo,each=length(manNamesEndo)),
      edge = "->",
      rhs = rep(manNamesEndo,times=length(latNamesEndo)),
      est = c(LY$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(LY[['std']]))
    {
      LYRAM[['std']] <- c(LY[['std']])
    }
    if (!is.null(LY[['par']]))
    {
      LYRAM[['par']] <- c(LY[['par']])
    }
    if (!is.null(LY[['fixed']]))
    {
      LYRAM[['fixed']] <- c(LY[['fixed']])
    }
  } else LYRAM <- dumdf
  
  
  # Define TE RAM:
  if (length(TE)>0)
  {
    if (!isSymmetric(TE$est)) stop("'TE' matrix must be symmetrical.")
    TE$est[upper.tri(TE$est)] <- 0
    
    TERAM <- data.frame(
      label = c(outer(1:nrow(TE$est),1:ncol(TE$est),function(x,y)paste0("theta[",x,y,"]^{(y)}"))), 
      lhs = rep(manNamesEndo,each=length(manNamesEndo)),
      edge = "<->",
      rhs = rep(manNamesEndo,times=length(manNamesEndo)),
      est = c(TE$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(TE[['std']]))
    {
      TERAM[['std']] <- c(TE[['std']])
    }
    if (!is.null(TE[['par']]))
    {
      TERAM[['par']] <- c(TE[['par']])
    }
    if (!is.null(TE[['fixed']]))
    {
      TERAM[['fixed']] <- c(TE[['fixed']])
    }
    
  } else TERAM <- dumdf

  
  
  # Define PS RAM:
  if (length(PS)>0)
  {
    if (!isSymmetric(PS$est)) stop("'PS' matrix must be symmetrical.")
    PS$est[upper.tri(PS$est)] <- 0
    
    PSRAM <- data.frame(
      label = c(outer(1:nrow(PS$est),1:ncol(PS$est),function(x,y)paste0("psi[",x,y,"]"))), 
      lhs = rep(latNamesEndo,each=length(latNamesEndo)),
      edge = "<->",
      rhs = rep(latNamesEndo,times=length(latNamesEndo)),
      est = c(PS$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(PS[['std']]))
    {
      PSRAM[['std']] <- c(PS[['std']])
    }
    if (!is.null(PS[['par']]))
    {
      PSRAM[['par']] <- c(PS[['par']])
    }
    if (!is.null(PS[['fixed']]))
    {
      PSRAM[['fixed']] <- c(PS[['fixed']])
    }
  } else PSRAM <- dumdf
  
  # Define BE RAM:
  if (length(BE)>0)
  {
    BERAM <- data.frame(
      label = c(outer(1:nrow(BE$est),1:ncol(BE$est),function(x,y)paste0("beta[",x,y,"]"))), 
      lhs = rep(latNamesEndo,each=length(latNamesEndo)),
      edge = "->",
      rhs = rep(latNamesEndo,times=length(latNamesEndo)),
      est = c(BE$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(BE[['std']]))
    {
      BERAM[['std']] <- c(BE[['std']])
    }
    if (!is.null(BE[['par']]))
    {
      BERAM[['par']] <- c(BE[['par']])
    }
    if (!is.null(BE[['fixed']]))
    {
      BERAM[['fixed']] <- c(BE[['fixed']])
    }
  } else BERAM <- dumdf

  # Define TY RAM:
  if (length(TY)>0)
  {
    TYRAM <- data.frame(
      label = paste0("tau[",1:length(TY$est),"]^{(y)}"), 
      lhs = "",
      edge = "int",
      rhs = manNamesEndo,
      est = TY$est,
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(TY[['std']]))
    {
      TYRAM[['std']] <- c(TY[['std']])
    }
    if (!is.null(TY[['par']]))
    {
      TYRAM[['par']] <- c(TY[['par']])
    }
    if (!is.null(TY[['fixed']]))
    {
      TYRAM[['fixed']] <- c(TY[['fixed']])
    }
  } else TYRAM <- dumdf
  
  # Define AL RAM:
  if (length(AL)>0)
  {
    ALRAM <- data.frame(
      label = paste0("alpha[",1:length(AL$est),"]"), 
      lhs = "",
      edge = "int",
      rhs = latNamesEndo,
      est = AL$est,
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(AL[['std']]))
    {
      ALRAM[['std']] <- c(AL[['std']])
    }
    if (!is.null(AL[['par']]))
    {
      ALRAM[['par']] <- c(AL[['par']])
    }
    if (!is.null(AL[['fixed']]))
    {
      ALRAM[['fixed']] <- c(AL[['fixed']])
    }
  } else ALRAM <- dumdf
  
    ### ExoGENOUS MODEL ###
  # If names missing, set default::
  if (missing(manNamesExo))
  {
    if (length(LX)>0) 
    {
      manNamesExo <- paste0("x[",1:nrow(LX$est),"]")
    } else if (length(TD)>0)
    {
      manNamesExo <- paste0("x[",1:nrow(TD$est),"]")
    } else if (length(TX)>0)
    {
      manNamesExo <- paste0("x[",1:length(TX$est),"]")
    } else manNamesExo <- character(0)
  }
  
  if (missing(latNamesExo))
  {
    if (length(LX)>0) 
    {
      latNamesExo <- paste0("xi[",1:ncol(LX$est),"]")
    } else if (length(PH)>0)
    {
      latNamesExo <- paste0("xi[",1:ncol(PH$est),"]")
    } else if (length(GA)>0)
    {
      latNamesExo <- paste0("xi[",1:ncol(GA$est),"]")
    } else  if (length(KA)>0)
    {
      latNamesExo <- paste0("xi[",1:length(KA$est),"]")
    } else latNamesExo <- character(0)
  }  
  
  # Define LX RAM:
  if (length(LX)>0)
  {
    LXRAM <- data.frame(
      label = c(outer(1:nrow(LX$est),1:ncol(LX$est),function(x,y)paste0("lambda[",x,y,"]^{(x)}"))), 
      lhs = rep(latNamesExo,each=length(manNamesExo)),
      edge = "->",
      rhs = rep(manNamesExo,times=length(latNamesExo)),
      est = c(LX$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(LX[['std']]))
    {
      LXRAM[['std']] <- c(LX[['std']])
    }
    if (!is.null(LX[['par']]))
    {
      LXRAM[['par']] <- c(LX[['par']])
    }
    if (!is.null(LX[['fixed']]))
    {
      LXRAM[['fixed']] <- c(LX[['fixed']])
    }
  } else LXRAM <- dumdf
  
  # Define TD RAM:
  if (length(TD)>0)
  {
    if (!isSymmetric(TD$est)) stop("'TD' matrix must be symmetrical.")
    TD[upper.tri(TD$est)] <- 0
    
    TDRAM <- data.frame(
      label = c(outer(1:nrow(TD$est),1:ncol(TD$est),function(x,y)paste0("theta[",x,y,"]^{(x)}"))), 
      lhs = rep(manNamesExo,each=length(manNamesExo)),
      edge = "<->",
      rhs = rep(manNamesExo,times=length(manNamesExo)),
      est = c(TD$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(TD[['std']]))
    {
      TDRAM[['std']] <- c(TD[['std']])
    }
    if (!is.null(TD[['par']]))
    {
      TDRAM[['par']] <- c(TD[['par']])
    }
    if (!is.null(TD[['fixed']]))
    {
      TDRAM[['fixed']] <- c(TD[['fixed']])
    }
  } else TDRAM <- dumdf

  # Define PH RAM:
  if (length(PH)>0)
  {
    if (!isSymmetric(PH$est)) stop("'PH' matrix must be symmetrical.")
    PH[upper.tri(PH$est)] <- 0
    
    PHRAM <- data.frame(
      label = c(outer(1:nrow(PH$est),1:ncol(PH$est),function(x,y)paste0("phi[",x,y,"]"))), 
      lhs = rep(latNamesExo,each=length(latNamesExo)),
      edge = "<->",
      rhs = rep(latNamesExo,times=length(latNamesExo)),
      est = c(PH$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(PH[['std']]))
    {
      PHRAM[['std']] <- c(PH[['std']])
    }
    if (!is.null(PH[['par']]))
    {
      PHRAM[['par']] <- c(PH[['par']])
    }
    if (!is.null(PH[['fixed']]))
    {
      PHRAM[['fixed']] <- c(PH[['fixed']])
    }
  } else PHRAM <- dumdf
  
  # Define GA RAM:
  if (length(GA)>0)
  {
    GARAM <- data.frame(
      label = c(outer(1:nrow(GA$est),1:ncol(GA$est),function(x,y)paste0("beta[",x,y,"]"))), 
      lhs = rep(latNamesExo,each=length(latNamesEndo)),
      edge = "->",
      rhs = rep(latNamesEndo,times=length(latNamesExo)),
      est = c(GA$est),
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(GA[['std']]))
    {
      GARAM[['std']] <- c(GA[['std']])
    }
    if (!is.null(GA[['par']]))
    {
      GARAM[['par']] <- c(GA[['par']])
    }
    if (!is.null(GA[['fixed']]))
    {
      GARAM[['fixed']] <- c(GA[['fixed']])
    }
  } else GARAM <- dumdf

  # Define TY RAM:
  if (length(TX)>0)
  {
    TXRAM <- data.frame(
      label = paste0("tau[",1:length(TX$est),"]^{(x)}"), 
      lhs = "",
      edge = "int",
      rhs = manNamesExo,
      est = TX$est,
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(TX[['std']]))
    {
      TXRAM[['std']] <- c(TX[['std']])
    }
    if (!is.null(TX[['par']]))
    {
      TXRAM[['par']] <- c(TX[['par']])
    }
    if (!is.null(TX[['fixed']]))
    {
      TXRAM[['fixed']] <- c(TX[['fixed']])
    }
  } else TXRAM <- dumdf
  
  # Define KA RAM:
  if (length(KA)>0)
  {
    KARAM <- data.frame(
      label = paste0("kappa[",1:length(KA$est),"]"), 
      lhs = "",
      edge = "int",
      rhs = latNamesExo,
      est = KA$est,
      std = NA,
      group = "",
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!is.null(KA[['std']]))
    {
      KARAM[['std']] <- c(KA[['std']])
    }
    if (!is.null(KA[['par']]))
    {
      KARAM[['par']] <- c(KA[['par']])
    }
    if (!is.null(KA[['fixed']]))
    {
      KARAM[['fixed']] <- c(KA[['fixed']])
    }
  } else KARAM <- dumdf
  
  
  #######
  
  # Combine RAMS:
  RAM <- rbind(LYRAM,TERAM,PSRAM,BERAM,LXRAM,TDRAM,PHRAM,GARAM,TYRAM,TXRAM,ALRAM,KARAM)
  
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
	setExo <- !(length(TD)>0 & length(LX)>0 & length(PH)>0 & length(GA)>0)
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


          