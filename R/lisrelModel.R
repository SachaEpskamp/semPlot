fixMatrix <- function(m)
{
  # If not a list (matrix itself added) put matrix in list (group 1) in list:
  if (!is.list(m))
  {
    if (is.matrix(m)|is.vector(m))
    {
      m <- list(list(est=m))
    } else stop("Wrong input for matrix")
  } else if ("est"%in%names(m)) {
    # Else if list, check if it is not a list of lists
    m <- list(m)
  }
  
  # Else check if empty list:
  if (length(m)>0)
  {
    # Assume multigroup. Check if all elements are list:
    if (!all(sapply(m,is.list))) stop("Not all elements are a list")
    
    # Clean each group:
    for (g in 1:length(m))
    {
      # Copy parSpec to par (lisrelToR compatibility)
      if (is.null(m[[g]][['par']]) & !is.null(m[[g]][['parSpec']]))
      {
        m[[g]][['par']] <- m[[g]][['parSpec']]
      }
      if (is.null(m[[g]][['fixed']]))
      {
        if (!is.null(m[[g]][['par']]))
        {
          m[[g]][['fixed']] <- m[[g]][['par']]==0
        } else if (!is.null(m[[g]][['parSpec']]))
        {
          m[[g]][['fixed']] <- m[[g]][['parSpec']]==0
        }
      }
      if (!is.null(m[[g]][['stdComp']]))
      {
        m[[g]][['std']] <- m[[g]][['stdComp']]
      } 
      if (is.null(m[[g]][['est']])) m[[g]] <- list()
    }
  }    
  
  return(m)
}


### SINGLE GROUP MODEL ###
lisrelModel <- function(LY,PS,BE,TE,TY,AL,manNamesEndo,latNamesEndo,LX,PH,GA,TD,TX,KA,manNamesExo,latNamesExo,ObsCovs,ImpCovs,setExo)
{
  # Input matrices either in matrix form or list containing  'est', 'std', ; fixed', and 'par' or 'parSpec' matrices. If 'stdComp' is in the list it overwrites 'std' (compatibility with 'lisrelToR' package):
  
  # Or a list of such lists for each group.
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
  
  ### NAMES ###
  # If names missing, set default::
  if (missing(manNamesEndo))
  {
    if (length(LY)>0 && !is.null(LY[[1]]$est)) 
    {
      manNamesEndo <- paste0("y[",1:nrow(LY[[1]]$est),"]")
    } else if (length(TE)>0 && !is.null(TE[[1]]$est))
    {
      manNamesEndo <- paste0("y[",1:nrow(TE[[1]]$est),"]")
    } else if (length(TY)>0 && !is.null(TY[[1]]$est))
    {
      manNamesEndo <- paste0("y[",1:length(TY[[1]]$est),"]")
    } else manNamesEndo <- character(0)
  }
  
  if (missing(latNamesEndo))
  {
    if (length(LY)>0 && !is.null(LY[[1]]$est)) 
    {
      latNamesEndo <- paste0("eta[",1:ncol(LY[[1]]$est),"]")
    } else if (length(PS)>0 && !is.null(PS[[1]]$est))
    {
      latNamesEndo <- paste0("eta[",1:ncol(PS[[1]]$est),"]")
    } else if (length(BE)>0 && !is.null(BE[[1]]$est))
    {
      latNamesEndo <- paste0("eta[",1:ncol(BE[[1]]$est),"]")
    } else if (length(AL)>0 && !is.null(AL[[1]]$est))
    {
      latNamesEndo <- paste0("eta[",1:length(AL[[1]]$est),"]")
    } else latNamesEndo <- character(0)
  }
  
  
  # If names missing, set default::
  if (missing(manNamesExo))
  {
    if (length(LX)>0 && !is.null(LX[[1]]$est)) 
    {
      manNamesExo <- paste0("x[",1:nrow(LX[[1]]$est),"]")
    } else if (length(TD)>0 && !is.null(TD[[1]]$est))
    {
      manNamesExo <- paste0("x[",1:nrow(TD[[1]]$est),"]")
    } else if (length(TX)>0 && !is.null(TX[[1]]$est))
    {
      manNamesExo <- paste0("x[",1:length(TX[[1]]$est),"]")
    } else manNamesExo <- character(0)
  }
  
  if (missing(latNamesExo))
  {
    if (length(LX)>0 && !is.null(LX[[1]]$est)) 
    {
      latNamesExo <- paste0("xi[",1:ncol(LX[[1]]$est),"]")
    } else if (length(PH)>0 && !is.null(PH[[1]]$est))
    {
      latNamesExo <- paste0("xi[",1:ncol(PH[[1]]$est),"]")
    } else if (length(GA)>0 && !is.null(GA[[1]]$est))
    {
      latNamesExo <- paste0("xi[",1:ncol(GA[[1]]$est),"]")
    } else  if (length(KA)>0 && !is.null(KA[[1]]$est))
    {
      latNamesExo <- paste0("xi[",1:length(KA[[1]]$est),"]")
    } else latNamesExo <- character(0)
  }  
  
  Len <- sapply(mats,function(x)length(get(x)))
  Len <- Len[Len>0]
  if (length(unique(Len))>1) stop("Number of groups are not equal across all given LISREL matrices.")
  Ng <- max(Len)
  
  RAMs <- list()
  dumRAM <- data.frame(
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
  for (g in 1:Ng)
  {
    # Extract matrices:
    if (length(LY)>0) LYRAM <- lisrelMat2RAM(LY[[g]],"->","lambda",symmetric=FALSE,vec=FALSE,latNamesEndo,manNamesEndo,group=paste("Group",g),exprsup="^{(y)}") else LYRAM <- dumRAM
    if (length(TE)>0) TERAM <- lisrelMat2RAM(TE[[g]],"<->","theta",symmetric=TRUE,vec=FALSE,manNamesEndo,manNamesEndo,group=paste("Group",g),exprsup="^{(epsilon)}")  else TERAM <- dumRAM
    if (length(PS)>0) PSRAM <- lisrelMat2RAM(PS[[g]],"<->","psi",symmetric=TRUE,vec=FALSE,latNamesEndo,latNamesEndo,group=paste("Group",g),exprsup="")  else PSRAM <- dumRAM
    if (length(BE)>0) BERAM <- lisrelMat2RAM(BE[[g]],"->","beta",symmetric=FALSE,vec=FALSE,latNamesEndo,latNamesEndo,group=paste("Group",g),exprsup="")  else BERAM <- dumRAM
    if (length(LX)>0) LXRAM <- lisrelMat2RAM(LX[[g]],"->","lambda",symmetric=FALSE,vec=FALSE,latNamesExo,manNamesExo,group=paste("Group",g),exprsup="^{(x)}")  else LXRAM <- dumRAM
    if (length(TD)>0) TDRAM <- lisrelMat2RAM(TD[[g]],"<->","theta",symmetric=TRUE,vec=FALSE,manNamesExo,manNamesExo,group=paste("Group",g),exprsup="^{(delta)}")  else TDRAM <- dumRAM
    if (length(PH)>0) PHRAM <- lisrelMat2RAM(PH[[g]],"<->","phi",symmetric=TRUE,vec=FALSE,latNamesExo,latNamesExo,group=paste("Group",g),exprsup="")  else PHRAM <- dumRAM
    if (length(GA)>0) GARAM <- lisrelMat2RAM(GA[[g]],"->","gamma",symmetric=FALSE,vec=FALSE,latNamesExo,latNamesEndo,group=paste("Group",g),exprsup="")  else GARAM <- dumRAM
    if (length(TY)>0) TYRAM <- lisrelMat2RAM(TY[[g]],"int","tau",symmetric=FALSE,vec=TRUE,"",manNamesEndo,group=paste("Group",g),exprsup="^{(y)}")  else TYRAM <- dumRAM
    if (length(TX)>0) TXRAM <- lisrelMat2RAM(TX[[g]],"int","tau",symmetric=FALSE,vec=TRUE,"",manNamesExo,group=paste("Group",g),exprsup="^{(x)}")  else TXRAM <- dumRAM
    if (length(AL)>0) ALRAM <- lisrelMat2RAM(AL[[g]],"int","alpha",symmetric=FALSE,vec=TRUE,"",latNamesEndo,group=paste("Group",g),exprsup="")  else ALRAM <- dumRAM
    if (length(KA)>0) KARAM <- lisrelMat2RAM(KA[[g]],"int","kappa",symmetric=FALSE,vec=TRUE,"",latNamesExo,group=paste("Group",g),exprsup="")  else KARAM <- dumRAM
    
    # Combine RAMS:
    RAMs[[g]] <- rbind(LYRAM,TERAM,PSRAM,BERAM,LXRAM,TDRAM,PHRAM,GARAM,TYRAM,TXRAM,ALRAM,KARAM)
    
    # Remove zeroes:
    RAMs[[g]] <- RAMs[[g]][RAMs[[g]]$est!=0,]
  }
  RAM <- do.call(rbind,RAMs)
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo),
    manifest = c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo)%in%c(manNamesEndo,manNamesExo),
    exogenous = rep(NA,length(c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo))),
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


