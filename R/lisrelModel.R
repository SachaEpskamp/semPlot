semPlotModel.lisrel <- function(object,...) 
{
  Res <- do.call(lisrelModel, c(object$matrices,list(...)))
  Res@Original <- list(object)
  return(Res)
}

InvEmp <- function(x)
{
  if (any(dim(x)==0)) 
  {
    return(array(0,dim=dim(x))) 
  } else {
    res <- tryCatch(solve(x), error = function(e) FALSE, silent = TRUE)
    if (is.matrix(res)) return(res) else 
    {
      res <- tryCatch(pseudoinverse(x), error = function(e) FALSE, silent = TRUE)
      if (is.matrix(res))
      {
        warning("Psuedoinverse used for singular matrix. Standardized solution might not be proper.")
        return(res) 
      } else 
      {
        warning("Uninvertable matrix found and psuedoinverse could not be computed. Standardized solutions probably not proper.")
        return(array(0, dim=dim(x)))
      }
    }
  }
}
  
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
    for (g in seq_along(m))
    {
      # Copy parSpec to par (lisrelToR compatibility)
      if (is.empty(m[[g]][['par']]) & !is.empty(m[[g]][['parSpec']]))
      {
        m[[g]][['par']] <- m[[g]][['parSpec']]
      }
      if (is.empty(m[[g]][['fixed']]))
      {
        if (!is.empty(m[[g]][['par']]))
        {
          m[[g]][['fixed']] <- m[[g]][['par']]==0
        } else if (!is.empty(m[[g]][['parSpec']]))
        {
          m[[g]][['fixed']] <- m[[g]][['parSpec']]==0
        }
      }
      if (!is.empty(m[[g]][['stdComp']]))
      {
        m[[g]][['std']] <- m[[g]][['stdComp']]
      } 
      if (is.empty(m[[g]][['est']])) m[[g]] <- list()
    }
  }    
  
  return(m)
}

is.empty <- function(x) is.null(x) || any(dim(x)==0)

### SINGLE GROUP MODEL ###
lisrelModel <- function(LY,PS,BE,TE,TY,AL,manNamesEndo,latNamesEndo,LX,PH,GA,TD,TX,KA,manNamesExo,latNamesExo,ObsCovs,ImpCovs,setExo,modelLabels = FALSE, reduce = TRUE)
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
    if (length(LY)>0 && !is.empty(LY[[1]]$est)) 
    {
      if (!is.null(rownames(LY[[1]]$est)) && !modelLabels)
      {
        manNamesEndo <- rownames(LY[[1]]$est)
      } else manNamesEndo <- paste0("y[",seq_len(nrow(LY[[1]]$est)),"]")
    } else if (length(TE)>0 && !is.empty(TE[[1]]$est))
    {
      if (!is.null(rownames(TE[[1]]$est)) && !modelLabels)
      {
        manNamesEndo <- rownames(TE[[1]]$est)
      } else manNamesEndo <- paste0("y[",seq_len(nrow(TE[[1]]$est)),"]")
    } else if (length(TY)>0 && !is.empty(TY[[1]]$est))
    {
      manNamesEndo <- paste0("y[",seq_along(TY[[1]]$est),"]")
    } else manNamesEndo <- character(0)
  }
  
  if (missing(latNamesEndo))
  {
    if (length(LY)>0 && !is.empty(LY[[1]]$est)) 
    {
      if (!is.null(colnames(LY[[1]]$est)) && !modelLabels)
      {
        latNamesEndo <- colnames(LY[[1]]$est)
      } else latNamesEndo <- paste0("eta[",1:ncol(LY[[1]]$est),"]")
    } else if (length(PS)>0 && !is.empty(PS[[1]]$est))
    {
      if (!is.null(colnames(PS[[1]]$est)) && !modelLabels)
      {
        latNamesEndo <- colnames(PS[[1]]$est)
      } else latNamesEndo <- paste0("eta[",1:ncol(PS[[1]]$est),"]")
    } else if (length(BE)>0 && !is.empty(BE[[1]]$est))
    {
      if (!is.null(colnames(BE[[1]]$est)) && !modelLabels)
      {
        latNamesEndo <- colnames(BE[[1]]$est)
      } else latNamesEndo <- paste0("eta[",1:ncol(BE[[1]]$est),"]")
    } else if (length(AL)>0 && !is.empty(AL[[1]]$est))
    {
      latNamesEndo <- paste0("eta[",seq_along(AL[[1]]$est),"]")
    } else latNamesEndo <- character(0)
  }
  
  
  # If names missing, set default::
  if (missing(manNamesExo))
  {
    if (length(LX)>0 && !is.empty(LX[[1]]$est)) 
    {
      if (!is.null(rownames(LX[[1]]$est)) && !modelLabels)
      {
        manNamesExo <- rownames(LX[[1]]$est)
      } else manNamesExo <- paste0("x[",seq_len(nrow(LX[[1]]$est)),"]")
    } else if (length(TD)>0 && !is.empty(TD[[1]]$est))
    {
      if (!is.null(rownames(TD[[1]]$est)) && !modelLabels)
      {
        manNamesExo <- rownames(TD[[1]]$est)
      } else manNamesExo <- paste0("x[",seq_len(nrow(TD[[1]]$est)),"]")
    } else if (length(TX)>0 && !is.empty(TX[[1]]$est))
    {
      manNamesExo <- paste0("x[",seq_along(TX[[1]]$est),"]")
    } else manNamesExo <- character(0)
  }
  
  if (missing(latNamesExo))
  {
    if (length(LX)>0 && !is.empty(LX[[1]]$est)) 
    {
      if (!is.null(colnames(LX[[1]]$est)) && !modelLabels)
      {
        latNamesExo <- colnames(LX[[1]]$est)
      } else latNamesExo <- paste0("xi[",1:ncol(LX[[1]]$est),"]")
    } else if (length(PH)>0 && !is.empty(PH[[1]]$est))
    {
      if (!is.null(colnames(PH[[1]]$est)) && !modelLabels)
      {
        latNamesExo <- colnames(PH[[1]]$est)
      } else latNamesExo <- paste0("xi[",1:ncol(PH[[1]]$est),"]")
    } else if (length(GA)>0 && !is.empty(GA[[1]]$est))
    {
      if (!is.null(colnames(GA[[1]]$est)) && !modelLabels)
      {
        latNamesExo <- colnames(GA[[1]]$est)
      } else latNamesExo <- paste0("xi[",1:ncol(GA[[1]]$est),"]")
    } else  if (length(KA)>0 && !is.empty(KA[[1]]$est))
    {
      latNamesExo <- paste0("xi[",seq_along(KA[[1]]$est),"]")
    } else latNamesExo <- character(0)
  }  
  
  # Check for duplicate names:
  if (!reduce)
  {
    redFun <- function(x,y,app)
    {
      x[x%in%y] <- paste0(x[x%in%y],app)
      return(x)
    }
    latNamesEndo <- redFun(latNamesEndo,c(latNamesExo,manNamesExo,manNamesEndo),"_Len")
    latNamesExo <- redFun(latNamesExo,c(latNamesEndo,manNamesEndo,manNamesExo),"_Lex")
    manNamesEndo <- redFun(manNamesEndo,c(manNamesExo,latNamesEndo,latNamesExo),"_Men")
    manNamesExo <- redFun(manNamesExo,c(manNamesEndo,latNamesEndo,latNamesExo),"_Mex")
    
  }
  
  
  
  Len <- sapply(mats,function(x)length(get(x)))
  Len <- Len[Len>0]
  if (length(unique(Len))>1) stop("Number of groups are not equal across all given LISREL matrices.")
  Ng <- max(Len)
  
  Parss <- list()
  dumPars <- data.frame(
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
  
  if (missing(ImpCovs))
  {
    modCovs <- list()
  }
  
  for (g in 1:Ng)
  {
    # Compute model implied covariance matrix and standardized matrices:
    # M is matrix list:
    M <- list()
    
    # Exogenous:
    if (length(LX)>0 && !is.empty(LX[[g]]$est))
    {
      M$LX <- LX[[g]]$est
    } else {
      M$LX <- matrix(,0,0)
    }
    
    if (length(PH)>0 && !is.empty(PH[[g]]$est))
    {
      M$PH <- PH[[g]]$est
    } else {
      M$PH <- diag(1,ncol(M$LX),ncol(M$LX))
    }
    
    if (length(TD)>0 && !is.empty(TD[[g]]$est))
    {
      M$TD <- TD[[g]]$est
    } else {
      M$TD <- matrix(0,nrow(M$LX),nrow(M$LX))
    }
    
    # Endogenous:
    if (length(LY)>0 && !is.empty(LY[[g]]$est))
    {
      M$LY <- LY[[g]]$est
    } else {
      M$LY <- matrix(,0,0)
    }
    
    if (length(PS)>0 && !is.empty(PS[[g]]$est))
    {
      M$PS <- PS[[g]]$est
    } else {
      M$PS <- diag(1,ncol(M$LY),ncol(M$LY))
    }
    
    if (length(TE)>0 && !is.empty(TE[[g]]$est))
    {
      M$TE <- TE[[g]]$est
    } else {
      M$TE <- matrix(0,nrow(M$LY),nrow(M$LY))
    }
    
    if (length(BE)>0 && !is.empty(BE[[g]]$est))
    {
      M$BE <- BE[[g]]$est
    } else {
      M$BE <- matrix(0,ncol(M$LY),ncol(M$LY))
    }
    
    if (length(GA)>0 && !is.empty(GA[[g]]$est))
    {
      M$GA <- GA[[g]]$est
    } else {
      M$GA <- matrix(0,ncol(M$LY),ncol(M$LX))
    }
    
    ImBinv <- InvEmp(diag(1,nrow(M$BE),ncol(M$BE)) - M$BE)
    
    # Implied covariances:
    XX <- with(M, LX %*% PH %*% t(LX) + TD)
    YY <- with(M, LY %*% ( ImBinv %*% (GA %*% PH %*% t(GA) + PS) %*% t(ImBinv)) %*% t(LY) + TE)
    XY <- with(M, LX %*% PH %*% t(GA) %*% t(ImBinv) %*% t(LY))
    
    if (missing(ImpCovs))
    { 
      modCovs[[g]] <- rbind(cbind(YY,t(XY)),
                            cbind(XY,XX)) 
      rownames(modCovs[[g]]) <- colnames(modCovs[[g]]) <- c(manNamesEndo,manNamesExo)
    }
    
    ## Standardize matrices
    # Diagonal matrices:
    EE <- with(M,  ( ImBinv %*% (GA %*% PH %*% t(GA) + PS) %*% t(ImBinv)) )
    
    M$De <- diag(sqrt(diag(EE)),nrow(EE),ncol(EE))
    KK <- with(M,  ( PH ) )
    M$Dk <- diag(sqrt(diag(KK)),nrow(KK),ncol(KK))
    M$Dx <- diag(sqrt(diag(XX)),nrow(XX),ncol(XX))
    M$Dy <- diag(sqrt(diag(YY)),nrow(YY),ncol(YY))
    # Inverses
    M$Dki <- InvEmp(M$Dk)
    M$Dei <- InvEmp(M$De)
    M$Dxi <- InvEmp(M$Dx)
    M$Dyi <- InvEmp(M$Dy)
    
    ## Standardize structural part:
    Mstd <- M
    # Exo:
    Mstd$LX <- M$LX %*% M$Dk
    Mstd$PH <- M$Dki %*% M$PH %*% M$Dki
    # Endo:
    Mstd$LY <- M$LY %*% M$De
    Mstd$PS <- M$Dei %*% M$PS %*% M$Dei
    Mstd$BE <- M$Dei %*% M$BE %*% M$De
    Mstd$GA <- M$Dei %*% M$GA %*% M$Dk
    
    ## Standardize measurment part:
    Mstd$LY <- M$Dyi %*% Mstd$LY
    Mstd$LX <- M$Dxi %*% Mstd$LX
    Mstd$TE <- M$Dyi %*% Mstd$TE %*% M$Dyi
    Mstd$TD <- M$Dxi %*% Mstd$TD %*% M$Dxi
    
    # Store matrices:
    if (length(LY) > 0 && !is.empty(LY[[g]]$est) && is.empty(LY[[g]]$std)) LY[[g]]$std <- Mstd$LY
    if (length(LX) > 0 && !is.empty(LX[[g]]$est) && is.empty(LX[[g]]$std)) LX[[g]]$std <- Mstd$LX
    if (length(TE) > 0 && !is.empty(TE[[g]]$est) && is.empty(TE[[g]]$std)) TE[[g]]$std <- Mstd$TE
    if (length(TD) > 0 && !is.empty(TD[[g]]$est) && is.empty(TD[[g]]$std)) TD[[g]]$std <- Mstd$TD
    if (length(PH) > 0 && !is.empty(PH[[g]]$est) && is.empty(PH[[g]]$std)) PH[[g]]$std <- Mstd$PH
    if (length(PS) > 0 && !is.empty(PS[[g]]$est) && is.empty(PS[[g]]$std)) PS[[g]]$std <- Mstd$PS
    if (length(GA) > 0 && !is.empty(GA[[g]]$est) && is.empty(GA[[g]]$std)) GA[[g]]$std <- Mstd$GA
    if (length(BE) > 0 && !is.empty(BE[[g]]$est) && is.empty(BE[[g]]$std)) BE[[g]]$std <- Mstd$BE
    
    
    # Extract matrices:
    if (length(LY)>0) LYPars <- modMat2Pars(LY[[g]],"->","lambda",symmetric=FALSE,vec=FALSE,latNamesEndo,manNamesEndo,group=paste("Group",g),exprsup="^{(y)}") else LYPars <- dumPars
    if (length(TE)>0) TEPars <- modMat2Pars(TE[[g]],"<->","theta",symmetric=TRUE,vec=FALSE,manNamesEndo,manNamesEndo,group=paste("Group",g),exprsup="^{(epsilon)}")  else TEPars <- dumPars
    if (length(PS)>0) PSPars <- modMat2Pars(PS[[g]],"<->","psi",symmetric=TRUE,vec=FALSE,latNamesEndo,latNamesEndo,group=paste("Group",g),exprsup="")  else PSPars <- dumPars
    if (length(BE)>0) BEPars <- modMat2Pars(BE[[g]],"->","beta",symmetric=FALSE,vec=FALSE,latNamesEndo,latNamesEndo,group=paste("Group",g),exprsup="")  else BEPars <- dumPars
    if (length(LX)>0) LXPars <- modMat2Pars(LX[[g]],"->","lambda",symmetric=FALSE,vec=FALSE,latNamesExo,manNamesExo,group=paste("Group",g),exprsup="^{(x)}")  else LXPars <- dumPars
    if (length(TD)>0) TDPars <- modMat2Pars(TD[[g]],"<->","theta",symmetric=TRUE,vec=FALSE,manNamesExo,manNamesExo,group=paste("Group",g),exprsup="^{(delta)}")  else TDPars <- dumPars
    if (length(PH)>0) PHPars <- modMat2Pars(PH[[g]],"<->","phi",symmetric=TRUE,vec=FALSE,latNamesExo,latNamesExo,group=paste("Group",g),exprsup="")  else PHPars <- dumPars
    if (length(GA)>0) GAPars <- modMat2Pars(GA[[g]],"->","gamma",symmetric=FALSE,vec=FALSE,latNamesExo,latNamesEndo,group=paste("Group",g),exprsup="")  else GAPars <- dumPars
    if (length(TY)>0) TYPars <- modMat2Pars(TY[[g]],"int","tau",symmetric=FALSE,vec=TRUE,"",manNamesEndo,group=paste("Group",g),exprsup="^{(y)}")  else TYPars <- dumPars
    if (length(TX)>0) TXPars <- modMat2Pars(TX[[g]],"int","tau",symmetric=FALSE,vec=TRUE,"",manNamesExo,group=paste("Group",g),exprsup="^{(x)}")  else TXPars <- dumPars
    if (length(AL)>0) ALPars <- modMat2Pars(AL[[g]],"int","alpha",symmetric=FALSE,vec=TRUE,"",latNamesEndo,group=paste("Group",g),exprsup="")  else ALPars <- dumPars
    if (length(KA)>0) KAPars <- modMat2Pars(KA[[g]],"int","kappa",symmetric=FALSE,vec=TRUE,"",latNamesExo,group=paste("Group",g),exprsup="")  else KAPars <- dumPars
    
    # Combine ParsS:
    Parss[[g]] <- rbind(LYPars,TEPars,PSPars,BEPars,LXPars,TDPars,PHPars,GAPars,TYPars,TXPars,ALPars,KAPars)
    
    # Remove zeroes:
    Parss[[g]] <- Parss[[g]][Parss[[g]]$est!=0,]
    
  }
  
  Pars <- do.call(rbind,Parss)
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo),
    manifest = c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo)%in%c(manNamesEndo,manNamesExo),
    exogenous = rep(NA,length(c(manNamesEndo,manNamesExo,latNamesEndo,latNamesExo))),
    stringsAsFactors=FALSE)

  # Remove duplicates plus factor loadings betwen mans and lats of same name:
  Vars <- Vars[!duplicated(Vars$name),]
  Pars <- Pars[!(Pars$lhs==Pars$rhs&Pars$edge!="<->"),]
  
  if (length(unique(Pars$group)) == 1) Pars$group <- ''
  
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
  semModel@Pars <- Pars
  semModel@Vars <- Vars
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
    semModel@ImpCovs <- modCovs
  }
  
  semModel@Computed <- length(semModel@ImpCovs) > 0
  
  return(semModel)
}


