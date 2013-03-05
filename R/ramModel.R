

### SINGLE GROUP MODEL ###
ramModel <- function(A,S,F,manNames,latNames,Names,ObsCovs,ImpCovs,modelLabels = FALSE)
{
  # Input matrices either in matrix form or list containing  'est', 'std', ; fixed', and 'par' or 'parSpec' matrices. If 'stdComp' is in the list it overwrites 'std' (compatibility with 'lisrelToR' package):
  
  # Or a list of such lists for each group.
  # Check input, replace matrices with list: 
  mats <- c("A","S","F")
  for (m in mats)
  {
    if (!do.call(missing,list(m)))
    {
      assign(m,fixMatrix(get(m)))
    } else {
      assign(m,list())
    }
  }
  #browser()
  ### Fix matrices:
  matList <- list(A,S,F)
  Ng <- max(sapply(matList,length))
  Nvar <- max(sapply(matList,function(x)sapply(x,function(y)ncol(y$est))))
  if (length(F)>0 && !is.null(F[[1]]$est))
  {
    Nman <- max(sapply(F,function(y)nrow(y$est)))
  } else 
  {
    if (!missing(manNames)) Nman <- length(manNames) else Nman <- Nvar
  }
  
  if (!missing(manNames) & !missing(latNames))
  {
    if (Nvar!=length(c(manNames,latNames))) stop("Number of variables in model not equal to given number of names")
  }
  
  if (!missing(manNames))
  {
    if (Nman!=length(manNames)) stop("Number of manifest variables in model not equal to given number of names")
  }
  
  # Fix A:
  if (length(A)==0)
  {
    A <- lapply(seq_len(Ng),function(x)list(est=matrix(0,Nvar,Nvar)))
  } else if (length(A) < Ng) A <- rep(A,length=Ng)
  
  # Fix S
  if (length(S)==0)
  {
    S <- lapply(seq_len(Ng),function(x)list(est=matrix(0,Nvar,Nvar)))
  } else if (length(S) < Ng) S <- rep(S,length=Ng)
  
  # Fix F:
  if (length(F)==0)
  {
    F <- lapply(seq_len(Ng),function(x)list(est=cbind(diag(1,Nman,Nman),matrix(0,Nman,Nvar-Nman))))
  } else if (length(F) < Ng) F <- rep(F,length=Ng)
  
  ### NAMES ###
  # If names missing, set default::
  if (missing(manNames))
  {
    if (length(F)>0 && !is.null(F[[1]]$est)) 
    {
      if (!is.null(colnames(F[[1]]$est)) && !modelLabels)
      {
        manNames <- colnames(F[[1]]$est)[colSums(F[[1]]$est)>0]
      } else manNames <- paste0(rep("m",Nman),seq_len(Nman))
    } else manNames <- paste0(rep("m",Nman),seq_len(Nman))
  }
    
  if (missing(latNames))
  {
    if (length(F)>0 && !is.null(F[[1]]$est)) 
    {
      if (!is.null(colnames(F[[1]]$est)) && !modelLabels)
      {
        latNames <- colnames(F[[1]]$est)[colSums(F[[1]]$est)==0]
      } else latNames <- paste0(rep("l",Nvar-Nman),seq_len(Nvar-Nman))
    } else latNames <- paste0(rep("l",Nvar-Nman),seq_len(Nvar-Nman))
  }
  
  if (missing(Names))
  {
    if (length(F)>0 && !is.null(F[[1]]$est)) 
    {
      if (!is.null(colnames(F[[1]]$est)) && !modelLabels)
      {
        Names <- colnames(F[[1]]$est)
      } else Names <- c(manNames,latNames)
    } else Names <- c(manNames,latNames)
  }
  
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
    M <- list(A=A[[g]]$est, S=S[[g]]$est, F=F[[g]]$est)    
    
    IminAinv <- InvEmp(diag(1,nrow(M$A),ncol(M$A)) - M$A)
    if (missing(ImpCovs))
    { 
      modCovs[[g]] <- with(M, F %*% IminAinv %*% S %*% t(IminAinv) %*% t(F))
        
      rownames(modCovs[[g]]) <- colnames(modCovs[[g]]) <- manNames
    }
    
    Mstd <- M
    ## Standardize matrices
    I <- diag(nrow(M$S))
    expCov <- IminAinv %*% M$S %*% t(IminAinv)
    invSDs <- 1/sqrt(diag(expCov))
    diag(I) <- invSDs
    # standardize the A, S and M matrices
    # A paths are value*sd(from)/sd(to) = I %*% A %*% solve(I)
    # S paths are value/(sd(from*sd(to))) = I %*% S %*% I
    Mstd$A <- I %*% M$A %*% solve(I)
    Mstd$S <- I %*% M$S %*% I
    
    # Store matrices:
    if (length(A) > 0 && !is.null(A[[g]]$est) && is.null(A[[g]]$std)) A[[g]]$std <- Mstd$A
    if (length(S) > 0 && !is.null(S[[g]]$est) && is.null(S[[g]]$std)) S[[g]]$std <- Mstd$S
    
    # Extract matrices:
    if (length(A)>0) APars <- modMat2Pars(A[[g]],"->","A",symmetric=FALSE,vec=FALSE,Names,Names,group=paste("Group",g),exprsup="") else APars <- dumPars
    if (length(S)>0) SPars <- modMat2Pars(S[[g]],"<->","S",symmetric=TRUE,vec=FALSE,Names,Names,group=paste("Group",g),exprsup="") else SPars <- dumPars
    
    # Combine ParsS:
    Parss[[g]] <- rbind(APars,SPars)
    
    # Remove zeroes:
    Parss[[g]] <- Parss[[g]][Parss[[g]]$est!=0,]
  }
  
  Pars <- do.call(rbind,Parss)
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = c(manNames,latNames),
    manifest = c(manNames,latNames)%in%manNames,
    exogenous = NA,
    stringsAsFactors=FALSE)

  # Remove duplicates plus factor loadings betwen mans and lats of same name:
  Vars <- Vars[!duplicated(Vars$name),]
  Pars <- Pars[!(Pars$lhs==Pars$rhs&Pars$edge!="<->"),]
  
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


