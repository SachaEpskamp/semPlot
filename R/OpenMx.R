### Path diagrams ###
# 
# semPaths_MxRAMModel <- function(object,...){
#   invisible(semPaths(semPlotModel(object),...))
# }
#           
# semPaths_MxModel <- function(object,...){
#   invisible(semPaths(semPlotModel(object),...))
# }
#  
### EXTRACT MODEL ###
          
### SINGLE GROUP ###
semPlotModel_MxRAMModel <- function(object){

  # Extract names:
  varNames <- object@manifestVars
  factNames <- object@latentVars
  if (!(length(varNames) || length(factNames))) 
    stop(as.character(substitute(object)), '@manifestVars (and ',
         as.character(substitute(object)),  '@latentVars if the model has ',
         'latent variables) must contain variable names. You can set them ',
         'using the manifestVars= and latentVars= arguments in mxModel().')

  # Standardized object:
  std <- OpenMx::mxStandardizeRAMpaths(object, SE = TRUE)

  # Extract directed paths:
#   Dirpaths <- which(t(object@matrices$A@free | object@matrices$A@values!=0),arr.ind=TRUE)
#   DirpathsFixed <- !t(object@matrices$A@free)[Dirpaths]
#   DirpathsValues <- t(object@matrices$A@values)[Dirpaths]
#   DirpathsLabels <- t(object@matrices$A@labels)[Dirpaths]
  
  # Extract symmetric paths:
#   Sympaths <- which(t(object@matrices$S@free | object@matrices$S@values!=0) & upper.tri(object@matrices$S@values,diag=TRUE),arr.ind=TRUE)
#   SympathsFixed <- !t(object@matrices$S@free)[Sympaths]
#   SympathsValues <- t(object@matrices$S@values)[Sympaths]
#   SympathsLabels <- t(object@matrices$A@labels)[Sympaths]
  
#   if (!is.null(object@matrices$M))
#   {
#     # Extract intercepts:
#     Means <- which(object@matrices$M@free | object@matrices$M@values!=0)
#     MeansFixed <- !object@matrices$M@free[Means]
#     MeansValues <- object@matrices$M@values[Means]
#     MeansLabels <- object@matrices$M@labels[Means]
#   } else
#   {
#     Means <- numeric(0)
#     MeansFixed <- logical(0)
#     MeansValues <- numeric(0)
#     MeansLabels <- character(0)
#   }
#   
#   ## Standardized
#   if (!length(object@output)==0)
#   {
#     # browser()
#     # Function by Ryne Estabrook (http://openmx.psyc.virginia.edu/thread/718)
#     
# standObj <- standardizeRAM(object,"model")
#     
#     # Extract directed paths:
#     # DirpathsValuesStd <- t(standObj@matrices$A@values)[Dirpaths]
#     # DirpathsValuesStd <- std$Std.Value[std$matrix=="A"]
#     
#       # Extract symmetric paths:
#     SympathsValuesStd <- t(standObj@matrices$S@values)[Sympaths]
#       
#       # Extract means:
#     
#if (!is.null(standObj@matrices$M))
#    {
#    MeansValuesStd <- standObj@matrices$S@values[Means]
#   } else    {
#      MeansValuesStd <- numeric(0)
#    }
#   } else 
#   {
#     DirpathsValuesStd <- rep(NA,nrow(Dirpaths)) 
#     SympathsValuesStd <- rep(NA,nrow(Sympaths))
#     MeansValuesStd <- rep(NA,length(Means))
#   }
#   
  # Vars dataframe:
  Vars <- data.frame(
    name = c(varNames,factNames),
    manifest = c(varNames,factNames)%in%varNames,
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  
  
  
 # standObj <- standardizeMx(object,free=T) # old semTools function, now in this file
  
  Edges <- std
  
  # Only edges in mats A and S:
  corMats <- Edges$matrix %in% c("A","S")
  
  # Define Pars:
  Pars <- data.frame(
    label = ifelse(is.na(Edges$label[corMats]),"",Edges$label[corMats]), 
    lhs = Edges$col[corMats],
    edge = ifelse(Edges$matrix[corMats]=="A","->","<->"),
    rhs = Edges$row[corMats],
    est = Edges$Raw.Value[corMats],
    std = Edges$Std.Value[corMats],
    group = '',
    fixed = Edges$Raw.SE[corMats]==0,
    par = 0,
    stringsAsFactors=FALSE)
  
  

  # Maybe remove ints?
  if (!is.null(object@matrices$M)) {
    MeanStd <- c(object@matrices$M$values)
    
    ## in case labels are NA, use variable names
    if (!is.null(colnames(object@matrices$M$values))) {
      v.names <- colnames(object@matrices$M$values)
      v.idx <- v.names
      names(MeanStd) <- v.names
    } else {
      #FIXME? Warn users that this assumes order is {all manifest, all latent}
      v.names <- c(varNames, factNames) 
      v.idx <- seq_along(v.names)
    }
    ## extract rows of std corresponding to the M matrix
    stdM <- std[std$matrix == "M", , drop = FALSE]
    ## loop over variable names that have a standardized estimate
    ## (only free parameters; assume others are fixed to zero)
    for (v in seq_along(stdM$col)) {
      MeanStd[ v.idx[v] ] <- stdM$Std.Value[stdM$col == v.idx[v] ]
    }
    ## old method (using deprecated semTools function, now at the bottom of this script)
    ## standardizeMx(object,free=T)[which(names(standardizeMx(object,free=T))%in%object@matrices$M$labels)]
    
    MeanEst <-data.frame(
      label = c(object@matrices$M$labels),
      ##### or, if they are NA, replace with variable names?
      # label = ifelse(!is.na(object@matrices$M$labels),
      #                yes = object@matrices$M$labels,
      #                no = v.names),
      lhs = '',
      rhs = v.names,
      edge = 'int',
      est = c(object@matrices$M$values),
      std = MeanStd,
      group = '',
      fixed = c(!object@matrices$M$free),
      par = 0,
      stringsAsFactors = FALSE )
    Pars <- rbind(Pars,MeanEst)
  }

  
  
  
  Pars$par[is.na(Pars$label)] <- seq_len(sum(is.na(Pars$label)))
  for (lbl in unique(Pars$label[!is.na(Pars$label)]))
  {
    Pars$par[Pars$label==lbl] <- max(Pars$par)+1
  }
#   
#   # Add standardized:
#   for (i in 1:nrow(standPars))
#   {
#     if (standPars$matrix[i] == "A")
#     {
#       Pars$std[Pars$lhs == standPars$col[i] & Pars$rhs == standPars$row[i] & Pars$edge == "->"] <- standPars[["Std. Estimate"]][i]
#     }
#     if (standPars$matrix[i] == "S")
#     {
#       Pars$std[Pars$lhs == standPars$col[i] & Pars$rhs == standPars$row[i] & Pars$edge == "<->"] <- standPars[["Std. Estimate"]][i]
#     }
#   }
  
  Pars$label[is.na(Pars$label)] <- ""
  
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- !length(object@output)==0
  semModel@Original <- list(object)
  
  if (!is.null(object@data))
  {
    if (object@data@type=="cov")
    {
      semModel@ObsCovs <- list(object@data@observed)
    } else if (object@data@type=="raw")
    {
      semModel@ObsCovs <- list(cov(object@data@observed))
    } else
    {
      semModel@ObsCovs <- list(NULL)
    }
  } else
  {
    semModel@ObsCovs <- list(NULL)
  }
  
  
  semModel@ImpCovs <- list(object@fitfunction@info$expCov)
  
  return(semModel)
}


semPlotModel_MxModel <- function(object){

  if (any(!"MxRAMModel"%in%sapply(object@submodels,class))) stop("Model or all submodels must be of class 'MxRAMModel'")
  for (i in 1:length(object@submodels)) object@submodels[[i]]@output <- list(TRUE)
  S4objects <- lapply(object@submodels,semPlotModel)
  
  semModel <- new("semPlotModel")
  semModel@Pars <- do.call("rbind",lapply(S4objects,slot,"Pars"))
  
  semModel@Pars$par <- 0
  semModel@Pars$par[semModel@Pars$label==""] <- seq_len(sum(semModel@Pars$label==""))
  for (lbl in unique(semModel@Pars$label[semModel@Pars$label!=""]))
  {
    semModel@Pars$par[semModel@Pars$label==lbl] <- max(semModel@Pars$par)+1
  }
  
  semModel@Vars <- S4objects[[1]]@Vars
  semModel@Computed <- !length(object@output)==0
  semModel@Original <- list(object)
  
  semModel@ObsCovs <- lapply(S4objects,function(x)x@ObsCovs[[1]])
  names(semModel@ObsCovs) <- sapply(object@submodels,slot,"name")
  
  
  semModel@ImpCovs <- lapply(S4objects,function(x)x@ImpCovs[[1]])
  names(semModel@ImpCovs) <- sapply(object@submodels,slot,"name")
  
  
  return(semModel)
}





## -----------------------------------------------------------------
## semTools function (no longer used, but can be borrowed if needed)
## -----------------------------------------------------------------

standardizeMx <- function(object, free = TRUE) {
  .Deprecated(msg = c("The standardizeMx function is deprecated, and it will",
                      " cease to be included in future versions of semTools.",
                      " See help('semTools-deprecated) for details."))
  # objectOrig <- object
  multigroup <- length(object@submodels) > 0
  if(multigroup) {
    defVars <- lapply(object@submodels, findDefVars)
    defVars <- do.call(c, defVars)
  } else {
    defVars <- findDefVars(object)
  }
  if(length(defVars) > 0) stop("The standardizeMx is not available for the model with definition variable.")
  if(multigroup) {
    object@submodels <- lapply(object@submodels, standardizeMxSingleGroup)
  } else {
    object <- standardizeMxSingleGroup(object)
  }
  vectorizeMx(object, free=free)
}

## Hidden functions

findDefVars <- function(object) {
  ## borrowed from OpenMx::imxIsDefinitionVariable
  imxSeparatorChar <- "."
  imxIsDefinitionVariable <- function (name) {
    if (is.na(name)) {
      return(FALSE)
    }
    components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
    if (length(components) == 2 && components[[1]] == "data") {
      return(TRUE)
    }
    else if (length(components) > 2 && components[[2]] == "data") {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
  }
  ## end borrowed code
  mat <- lapply(object@matrices, slot, "labels")
  defvars <- sapply(mat, function(x) x[apply(x, c(1,2), imxIsDefinitionVariable)])
  Reduce("c", defvars)
}

vectorizeMx <- function(object, free = TRUE) {
  multigroup <- length(object@submodels) > 0
  if(multigroup) {
    object <- object@submodels
  } else {
    object <- list(object)
  }
  result <- NULL
  for(i in seq_along(object)) {
    name <- ""
    if(multigroup) name <- paste0(object[[i]]@name, ".")
    mat <- object[[i]]@matrices
    for(j in seq_along(mat)) {
      tempname <- paste0(name, mat[[j]]@name)
      lab <- mat[[j]]@labels
      tempfree <- as.vector(mat[[j]]@free)
      madeLab <- paste0(tempname, "[", row(lab), ",", col(lab), "]")
      lab <- as.vector(lab)
      madeLab[!is.na(lab)] <- lab[!is.na(lab)]
      if(!free) tempfree <- rep(TRUE, length(tempfree))
      temp <- mat[[j]]@values[tempfree]
      names(temp) <- madeLab[tempfree]
      result <- c(result, temp)
    }
  }
  
  result[!duplicated(names(result))]
}

standardizeMxSingleGroup <- function(object) {
  if (!is(object@expectation, "MxExpectationRAM"))
    stop("The standardizeMx function is available for the MxExpectationRAM only.")
  A <- object@matrices$A@values
  I <- diag(nrow(A))
  S <- object@matrices$S@values
  # F <- object@matrices$F@values
  Z <- solve(I - A)
  impliedCov <- Z %*% S %*% t(Z)
  temp <- sqrt(diag(impliedCov))
  if (length(temp) == 1) {
    ImpliedSd <- as.matrix(temp)
  } else {
    ImpliedSd <- diag(temp)
  }
  ImpliedInvSd <- solve(ImpliedSd)
  object@matrices$S@values <- ImpliedInvSd %*% S %*% ImpliedInvSd
  object@matrices$A@values <- ImpliedInvSd %*% A %*% ImpliedSd
  if (!is.null(object@matrices$M)) {
    M <- object@matrices$M@values
    object@matrices$M@values <- M %*% ImpliedInvSd
  }
  object
}



