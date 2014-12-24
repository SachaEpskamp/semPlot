
semPlotModel.lm <- function(object, ...)
{
  coef <- as.matrix(coef(object))
  Nr <- nrow(coef)
  Nc <- ncol(coef)
  
  combLetters <- function(x)
  {
    if (length(x)>1) return(sapply(x,combLetters))
    
    f <- function(x)
    {  
      if (x[1]>26)  c(f(floor(x/26)),x%%26 + 1) else x
    }
    
    paste(LETTERS[f(x)],collapse="")
  }
  
  if (is.null(rownames(coef)))
  {
    rownames(coef) <- names(object$model)[(Nc+1):length(object$model)] 
  }
  
  if (is.null(colnames(coef)))
  {
    colnames(coef) <- names(object$model)[1:Nc]
  }
  
  namesCoef <- rownames(coef)
  stdCoef <- coef(standardize(object))
  names(stdCoef) <- gsub("`","",names(stdCoef))
  
  NamesR <- rownames(coef)
  NamesC <- colnames(coef)

  
  Pars  <- data.frame(
    label = "", 
    lhs = rep(NamesR,times=Nc),
    edge = "->",
    rhs = rep(NamesC,each=Nr),
    est = c(coef),
    std = unname(c(stdCoef[paste0(namesCoef,"s")])),
    group = "",
    fixed = FALSE,
    par = 1:(Nr*Nc),
    knot = 0,
    stringsAsFactors=FALSE)
  
  ## Split interactions:
  if (any(grepl(":",Pars$lhs)))
  {
    colons <- grep(":",Pars$lhs)
    for (i in seq_along(colons))
    {
      labs <- strsplit(Pars$lhs[colons[i]],split=":")[[1]]
      Pars$lhs[colons[i]] <- labs[1]
      Pars$knot[colons[i]] <- i
      for (j in 2:length(labs))
      {
        Pars <- rbind(Pars,Pars[colons[i],])
        Pars$lhs[nrow(Pars)] <- labs[j]
      }
    }
  }
  
  Pars$edge[grepl("intercept",Pars$lhs,ignore.case=TRUE)] <- "int"
  Pars$lhs[grepl("intercept",Pars$lhs,ignore.case=TRUE)] <- ""
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = unique(c(Pars$lhs,Pars$rhs)),
    manifest = TRUE,
    exogenous = NA,
    stringsAsFactors=FALSE)
  Vars <- Vars[Vars$name!="",]
  
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  
  return(semModel)
}