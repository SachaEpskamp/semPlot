
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
  # standardize() fails for non-gaussian GLMs (e.g. binomial) because
  # standardizing the response breaks the link function constraints.
  # Fall back to raw coefficients in that case.
  std_ok <- TRUE
  if (requireNamespace("rockchalk", quietly = TRUE))
  {
    stdCoef <- tryCatch(
      coef(rockchalk::standardize(object)),
      error = function(e) {
        warning("Could not compute standardized coefficients: ", e$message,
                ". Using raw coefficients instead.")
        std_ok <<- FALSE
        coef(object)
      }
    )
  } else {
    warning("The 'rockchalk' package is not installed; using raw coefficients instead of standardized ones.")
    std_ok <- FALSE
    stdCoef <- coef(object)
  }
  names(stdCoef) <- gsub("`","",names(stdCoef))
  
  NamesR <- rownames(coef)
  NamesC <- colnames(coef)

  # standardize() renames every variable with an "s" suffix (e.g. "x" -> "xs"),
  # so on success we look up "<name>s"; in the raw-coefficient fallback the names
  # are unchanged, so we look up "<name>" directly. rockchalk drops the intercept
  # from the standardized fit (no "(Intercept)s"), yielding NA there, so we force
  # the intercept to NA in the fallback too to stay consistent.
  if (std_ok)
  {
    stdVals <- stdCoef[paste0(namesCoef,"s")]
  } else {
    stdVals <- stdCoef[namesCoef]
    stdVals[grepl("intercept", namesCoef, ignore.case = TRUE)] <- NA
  }

  Pars  <- data.frame(
    label = "",
    lhs = rep(NamesR,times=Nc),
    edge = "->",
    rhs = rep(NamesC,each=Nr),
    est = c(coef),
    std = unname(c(stdVals)),
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
  
  isInt <- Pars$lhs == "(Intercept)" | tolower(Pars$lhs) == "intercept"
  Pars$edge[isInt] <- "int"
  Pars$lhs[isInt] <- ""
  
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