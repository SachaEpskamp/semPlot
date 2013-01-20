
semPlotModel.lm <- function(object)
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
  
  NamesR <- rownames(coef)
  NamesC <- colnames(coef)
  
  RAM  <- data.frame(
    label = "", 
    lhs = rep(NamesR,times=Nc),
    edge = "->",
    rhs = rep(NamesC,each=Nr),
    est = c(coef),
    std = c(coef(standardize(object))),
    group = "",
    fixed = FALSE,
    par = 1:(Nr*Nc),
    knot = 0,
    stringsAsFactors=FALSE)
  
  ## Split interactions:
  if (any(grepl(":",RAM$lhs)))
  {
    colons <- grep(":",RAM$lhs)
    for (i in seq_along(colons))
    {
      labs <- strsplit(RAM$lhs[colons[i]],split=":")[[1]]
      RAM$lhs[colons[i]] <- labs[1]
      RAM$knot[colons[i]] <- i
      for (j in 2:length(labs))
      {
        RAM <- rbind(RAM,RAM[colons[i],])
        RAM$lhs[nrow(RAM)] <- labs[j]
      }
    }
  }
  
  RAM$edge[grepl("intercept",RAM$lhs,ignore.case=TRUE)] <- "int"
  RAM$lhs[grepl("intercept",RAM$lhs,ignore.case=TRUE)] <- ""
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = unique(c(RAM$lhs,RAM$rhs)),
    manifest = TRUE,
    exogenous = NA,
    stringsAsFactors=FALSE)
  Vars <- Vars[Vars$name!="",]
  
  semModel <- new("semPlotModel")
  semModel@RAM <- RAM
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  
  return(semModel)
}