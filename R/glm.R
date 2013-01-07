
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
  
  if (is.null(rownames(coef))) rownames(coef) <- combLetters(1:Nr)
  if (is.null(colnames(coef))) colnames(coef) <- combLetters(1:Nc)
  
  NamesR <- rownames(coef)
  NamesC <- colnames(coef)
  
  RAM  <- data.frame(
    label = "", 
    lhs = rep(NamesR,times=Nc),
    edge = "->",
    rhs = rep(NamesC,each=Nr),
    est = c(coef),
    std = NA,
    group = "",
    fixed = FALSE,
    par = 1:(Nr*Nc),
    stringsAsFactors=FALSE)
  
  RAM$edge[grepl("intercept",RAM$lhs,ignore.case=TRUE)] <- "int"
  RAM$lhs[grepl("intercept",RAM$lhs,ignore.case=TRUE)] <- ""
  
  ## Remove interactions plus warning:
  if (any(grepl("\\:",RAM$lhs)))
  {
    warning("Interaction terms are not yet supported by semPlot, and are removed.")
    RAM <- RAM[!grepl("\\:",RAM$lhs),]  
  }
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = c(NamesR[!grepl("intercept",NamesR,ignore.case=TRUE)],NamesC),
    manifest = TRUE,
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  semModel <- new("semPlotModel")
  semModel@RAM <- RAM
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  
  return(semModel)
}