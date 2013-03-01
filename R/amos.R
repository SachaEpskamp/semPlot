
semPlotModel_Amos <- function(object)
{
  ## Warnings:
  warning("(Residual) variances of Amos model is not yet supported")
  
  # Read characters:
  str <- readChar(object,nchars=file.info(object)$size)
  
  # Extract Estimates section:
  estLocs <- gregexpr('<div ntype="estimates" nodecaption="Estimates" caption="Estimates">',str)[[1]]
  nModel <- length(estLocs)
  Parss <- list()
  
  # Open and close div:
  open <- gregexpr("<div",str)[[1]]
  close <- gregexpr("</div>",str)[[1]]
  
  for (mod in 1:nModel)
  {  
    startSect <- which(open==estLocs[mod])
    # Find title:
    titleString <- substring(str,open[startSect+1],close[which(close>open[startSect+1])[1]])
    modName <- regmatches(titleString,regexpr("(?<=<h5>).*?(?=</h5>)",titleString,perl=TRUE))
    
    # Find close:
    nest <- 1
    curOpen <- startSect + 1
    curClose <- which(close>open[curOpen])[1]
    repeat{
      # If next is opened:
      if (open[curOpen] < close[curClose])
      {
        nest <- nest + 1
        curOpen <- curOpen + 1
      } else {
        # If next is closed:
        nest <- nest - 1
        if (nest==0) break
        curClose <- curClose + 1
      }
    }
    EstTabs <- substring(str,open[startSect],close[curClose] + 5)
    
    # Extract tables:
    Tabs <- readHTMLTable(EstTabs)
    
    # Find names of tables;
    Tabspl <- strsplit(EstTabs,split="<table>")[[1]]
    Names <- regmatches(Tabspl,gregexpr('(?<=nodecaption=").*(?=">)',Tabspl,perl=TRUE))[-length(Tabspl)]
    Names <- sapply(Names,function(x)x[length(x)])
    
    names(Tabs) <- Names
    
    # Regression weights:
    Reg <- Tabs[[which(grepl("regression",names(Tabs),ignore.case=TRUE))[1]]]
    Reg <-  as.data.frame(lapply(Reg,as.character),stringsAsFactors=FALSE)
    if (is.null(Reg$Estimate)) Reg$Estimate <- 1
    if (is.null(Reg$Label)) Reg$Label <- ""
    
    # Make Pars:
    # Define Pars:
    Pars <- data.frame(
      label = Reg$Label, 
      lhs = Reg[,3],
      edge = "->",
      rhs = Reg[,1],
      est = as.numeric(gsub(",",".",Reg$Estimate)),
      std = NA,
      group = modName,
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
#     Pars$par <- 1:nrow(Pars)
    Pars$label[is.na(Pars$label)] <- ""
    
#     # Fix edges:
#     Pars$edge[Reg[,2]=="<---"] <- "->"
#     Pars$edge[Reg[,2]=="<-->"] <- "<->"
    
    # Test for fixed:
    if (!is.null(Reg$P)) Pars$fixed <- is.na(Reg$P)
    
    # Standardized values:
    if (any(grepl("standardized",names(Tabs),ignore.case=TRUE)))
    {
      Std <- Tabs[[which(grepl("standardized",names(Tabs),ignore.case=TRUE))[1]]]
      Std <- as.data.frame(lapply(Std,as.character),stringsAsFactors=FALSE)
      if (is.null(Std$Estimate)) Std$Estimate <- 1
      Pars$std <- as.numeric(gsub(",",".",Std$Estimate))
    }
    
    # Add covariances:
    if (any(grepl("covariance",names(Tabs),ignore.case=TRUE)))
    {
      Cov <- Tabs[[which(grepl("covariance",names(Tabs),ignore.case=TRUE))[1]]]
      Cov <- as.data.frame(lapply(Cov,as.character),stringsAsFactors=FALSE)
      if (is.null(Cov$Estimate)) Cov$Estimate <- 1
      if (is.null(Cov$Label)) Cov$Label <- ""
      
      covPars <- data.frame(
        label = Cov$Label, 
        lhs = Cov[,3],
        edge = "<->",
        rhs = Cov[,1],
        est = as.numeric(gsub(",",".",Cov$Estimate)),
        std = NA,
        group = modName,
        fixed = FALSE,
        par = 0,
        stringsAsFactors=FALSE)
      
        if (!is.null(Cov$P)) covPars$fixed <- is.na(Cov$P)
      
        # Check cors:
      if (any(grepl("correlation",names(Tabs),ignore.case=TRUE)))
      {
        Cor <- Tabs[[which(grepl("correlation",names(Tabs),ignore.case=TRUE))[1]]]
        Cor <- as.data.frame(lapply(Cor,as.character),stringsAsFactors=FALSE)
        if (is.null(Cor$Estimate)) Cor$Estimate <- 1
        covPars$std <- Cor$Estimate
      }
      Pars <- rbind(Pars,covPars)
    }    
    
    Parss[[mod]] <- Pars
  }
  
  Pars <- do.call(rbind,Parss)
  Pars$par <- 1:nrow(Pars)
  
  ## Extract variable info:
  startSect <- which(open==gregexpr('<div nodecaption="Variable list"',str)[[1]])
  nest <- 1
  curOpen <- startSect + 1
  curClose <- which(close>open[curOpen])[1]
  repeat{
    # If next is opened:
    if (open[curOpen] < close[curClose])
    {
      nest <- nest + 1
      curOpen <- curOpen + 1
    } else {
      nest <- nest - 1
      if (nest==0) break
      curClose <- curClose + 1
    }
  }
  VarList <- substring(str,open[startSect],close[curClose] + 5)
  # Reove html tags:
  VarList <- gsub("<(.|\n)*?>","",VarList)
  # Per line:
  VarList <-  scan(text=VarList,what="character",sep="\n")
  # Remove leading and trailing whitespace:
  VarList <- gsub("^\\s*","",VarList)
  VarList <- gsub("\\s*$","",VarList)
  
  AllVars <- unique(c(Pars$lhs,Pars$rhs))
  AllVars <- AllVars[AllVars!=""]

  # Location of indicators:
  # manEndo - latEndo - manExo - latExo
  grepRep <- function(...)
  {
    res <- grep(...)
    if (length(res)==0) res <- -1
    return(res[1])
  }
  ind <- c(
    grepRep("Observed, endogenous variables",VarList),
    grepRep("Unobserved, endogenous variables",VarList),
    grepRep("Observed, exogenous variables",VarList),
    grepRep("Unobserved, exogenous variables",VarList))
  
  if (ind[1] > 0)
  {
    manEndo <- VarList[(ind[1]+1):(which((1:length(VarList) == length(VarList)) | (1:length(VarList) > ind[1]+1 & 1:length(VarList) %in% ind))[1] - 1)]
  } else manEndo <- character(0)
  
  if (ind[2] > 0)
  {
    latEndo <- VarList[(ind[2]+1):(which((1:length(VarList) == length(VarList)) | (1:length(VarList) > ind[2]+1 & 1:length(VarList) %in% ind))[1] - 1)]
  } else latEndo <- character(0)
  
  if (ind[3] > 0)
  {
    manExo <- VarList[(ind[3]+1):(which((1:length(VarList) == length(VarList)) | (1:length(VarList) > ind[3]+1 & 1:length(VarList) %in% ind))[1] - 1)]
  } else manExo <- character(0)
  
  if (ind[4] > 0)
  {
    latExo <- VarList[(ind[4]+1):(which((1:length(VarList) == length(VarList)) | (1:length(VarList) > ind[4]+1 & 1:length(VarList) %in% ind))[1] - 1)]
  } else latExo <- character(0)
  
  Vars <- data.frame(
    name = c(manEndo,manExo,latEndo,latExo),
    manifest = c(rep(TRUE,length(c(manEndo,manExo))),rep(FALSE,length(c(latEndo,latExo)))),
    exogenous = c(rep(FALSE,length(manEndo)),rep(TRUE,length(manExo)),rep(FALSE,length(latEndo)),rep(TRUE,length(latExo))),
    stringsAsFactors=FALSE)
  
  
  
  # Return:
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- FALSE
  semModel@Original <- list(str)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  # semModel@Thresholds <- Thresh
  
  return(semModel)
}