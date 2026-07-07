
semPlotModel_Onyx <- function(object)
{
  if (!requireNamespace("XML", quietly = TRUE))
  {
    stop("The 'XML' package is required to import Onyx models. Install it with install.packages('XML').")
  }
  # Parse Onyx model:
  doc <- XML::xmlParse(object)
  
  # Get Nodes and Edges:
  Nodes <- XML::getNodeSet(doc, "/model/graph/node")
  Edges <- XML::getNodeSet(doc, "/model/graph/edge")
  Const <- as.logical(sapply(Nodes, function(n) XML::xmlGetAttr(n, "constant")))
  
  # Get NodeNames:
  NodeNames <- sapply(Nodes, function(n) XML::xmlGetAttr(n, "caption"))
  NodeNames[Const] <- ""
  
  # Get edgelist:
  Edgelist <- data.frame(
    From = as.numeric(as.character(sapply(Edges, function(n) XML::xmlGetAttr(n, "sourceNodeId")))),
    To = as.numeric(as.character(sapply(Edges, function(n) XML::xmlGetAttr(n, "targetNodeId")))),  
    stringsAsFactors=FALSE) + 1
  
  # Define Pars:
  Pars <- data.frame(
    label = sapply(Edges, function(n) XML::xmlGetAttr(n, "parameterName")), 
    lhs = NodeNames[Edgelist$From],
    edge = ifelse(as.logical(sapply(Edges, function(n) XML::xmlGetAttr(n, "doubleHeaded"))),"<->","->"),
    rhs = NodeNames[Edgelist$To],
    est = as.numeric(as.character(sapply(Edges, function(n) XML::xmlGetAttr(n, "value")))),
    std = NA,
    group = "",
    fixed = as.logical(sapply(Edges, function(n) XML::xmlGetAttr(n, "fixed"))),
    par = 0,
    stringsAsFactors=FALSE)
  
  Pars$edge[Pars$lhs==""] <- "int"
  Pars$par <- 1:nrow(Pars)
  
  # Vars:
  Vars <- data.frame(
    name = NodeNames,
    manifest = !as.logical(sapply(Nodes, function(n) XML::xmlGetAttr(n, "latent"))),
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  Vars <- Vars[c(which(Vars$manifest),which(!Vars$manifest)),]
  
  # Return:
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- FALSE
  semModel@Original <- list(doc)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  # semModel@Thresholds <- Thresh
  
  return(semModel)
}