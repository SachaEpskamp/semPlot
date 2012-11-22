
# object <- readModels(file.choose())
semPaths.mplus.model <- function(object,...) 
{
  invisible(semPaths(semPlotModel(object),...))
}

semPlotModel.mplus.model <- function(object)
{
  if (length(object$parameters)==0) stop("No parameters detected in mplus output.")
  
  parsUS <- object$parameters$unstandardized
  if (is.null(parsUS$Group)) parsUS$Group <- ""
  if (is.null(parsUS$BetweenWithin)) parsUS$BetweenWithin <- ""
  
  if (any(grepl("\\|",parsUS$paramHeader)))
  {
    parsUS <- parsUS[!grepl("\\|",parsUS$paramHeader),]
    warning("'|' operator is not yet supported by semPlot. Parameters using this operator will not be shown and unexpected results might occur.")
  }
  
  if (any(grepl("New.Additional.Parameters",parsUS$paramHeader)))
  {
    parsUS <- parsUS[!grepl("New.Additional.Parameters",parsUS$paramHeader),]
    warning("'New.Additional.Parameters' is not yet supported by semPlot. Parameters will not be shown and unexpected results might occur.")
  }
  
  noPars <- FALSE
  # Temporary fix for EFA:
  if (is.null(parsUS$est))
  {
    if (!is.null(parsUS$average)) 
    {
      parsUS$est <- parsUS$average
      parsUS$se <- parsUS$average_se
    } else 
    {
      parsUS$est <- 0
      parsUS$se <- 0
      noPars <- TRUE
    }
  }
  
  # Define RAM:
  RAM <- data.frame(
    label = "", 
    lhs = "",
    edge = "--",
    rhs = parsUS$param,
    est = parsUS$est,
    std = NA,
    group = paste(parsUS$Group, parsUS$BetweenWithin, sep=""),
    fixed = parsUS$se==0,
    par = 0,
    stringsAsFactors=FALSE)
  
  if (!noPars)
  {
    parNums <- dlply(cbind(sapply(parsUS[c("est","se")],function(x)round(as.numeric(x),10)),data.frame(num=1:nrow(parsUS))),c("est","se"),'[[',"num")
    for (i in 1:length(parNums)) RAM$par[parNums[[i]]] <- i
    RAM$par[RAM$fixed] <- 0  
  } else RAM$par <- 1:nrow(RAM)
  
#   
#   c <- 1
#   for (i in 1:nrow(RAM))
#   {
#     if (!isTRUE(RAM$fixed[i]) & RAM$par[i]==0)
#     {
#       par <- sapply(1:nrow(parsUS),function(j)isTRUE(all.equal(unlist(parsUS[j,c("est","se","est_se","pval")]),unlist(parsUS[i,c("est","se","est_se","pval")]))))
#       RAM$par[par] <- c
#       c <- c+1
#     }
#   }
  
  
  if (!is.null(object$parameters$std.standardized))
  {
    RAM$std <- object$parameters$std.standardized$est
  }
  
  RAM$lhs[grepl("BY",parsUS$paramHeader)] <- gsub("\\.BY","",parsUS$paramHeader[grepl("BY",parsUS$paramHeader)])
  RAM$edge[grepl("BY",parsUS$paramHeader)] <- "->"
  
  RAM$lhs[grepl("ON",parsUS$paramHeader)] <- gsub("\\.ON","",parsUS$paramHeader[grepl("ON",parsUS$paramHeader)])
  RAM$edge[grepl("ON",parsUS$paramHeader)] <- "~>"
  RAM[grepl("ON",parsUS$paramHeader),c("lhs","rhs")] <- RAM[grepl("ON",parsUS$paramHeader),c("rhs","lhs")]
  
  RAM$lhs[grepl("WITH",parsUS$paramHeader)] <- gsub("\\.WITH","",parsUS$paramHeader[grepl("WITH",parsUS$paramHeader)])
  RAM$edge[grepl("WITH",parsUS$paramHeader)] <- "<->"
  
  RAM$lhs[grepl("Variances",parsUS$paramHeader)] <- RAM$rhs[grepl("Variances",parsUS$paramHeader)]
  RAM$edge[grepl("Variances",parsUS$paramHeader)] <- "<->"
  
  RAM$edge[grepl("Means|Intercepts",parsUS$paramHeader)] <- "int"
  
  if (!is.null(object$parameters$standardized)) RAM$std <- object$parameters$standardized$est
  
  
  # Extract threshold model:
  Thresh <- RAM[grepl("Thresholds",parsUS$paramHeader),-(3:4)]
  Thresh$lhs <- gsub("\\$.*","",RAM$rhs[grepl("Thresholds",parsUS$paramHeader)])
  RAM <- RAM[!grepl("Thresholds",parsUS$paramHeader),]
  
  # Detect latent/manifest:
  Latents <- unique(gsub("\\.BY","",parsUS$paramHeader[grepl("BY",parsUS$paramHeader)]))
  var <- unique(unlist(RAM[c("lhs","rhs")]))
  var <- var[var!=""]
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = var,
    manifest = !var%in%Latents,
    stringsAsFactors=FALSE)
  

  
  semModel <- new("semPlotModel")
  semModel@RAM <- RAM
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  semModel@Thresholds <- Thresh
  
  return(semModel)
}