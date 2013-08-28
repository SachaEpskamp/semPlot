# 
# # object <- readModels(file.choose())
# semPaths.mplus.model <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }

semPlotModel.mplus.model <- function(object, ...)
{
  # Check for mplusAutomation:
  if (!require("MplusAutomation")) stop("'MplusAutomation' package must be installed to read Mplus output.")
  
  addInteractions <- FALSE
  if (is.character(object))
  {
    modfile <- object
    object <- readModels(object)
    
    mod <- readLines(modfile)
    # Find XWITH:
    xs <- grep("XWITH",mod)
    if (length(xs)>0)
    {
      # Split:
      spl <- strsplit(mod[xs],split="\\|")
      # Find vars that interact:
      vars <- lapply(spl,function(x)strsplit(x[2],split="XWITH")[[1]])
      # Extract Vars
      newvars <- sapply(spl,'[',1)
      # sanitize:
      newvars <- gsub("\\W","",newvars)
      vars <- lapply(vars,gsub,pattern="\\W",replacement="")
      
      addInteractions <- TRUE
    }
    
  } else warning("Interactions are ommited. Use semPlotModel.mplus.model on the path to mplus output file for semPlot to attempt to find assigned interactions.")
  
  if (length(object$parameters)==0) stop("No parameters detected in mplus output.")
  
  parsUS <- object$parameters$unstandardized
  if (is.null(parsUS$Group)) parsUS$Group <- ""
  if (is.null(parsUS$BetweenWithin)) parsUS$BetweenWithin <- ""
  
  if (any(grepl("\\|",parsUS$paramHeader)))
  {
    parsUS$paramHeader <- gsub("\\|", "BY", parsUS$paramHeader)
    warning("'|' operator replaced by BY operator.")
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
  
  # Define Pars:
  Pars <- data.frame(
    label = "", 
    lhs = "",
    edge = "--",
    rhs = parsUS$param,
    est = parsUS$est,
    std = NA,
    group = parsUS$Group,
    fixed = parsUS$se==0,
    par = 0,
    BetweenWithin = parsUS$BetweenWithin,
    stringsAsFactors=FALSE)
  
  if (!noPars)
  {
    parNums <- dlply(cbind(sapply(parsUS[c("est","se")],function(x)round(as.numeric(x),10)),data.frame(num=1:nrow(parsUS))),c("est","se"),'[[',"num")
    for (i in 1:length(parNums)) Pars$par[parNums[[i]]] <- i
    Pars$par[Pars$fixed] <- 0  
  } else Pars$par <- 1:nrow(Pars)
  
#   
#   c <- 1
#   for (i in 1:nrow(Pars))
#   {
#     if (!isTRUE(Pars$fixed[i]) & Pars$par[i]==0)
#     {
#       par <- sapply(1:nrow(parsUS),function(j)isTRUE(all.equal(unlist(parsUS[j,c("est","se","est_se","pval")]),unlist(parsUS[i,c("est","se","est_se","pval")]))))
#       Pars$par[par] <- c
#       c <- c+1
#     }
#   }
  
  
  if (!is.null(object$parameters$std.standardized))
  {
    Pars$std <- object$parameters$std.standardized$est
  }
  
  Pars$lhs[grepl(".BY$",parsUS$paramHeader)] <- gsub("\\.BY$","",parsUS$paramHeader[grepl(".BY$",parsUS$paramHeader)])
  Pars$edge[grepl(".BY$",parsUS$paramHeader)] <- "->"
  
  Pars$lhs[grepl(".ON$",parsUS$paramHeader)] <- gsub("\\.ON$","",parsUS$paramHeader[grepl(".ON$",parsUS$paramHeader)])
  Pars$edge[grepl(".ON$",parsUS$paramHeader)] <- "~>"
  Pars[grepl(".ON$",parsUS$paramHeader),c("lhs","rhs")] <- Pars[grepl(".ON$",parsUS$paramHeader),c("rhs","lhs")]
  
  Pars$lhs[grepl(".WITH$",parsUS$paramHeader)] <- gsub("\\.WITH$","",parsUS$paramHeader[grepl(".WITH$",parsUS$paramHeader)])
  Pars$edge[grepl(".WITH$",parsUS$paramHeader)] <- "<->"
  
  Pars$lhs[grepl("Variances",parsUS$paramHeader)] <- Pars$rhs[grepl("Variances",parsUS$paramHeader)]
  Pars$edge[grepl("Variances",parsUS$paramHeader)] <- "<->"
  
  Pars$edge[grepl("Means|Intercepts",parsUS$paramHeader)] <- "int"
  
  if (!is.null(object$parameters$standardized)) Pars$std <- object$parameters$standardized$est
  
  
  # Extract threshold model:
  Thresh <- Pars[grepl("Thresholds",parsUS$paramHeader),-(3:4)]
  Thresh$lhs <- gsub("\\$.*","",Pars$rhs[grepl("Thresholds",parsUS$paramHeader)])
  Thresh$BetweenWithin[Thresh$BetweenWithin == "Between"] <- "Within"
  Pars <- Pars[!grepl("Thresholds",parsUS$paramHeader),]
  
  # Detect latent/manifest:
  Latents <- unique(gsub("\\.BY$","",parsUS$paramHeader[grepl(".BY$",parsUS$paramHeader)]))
  var <- unique(unlist(Pars[c("lhs","rhs")]))
  var <- var[var!=""]
  
  # Variable dataframe: 
  Vars <- data.frame(
    name = var,
    manifest = !var%in%Latents,
    exogenous = NA,
    stringsAsFactors=FALSE)
  
  
  ### Add interactions and remove dummy variables:
  if (addInteractions)
  {
    Vars <- Vars[!tolower(Vars$name)%in%tolower(newvars),]
    
    Pars$knot <- 0
    k <- 1
    for (i in rev(seq_along(newvars)))
    {
      varlocs <- which(tolower(Pars$lhs)==tolower(newvars[i])|tolower(Pars$rhs)==tolower(newvars[i]))
      for (v in seq_along(varlocs))
      {
        for (j in 1:length(vars[[i]]))
        {
          Parsnew <- Pars[varlocs[v],]
          Parsnew$lhs[tolower(Parsnew$lhs)==tolower(newvars[i])] <- Vars$name[match(tolower(vars[[i]][j]),tolower(Vars$name))]
          Parsnew$rhs[tolower(Parsnew$rhs)==tolower(newvars[i])] <- Vars$name[match(tolower(vars[[i]][j]),tolower(Vars$name))]
          if (Parsnew$knot==0)
          {
            Parsnew$knot <- k
          }
          Pars <- rbind(Pars,Parsnew)
        } 
        if (any(Pars$knot==k)) k <- k + 1
      }
      Pars <- Pars[-varlocs,]
    }
    
  }

  
  # Abbreviate names with more than 8 characters:
  Pars$lhs <- substring(Pars$lhs,1,8)
  Pars$rhs <- substring(Pars$rhs,1,8)
  Vars$name <- substring(Vars$name,1,8)
  Vars <- Vars[!duplicated(Vars),]
  
  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  semModel@Thresholds <- Thresh
  
  return(semModel)
}