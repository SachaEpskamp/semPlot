## This function is commented out because semspec is not yet on CRAN. For the full version of this function please see
## www.sachaepskamp.com/semPlot

semPlotModel.semspec <- function(object)
{
  
  stop("This function is not included in the CRAN release because semspec is not on CRAN. Please see www.sachaepskamp.com for the function")
#   
#   # Load 'semspec':
#   if (!require("semspec")) stop('semspec is required: install.packages("semspec", repos="http://R-Forge.R-project.org")')
#   
#   semreprObject <- semrepr(object)
#   sumObject <- summary(object)
# 
#   # Define Pars:
#   Pars <- data.frame(
#     label = "", 
#     lhs = semreprObject$lhs,
#     edge = "--",
#     rhs = semreprObject$rhs,
#     est = NA,
#     std = NA,
#     group = ifelse(is.na(semreprObject$group),"",semreprObject$group),
#     fixed = FALSE,
#     par = 0,
#     stringsAsFactors=FALSE)
#   
#   # Label:
#   if (!is.null(semreprObject$param)) Pars$label <- semreprObject$param
#   
#   
#   # Fixed:
# #   if (!is.null(semreprObject$free)) Pars$fixed <- !semreprObject$free
#   if (length(sumObject$constraints$details$Constraint)>0)
#   {
#     spl <- strsplit(sumObject$constraints$details$Constraint,split=" == ")[grepl("==",sumObject$constraints$details$Constraint)]
#     parNum <- sapply(spl,function(x)sum(x%in%Pars$label))
#     parIt <- 1
#     for (p in 1:length(spl))
#     {
#       if (parNum[p]==1)
#       {
#         Pars$fixed[Pars$label%in%spl[[p]]] <- TRUE
#       } else if (parNum[p]==2)
#       {
#         Pars$par[Pars$label%in%spl[[p]]] <- parIt
#         parIt <- parIt + 1
#       } else warning("Error in computation of equality constraints.")
#     }
#   }
#   
#   if (max(Pars$par) < nrow(Pars))
#   {
#     Pars$par[Pars$par==0] <- max(Pars$par)+(1:sum(Pars$par==0))
#   }
#   
#   # Extract parameter estimates:
#   Pars$est[object$ram[,4]!=0] <- object$coef[object$ram[,4]]
#   
#   # Switch sides in regression:
#   Pars[c("lhs","rhs")][semreprObject$type=="regression",] <- Pars[c("rhs","lhs")][semreprObject$type=="regression",]
#   
#   # Set edges:
#   Pars$edge[semreprObject$type=="regression"] <- "~>"
#   Pars$edge[semreprObject$type=="latent"] <- "->"
#   Pars$edge[semreprObject$type=="covariance"] <- "<->"
#   Pars$edge[semreprObject$type=="intercept"] <- "int"
#   
#   # Variable dataframe: 
#   Vars <- data.frame(
#     name = sumObject$variables$details$Variable,
#     manifest = sumObject$variables$details$Type == "Manifest",
#     exogenous = NA,
#     stringsAsFactors=FALSE)
#   
#   # If all are latent, make guess at which are latent:
#   if (all(!Vars$manifest))
#   {
#     for (i in 1:nrow(Vars))
#     {
#       Vars$manifest[i] <- !any(semreprObject$type[semreprObject$lhs==Vars$name[i]]=="latent")
#     }
#   }
#   
#   semModel <- new("semPlotModel")
#   semModel@Pars <- Pars
#   semModel@Vars <- Vars
#   semModel@Computed <- FALSE
#   semModel@Original <- list()
#   semModel@ObsCovs <- list()
#   semModel@ImpCovs <- list()
#   
#   return(semModel)
}

