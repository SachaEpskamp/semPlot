### Path diagrams ###
# 
# setMethod("semPaths.S4",signature("lavaan"),function(object,...){
#   invisible(semPaths(semPlotModel(object),...))
# })
# 


## EXTRACT MODEL ###
setMethod("semPlotModel_S4",signature("lavaan"),function(object){

  if (is(object,"blavaan")) class(object) <- 'lavaan'
  if (!is(object,"lavaan")) stop("Input must me a 'lavaan' object")

  
  # Extract parameter estimates:
  pars <- parameterEstimates(object,standardized=TRUE)
  list <- inspect(object,"list")
  
  # Remove mean structure (TEMP SOLUTION)
  # meanstructure <- pars$op=="~1"
  # pars <- pars[!meanstructure,]
  
  # Extract variable and factor names:
  # varNames <- fit@Model@dimNames$lambda[[1]]
  # factNames <- fit@Model@dimNames$lambda[[2]]
#   Lambda <- inspect(object,"coef")$lambda
  varNames <- lavaanNames(object, type="ov")
  factNames <- lavaanNames(object, type="lv")
#   rm(Lambda)
  
  factNames <- factNames[!factNames%in%varNames]
  
  # Extract number of variables and factors
  n <- length(varNames)
  k <- length(factNames)
  
  # Extract parameter names:
  if (is.null(pars$label)) pars$label <- rep("",nrow(pars))
  
  semModel <- new("semPlotModel")

  if (is.null(pars$group)) pars$group <- ""

  # Create edges dataframe
  semModel@Pars <- data.frame(
    label = pars$label,
    lhs = ifelse(pars$op=="~"|pars$op=="~1",pars$rhs,pars$lhs),
    edge = "--",
    rhs = ifelse(pars$op=="~"|pars$op=="~1",pars$lhs,pars$rhs),
    est = pars$est,
    std = pars$std.all,
    group = pars$group,
    fixed = list$free[list$op!="=="]==0,
    par = list$free[list$op!="=="],
    stringsAsFactors=FALSE)


  semModel@Pars$edge[pars$op=="~~"] <- "<->"  
  semModel@Pars$edge[pars$op=="~*~"] <- "<->"  
  semModel@Pars$edge[pars$op=="~"] <- "~>"
  semModel@Pars$edge[pars$op=="=~"] <- "->"
  semModel@Pars$edge[pars$op=="~1"] <- "int"
  semModel@Pars$edge[grepl("\\|",pars$op)] <- "|"
  
  # Move thresholds to Thresholds slot:
  semModel@Thresholds <- semModel@Pars[grepl("\\|",semModel@Pars$edge),-(3:4)]
  
  # Remove constraints and weird stuff:
  semModel@Pars  <- semModel@Pars[!pars$op %in% c('<', '>',':=','<','>','==','|'),]
  
  # Remove thresholds from Pars:
#   semModel@Pars <- semModel@Pars[!grepl("\\|",semModel@Pars$edge),]
  
  semModel@Vars <- data.frame(
    name = c(varNames,factNames),
    manifest = c(varNames,factNames)%in%varNames,
    exogenous = NA,
    stringsAsFactors=FALSE)
    
  # res.cov <- lavTech(object, "sampstat")$res.cov
  # lavTech(object, "sampstat")$cov
  # if (!is.null(res.cov) && !length(res.cov) == 0){
      # if (!is.null(res.cov[[1]])){
      #   semModel@ObsCovs <- object@SampleStats@res.cov    
      # } else {
      #   semModel@ObsCovs <- object@SampleStats@cov
      # }    
  # } else {
  #   semModel@ObsCovs <- list(matrix(NA,
  #          length(varNames),length(varNames)))
  # } 
  
  # Use add.labels=TRUE so lavTech returns named matrices (handles multigroup with different vars)
  if (lavInspect(object, "options")$conditional.x){
    semModel@ObsCovs <- lapply(lavTech(object, "sampstat", add.labels = TRUE),"[[","res.cov")
  } else {
    semModel@ObsCovs <- lapply(lavTech(object, "sampstat", add.labels = TRUE),"[[","cov")
  }
  names(semModel@ObsCovs) <- lavInspect(object, "group.label")

  if (lavInspect(object, "options")$conditional.x){
    semModel@ImpCovs <- lapply(lavTech(object, "implied", add.labels = TRUE), "[[", "res.cov")
  } else {
    semModel@ImpCovs <- lapply(lavTech(object, "implied", add.labels = TRUE), "[[", "cov")
  }
  names(semModel@ImpCovs) <- lavInspect(object, "group.label") # object@Data@group.label
  
  semModel@Computed <- TRUE
  
  semModel@Original <- list(object)
  
  return(semModel)
})



