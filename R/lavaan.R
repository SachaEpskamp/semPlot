### Path diagrams ###
# 
# setMethod("semPaths.S4",signature("lavaan"),function(object,...){
#   invisible(semPaths(semPlotModel(object),...))
# })
# 


## EXTRACT MODEL ###
setMethod("semPlotModel_S4",signature("lavaan"),function(object){

  # blavaan objects extend lavaan and dispatch here via inheritance; all
  # accessors used below (parameterEstimates, parTable, lavNames, lavInspect,
  # lavTech) work on them natively, so no class manipulation is needed (an
  # older version overwrote class(object), discarding Bayesian information).
  if (!is(object,"lavaan")) stop("Input must be a 'lavaan' object")

  
  # Extract parameter estimates and the parameter table:
  pars <- parameterEstimates(object,standardized=TRUE)
  pt <- parTable(object)
  pt <- pt[!pt$op %in% c("==","<",">"), , drop = FALSE]

  # Align pars with the parameter table by matching keys rather than
  # positionally (robust against row-set differences between lavaan
  # versions). Key columns are those present in both tables:
  keyCols <- intersect(c("lhs","op","rhs","group","block"),
                       intersect(names(pars), names(pt)))
  keyOf <- function(d) do.call(paste, c(d[keyCols], sep = "\r"))
  ptIdx <- match(keyOf(pars), keyOf(pt))
  if (any(is.na(ptIdx)))
  {
    stop("semPlot: could not align parameterEstimates() with parTable(). Please report this issue, including your lavaan version.")
  }

  # Remove mean structure (TEMP SOLUTION)
  # meanstructure <- pars$op=="~1"
  # pars <- pars[!meanstructure,]
  
  # Extract variable and factor names:
  # varNames <- fit@Model@dimNames$lambda[[1]]
  # factNames <- fit@Model@dimNames$lambda[[2]]
#   Lambda <- inspect(object,"coef")$lambda
  varNames <- lavNames(object, type="ov")
  factNames <- lavNames(object, type="lv")
#   rm(Lambda)
  
  factNames <- factNames[!factNames%in%varNames]
  
  # Extract number of variables and factors
  n <- length(varNames)
  k <- length(factNames)
  
  # Extract parameter names:
  if (is.null(pars$label)) pars$label <- rep("",nrow(pars))
  
  semModel <- new("semPlotModel")

  if (is.null(pars$group)) pars$group <- ""

  # Create edges dataframe. Regressions ('~'), intercepts ('~1') and
  # composite/formative loadings ('<~') store the arrow origin in rhs:
  semModel@Pars <- data.frame(
    label = pars$label,
    lhs = ifelse(pars$op%in%c("~","~1","<~"),pars$rhs,pars$lhs),
    edge = "--",
    rhs = ifelse(pars$op%in%c("~","~1","<~"),pars$lhs,pars$rhs),
    est = pars$est,
    std = pars$std.all,
    group = pars$group,
    fixed = pt$free[ptIdx]==0,
    par = pt$free[ptIdx],
    stringsAsFactors=FALSE)


  semModel@Pars$edge[pars$op=="~~"] <- "<->"
  semModel@Pars$edge[pars$op=="~*~"] <- "<->"
  semModel@Pars$edge[pars$op=="~"] <- "~>"
  semModel@Pars$edge[pars$op=="=~"] <- "->"
  semModel@Pars$edge[pars$op=="<~"] <- "~>"
  semModel@Pars$edge[pars$op=="~1"] <- "int"
  semModel@Pars$edge[grepl("\\|",pars$op)] <- "|"

  # Mark composite indicator edges so semSyntax() can re-emit them as '<~':
  semModel@Pars$composite <- pars$op == "<~"

  # Multilevel models: map lavaan's level column onto the Within/Between
  # panel mechanism of semPaths (level 1 = Within, higher levels = Between):
  nLevels <- lavInspect(object, "nlevels")
  if (nLevels > 1 && !is.null(pars$level))
  {
    if (nLevels > 2) warning("Models with more than two levels are displayed as Within (level 1) versus Between (all higher levels).")
    semModel@Pars$BetweenWithin <- ifelse(pars$level == 1, "Within", "Between")
  }

  # Move thresholds to Thresholds slot:
  semModel@Thresholds <- semModel@Pars[grepl("\\|",semModel@Pars$edge), !names(semModel@Pars) %in% c("edge","rhs")]
  
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
  covName <- if (isTRUE(lavInspect(object, "options")$conditional.x)) "res.cov" else "cov"
  semModel@ObsCovs <- lapply(lavTech(object, "sampstat", add.labels = TRUE), "[[", covName)
  semModel@ImpCovs <- lapply(lavTech(object, "implied", add.labels = TRUE), "[[", covName)

  # Name the blocks: group labels, crossed with Within/Between for
  # multilevel models (lavaan blocks are ordered group-major):
  grpLab <- lavInspect(object, "group.label")
  if (length(grpLab) == 0) grpLab <- ""
  if (nLevels > 1)
  {
    lvlLab <- c("Within", rep("Between", nLevels - 1))
    blockNames <- as.vector(t(outer(grpLab, lvlLab, function(g, l) trimws(paste(g, "-", l)))))
    blockNames <- sub("^- ", "", blockNames)
  } else {
    blockNames <- grpLab
  }
  if (length(semModel@ObsCovs) == length(blockNames) && any(blockNames != "")) names(semModel@ObsCovs) <- blockNames
  if (length(semModel@ImpCovs) == length(blockNames) && any(blockNames != "")) names(semModel@ImpCovs) <- blockNames
  
  semModel@Computed <- TRUE
  
  semModel@Original <- list(object)
  
  return(semModel)
})



# lavaan's efa() returns an 'efaList': a plain (non-S4) list of complete lavaan
# fits, one per requested number of factors, named e.g. "nf2", "nf3". Import a
# selected solution by delegating to the existing lavaan (S4) importer.
semPlotModel.efaList <- function(object, which, ...)
{
  fits <- unclass(object)

  # Besides the fits (named e.g. "nf2"), recent lavaan versions store extra
  # elements (such as "loadings") in the list -- keep only the lavaan fits:
  fits <- Filter(function(x) is(x, "lavaan"), fits)
  if (length(fits) == 0) stop("This efaList contains no lavaan fits.")

  # Default to the last solution (the highest number of factors):
  if (missing(which)) which <- length(fits)

  # Validate the selection and resolve a display name:
  nms <- names(fits)
  if (is.character(which)) {
    if (!which %in% nms) {
      stop("'which' must be one of: ", paste(nms, collapse = ", "))
    }
    selName <- which
  } else {
    if (!is.numeric(which) || length(which) != 1 ||
        which < 1 || which > length(fits)) {
      stop("'which' must be a single index between 1 and ", length(fits),
           ", or one of: ", paste(nms, collapse = ", "))
    }
    selName <- if (is.null(nms) || is.na(nms[[which]])) as.character(which) else nms[[which]]
  }

  if (length(fits) > 1) {
    available <- if (is.null(nms)) seq_along(fits) else nms
    message("efaList contains ", length(fits), " solutions (",
            paste(available, collapse = ", "), "); importing '", selName,
            "'. Use 'which' to select another.")
  }

  semPlotModel(fits[[which]], ...)
}
