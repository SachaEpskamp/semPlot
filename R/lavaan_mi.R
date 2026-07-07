# Importer for lavaan.mi objects: lavaan models pooled over multiple
# imputations (package 'lavaan.mi', the successor of semTools::runMI).
# The object is S4 (contains lavaanList) but the class is only defined when
# lavaan.mi is installed, so dispatch happens by class name from
# semPlotModel() rather than via an S4 method.

semPlotModel_lavaan_mi <- function(object, ...)
{
  if (!requireNamespace("lavaan.mi", quietly = TRUE))
  {
    stop("The 'lavaan.mi' package is required to import pooled multiple-imputation models. Install it with install.packages('lavaan.mi').")
  }

  # Pooled (Rubin's rules) estimates and standardized solution:
  pars <- lavaan.mi::parameterEstimates.mi(object)
  std <- tryCatch(lavaan.mi::standardizedSolution.mi(object), error = function(e) NULL)

  # Model skeleton (free/fixed, parameter numbers) from the parameter table:
  pt <- parTable(object)
  pt <- pt[!pt$op %in% c("==","<",">"), , drop = FALSE]

  # Match rows by key, never positionally:
  keyCols <- intersect(c("lhs","op","rhs","group","block"),
                       intersect(names(pars), names(pt)))
  keyOf <- function(d) do.call(paste, c(d[keyCols], sep = "\r"))
  ptIdx <- match(keyOf(pars), keyOf(pt))
  if (any(is.na(ptIdx)))
  {
    stop("semPlot: could not align parameterEstimates.mi() with parTable(). Please report this issue, including your lavaan.mi version.")
  }
  stdVals <- rep(NA_real_, nrow(pars))
  if (!is.null(std))
  {
    stdCols <- intersect(keyCols, names(std))
    keyOf2 <- function(d) do.call(paste, c(d[stdCols], sep = "\r"))
    stdIdx <- match(keyOf2(pars), keyOf2(std))
    stdVals <- std$est.std[stdIdx]
  }

  varNames <- lavNames(object, type = "ov")
  factNames <- lavNames(object, type = "lv")
  factNames <- factNames[!factNames %in% varNames]

  if (is.null(pars$label)) pars$label <- rep("", nrow(pars))
  if (is.null(pars$group)) pars$group <- ""

  semModel <- new("semPlotModel")

  semModel@Pars <- data.frame(
    label = pars$label,
    lhs = ifelse(pars$op %in% c("~","~1","<~"), pars$rhs, pars$lhs),
    edge = "--",
    rhs = ifelse(pars$op %in% c("~","~1","<~"), pars$lhs, pars$rhs),
    est = pars$est,
    std = stdVals,
    group = pars$group,
    fixed = pt$free[ptIdx] == 0,
    par = pt$free[ptIdx],
    stringsAsFactors = FALSE)

  semModel@Pars$edge[pars$op == "~~"] <- "<->"
  semModel@Pars$edge[pars$op == "~*~"] <- "<->"
  semModel@Pars$edge[pars$op == "~"] <- "~>"
  semModel@Pars$edge[pars$op == "=~"] <- "->"
  semModel@Pars$edge[pars$op == "<~"] <- "~>"
  semModel@Pars$edge[pars$op == "~1"] <- "int"
  semModel@Pars$edge[grepl("\\|", pars$op)] <- "|"
  semModel@Pars$composite <- pars$op == "<~"

  semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge),
                                       !names(semModel@Pars) %in% c("edge","rhs")]
  semModel@Pars <- semModel@Pars[!pars$op %in% c("<",">",":=","==","|"), ]

  semModel@Vars <- data.frame(
    name = c(varNames, factNames),
    manifest = c(varNames, factNames) %in% varNames,
    exogenous = NA,
    stringsAsFactors = FALSE)

  # Pooled sample statistics are not uniquely defined across imputations:
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()

  semModel@Computed <- TRUE
  semModel@Original <- list(object)

  return(semModel)
}
