# Importers for composite-based (PLS) SEM: seminr and cSEM.
# PLS estimates are standardized; est and std carry the same values.
# Display conventions: reflective(-ly displayed) constructs get directed
# loadings construct -> item; composite (mode B / formative) constructs get
# weight arrows item -> construct.

semPlotModel.seminr_model <- function(object, ...)
{
  mm <- object$mmMatrix
  if (is.null(mm) || is.null(object$path_coef))
  {
    stop("This seminr model does not contain the expected mmMatrix/path_coef elements.")
  }

  constructs <- unique(mm[, "construct"])
  # Only measurement rows whose item actually appears in the loading matrix
  # (guards against special construct types, e.g. interactions):
  mm <- mm[mm[, "measurement"] %in% rownames(object$outer_loadings), , drop = FALSE]
  items <- unique(mm[, "measurement"])

  Pars <- list()
  addPar <- function(lhs, edge, rhs, est)
  {
    Pars[[length(Pars) + 1]] <<- data.frame(
      label = "", lhs = lhs, edge = edge, rhs = rhs,
      est = est, std = est, group = "", fixed = FALSE, par = 0,
      stringsAsFactors = FALSE)
  }

  # Measurement model: mode B ("B") constructs get formative weight arrows,
  # everything else (mode A composites "A", reflective "C") gets loadings:
  for (i in seq_len(nrow(mm)))
  {
    cons <- mm[i, "construct"]
    item <- mm[i, "measurement"]
    if (unname(mm[i, "type"]) == "B")
    {
      addPar(item, "~>", cons, object$outer_weights[item, cons])
    } else
    {
      addPar(cons, "->", item, object$outer_loadings[item, cons])
    }
  }

  # Structural model: path_coef is [from, to]:
  pc <- object$path_coef
  idx <- which(pc != 0, arr.ind = TRUE)
  for (k in seq_len(nrow(idx)))
  {
    addPar(rownames(pc)[idx[k, 1]], "~>", colnames(pc)[idx[k, 2]], pc[idx[k, 1], idx[k, 2]])
  }

  Pars <- do.call(rbind, Pars)
  Pars$par <- seq_len(nrow(Pars))

  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- data.frame(
    name = c(items, constructs),
    manifest = c(rep(TRUE, length(items)), rep(FALSE, length(constructs))),
    exogenous = NA,
    stringsAsFactors = FALSE)
  semModel@Thresholds <- data.frame()
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()

  return(semModel)
}

semPlotModel.cSEMResults <- function(object, ...)
{
  if (!requireNamespace("cSEM", quietly = TRUE))
  {
    stop("The 'cSEM' package is required to import cSEM models. Install it with install.packages('cSEM').")
  }
  if (inherits(object, "cSEMResults_multi"))
  {
    stop("Multi-group/multi-dataset cSEM results are not yet supported; import a single group, e.g. semPlotModel(object[[1]]).")
  }

  sm <- cSEM::summarize(object)
  est <- sm$Estimates
  cType <- object$Information$Model$construct_type

  # Parse "lhs op rhs" name strings from the tidy estimate tables:
  parseName <- function(name, op)
  {
    parts <- strsplit(name, paste0(" ", op, " "), fixed = TRUE)
    do.call(rbind, lapply(parts, function(x) data.frame(lhs = x[1], rhs = x[2], stringsAsFactors = FALSE)))
  }

  Pars <- list()
  addPars <- function(lhs, edge, rhs, ests)
  {
    Pars[[length(Pars) + 1]] <<- data.frame(
      label = "", lhs = lhs, edge = edge, rhs = rhs,
      est = ests, std = ests, group = "", fixed = FALSE, par = 0,
      stringsAsFactors = FALSE)
  }

  # Structural paths ("lhs ~ rhs" means rhs predicts lhs):
  pe <- est$Path_estimates
  if (!is.null(pe) && nrow(pe) > 0)
  {
    nm <- parseName(pe$Name, "~")
    addPars(nm$rhs, "~>", nm$lhs, pe$Estimate)
  }

  # Measurement: loadings for common factors, weights for composites:
  le <- est$Loading_estimates
  if (!is.null(le) && nrow(le) > 0)
  {
    nm <- parseName(le$Name, "=~")
    keep <- cType[nm$lhs] != "Composite"
    if (any(keep)) addPars(nm$lhs[keep], "->", nm$rhs[keep], le$Estimate[keep])
  }
  we <- est$Weight_estimates
  if (!is.null(we) && nrow(we) > 0)
  {
    nm <- parseName(we$Name, "<~")
    keep <- cType[nm$lhs] == "Composite"
    if (any(keep)) addPars(nm$rhs[keep], "~>", nm$lhs[keep], we$Estimate[keep])
  }

  Pars <- do.call(rbind, Pars)
  Pars$par <- seq_len(nrow(Pars))

  constructs <- names(cType)
  items <- unique(Pars$lhs[!Pars$lhs %in% constructs])
  items <- unique(c(items, Pars$rhs[!Pars$rhs %in% constructs]))

  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- data.frame(
    name = c(items, constructs),
    manifest = c(rep(TRUE, length(items)), rep(FALSE, length(constructs))),
    exogenous = NA,
    stringsAsFactors = FALSE)
  semModel@Thresholds <- data.frame()
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()

  return(semModel)
}
