# Importer for psych-package factor-analytic objects.
#
# psych::fa() returns an object of class c("psych","fa"); its rotated,
# standardized loadings are in $loadings (a 'loadings' matrix), factor
# correlations in $Phi (present only for oblique rotations) and item
# uniquenesses in $uniquenesses. Because fa() estimates a proper factor
# model, loadings are shown as DIRECTED edges (latent -> indicator),
# unlike the exploratory loadings importer which uses undirected "--".
#
# psych::omega() returns c("psych","omega"); the Schmid-Leiman bifactor
# solution lives in $schmid$sl, a matrix with a general factor column "g",
# group-factor columns "F1*","F2*",... (trailing "*"), plus h2/u2/p2/com.
# All factors are orthogonal, so no covariance edges are drawn. Pair the
# result with semPaths(..., bifactor = "g", layout = "tree2").
#
# psych::principal() is delegated to the existing principal method.
#
# psych is not a hard dependency: this importer only reads fields of the
# supplied object, so no psych functions are called.

semPlotModel.psych <- function(object, cut = 0, ...)
{
  if (inherits(object, "principal"))
  {
    return(semPlotModel.principal(object, ...))
  }

  if (inherits(object, "fa"))
  {
    L <- unclass(object$loadings)
    items <- rownames(L)
    facs <- colnames(L)

    # Blocks of parameters; par index runs across all of them:
    lhs <- character(0)
    rhs <- character(0)
    edge <- character(0)
    est <- numeric(0)

    # 1. Loadings (latent -> indicator), keeping those above the cutoff:
    for (f in facs)
    {
      keep <- abs(L[, f]) >= cut
      lhs <- c(lhs, rep(f, sum(keep)))
      rhs <- c(rhs, items[keep])
      edge <- c(edge, rep("->", sum(keep)))
      est <- c(est, L[keep, f])
    }

    # 2. Factor correlations (oblique rotations only), upper triangle:
    Phi <- object$Phi
    if (is.matrix(Phi))
    {
      ind <- which(upper.tri(Phi), arr.ind = TRUE)
      if (nrow(ind) > 0)
      {
        lhs <- c(lhs, facs[ind[, 1]])
        rhs <- c(rhs, facs[ind[, 2]])
        edge <- c(edge, rep("<->", nrow(ind)))
        est <- c(est, Phi[ind])
      }
    }

    # 3. Residual variances (uniquenesses) as self-loops:
    uniq <- object$uniquenesses
    lhs <- c(lhs, items)
    rhs <- c(rhs, items)
    edge <- c(edge, rep("<->", length(items)))
    est <- c(est, as.numeric(uniq[items]))

    Pars <- data.frame(
      label = "",
      lhs = lhs,
      edge = edge,
      rhs = rhs,
      est = est,
      std = est,
      group = "",
      fixed = FALSE,
      par = seq_along(est),
      stringsAsFactors = FALSE)

    Vars <- data.frame(
      name = c(items, facs),
      manifest = c(rep(TRUE, length(items)), rep(FALSE, length(facs))),
      exogenous = NA,
      stringsAsFactors = FALSE)

    ObsCovs <- if (is.matrix(object$r)) list(object$r) else list()
    ImpCovs <- if (is.matrix(object$model) &&
                   all(dim(object$model) == length(items))) list(object$model) else list()

    semModel <- new("semPlotModel")
    semModel@Pars <- Pars
    semModel@Vars <- Vars
    semModel@Computed <- TRUE
    semModel@Original <- list(object)
    semModel@ObsCovs <- ObsCovs
    semModel@ImpCovs <- ImpCovs

    return(semModel)
  }

  if (inherits(object, "omega"))
  {
    sl <- object$schmid$sl
    items <- rownames(sl)
    slCols <- colnames(sl)

    # General factor "g" and group factors (colnames ending in "*"):
    grpCols <- slCols[grepl("\\*$", slCols)]
    grpNames <- sub("\\*$", "", grpCols)
    facs <- c("g", grpNames)

    lhs <- character(0)
    rhs <- character(0)
    edge <- character(0)
    est <- numeric(0)

    # 1. General-factor loadings:
    gLoad <- sl[, "g"]
    keep <- abs(gLoad) >= cut
    lhs <- c(lhs, rep("g", sum(keep)))
    rhs <- c(rhs, items[keep])
    edge <- c(edge, rep("->", sum(keep)))
    est <- c(est, gLoad[keep])

    # 2. Group-factor loadings:
    for (j in seq_along(grpCols))
    {
      col <- sl[, grpCols[j]]
      keep <- abs(col) >= cut
      lhs <- c(lhs, rep(grpNames[j], sum(keep)))
      rhs <- c(rhs, items[keep])
      edge <- c(edge, rep("->", sum(keep)))
      est <- c(est, col[keep])
    }

    # 3. Residual variances from the "u2" column, as self-loops:
    u2 <- sl[, "u2"]
    lhs <- c(lhs, items)
    rhs <- c(rhs, items)
    edge <- c(edge, rep("<->", length(items)))
    est <- c(est, as.numeric(u2))

    Pars <- data.frame(
      label = "",
      lhs = lhs,
      edge = edge,
      rhs = rhs,
      est = est,
      std = est,
      group = "",
      fixed = FALSE,
      par = seq_along(est),
      stringsAsFactors = FALSE)

    # Factors: "g" first among the latents (for bifactor = "g" layouts):
    Vars <- data.frame(
      name = c(items, facs),
      manifest = c(rep(TRUE, length(items)), rep(FALSE, length(facs))),
      exogenous = NA,
      stringsAsFactors = FALSE)

    semModel <- new("semPlotModel")
    semModel@Pars <- Pars
    semModel@Vars <- Vars
    semModel@Computed <- TRUE
    semModel@Original <- list(object)
    semModel@ObsCovs <- list()
    semModel@ImpCovs <- list()

    return(semModel)
  }

  stop("psych objects of class '", paste(class(object), collapse = "/"),
       "' are not supported by semPlot (supported: fa, omega, principal).")
}
