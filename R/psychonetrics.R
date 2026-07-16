# Importer for psychonetrics models (github.com/SachaEpskamp/psychonetrics).
#
# Supports the 'lvm' and 'varcov' frameworks, including every subtype of
# their (co)variance structures: "cov" (covariances), "ggm" (undirected
# network Omega + scaling Delta), "chol" (Cholesky) and "prec" (precision).
# GGM blocks are displayed as undirected ("--") edges; the Delta scaling
# parameters are, by default, replaced by the implied variances
# (diag of Sigma = Delta (I - Omega)^-1 Delta), which are directly
# interpretable. "chol"/"prec" blocks are displayed as implied covariances.
#
# psychonetrics is a Suggests dependency; this function is called from
# semPlotModel() via class dispatch and must swallow unused arguments
# passed through modelOpts (e.g. mplusStd) in '...'.

semPlotModel_psychonetrics <- function(object, delta = c("variance","ignore","show"), ...)
{
  if (!requireNamespace("psychonetrics", quietly = TRUE))
  {
    stop("The 'psychonetrics' package is required to import psychonetrics models. Install it with install.packages('psychonetrics').")
  }
  delta <- match.arg(delta)

  framework <- object@model
  if (!framework %in% c("lvm","varcov"))
  {
    stop("psychonetrics framework '", framework, "' is not (yet) supported by semPlot. ",
         "Currently supported: 'lvm' and 'varcov' (including the ggm/precision wrappers). ",
         "For temporal / multi-level frameworks, consider plotting the relevant matrix directly, ",
         "e.g. qgraph::qgraph(psychonetrics::getmatrix(object, \"omega_zeta_within\")).")
  }

  pars <- as.data.frame(object@parameters)
  # Select columns by name; newer psychonetrics versions append many extra columns:
  needed <- c("var1","op","var2","est","std","se","matrix","row","col","par","group","group_id","fixed")
  missingCols <- needed[!needed %in% names(pars)]
  if (length(missingCols) > 0)
  {
    stop("Unexpected psychonetrics parameter table: missing column(s) ", paste(missingCols, collapse = ", "),
         ". Please report this together with your psychonetrics version.")
  }
  pars <- pars[needed]

  # Group labels in group_id order ('fullsample' for single-group models):
  grpLabels <- unique(pars$group[order(pars$group_id)])
  singleGroup <- length(grpLabels) == 1

  # Variable names (order matches psychonetrics' matrix rows). Latent names are
  # taken in 'lambda' column order, so that latNames[j] is the latent belonging
  # to column j of lambda/beta and can be used to index the standardizing SDs
  # below. This runs before the structural-zero rows are dropped, so every
  # column is still represented.
  manNames <- object@sample@variables$label
  lambdaRows <- pars$matrix == "lambda"
  latNames <- unique(pars$var2[lambdaRows][order(pars$col[lambdaRows])])
  latNames <- latNames[!is.na(latNames)]

  # Thresholds not (yet) supported:
  if (any(pars$matrix == "tau"))
  {
    message("Thresholds ('tau') of psychonetrics models are not yet supported by semPlot and are not shown.")
    pars <- pars[pars$matrix != "tau", ]
  }

  # How each block is parameterized is a property of the model, so it is read
  # from the full parameter table, before the structural-zero drop below. An
  # empty network (every omega fixed to zero, as for a default residual = "ggm")
  # would otherwise disappear in that drop, leaving the block looking
  # non-GGM and its delta rows unhandled.
  allMatrices <- unique(pars$matrix)

  # Drop structurally-zero placeholder rows (e.g. excluded loadings):
  pars <- pars[!(pars$fixed & pars$par == 0 & pars$est == 0), ]

  # Normalize a getmatrix() result to a named per-group list:
  getBlockMatrix <- function(name)
  {
    m <- tryCatch(psychonetrics::getmatrix(object, name), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    if (!is.list(m)) m <- list(m)
    if (is.null(names(m)) || any(names(m) == "")) names(m) <- grpLabels[seq_along(m)]
    m
  }

  # (Co)variance blocks and the node names their matrices refer to:
  blockSuffix <- if (framework == "lvm") c("_zeta","_epsilon") else ""
  blockNodes <- function(sfx)
  {
    if (sfx == "_zeta") latNames else manNames
  }

  # ---- Standardizing standard deviations --------------------------------
  # psychonetrics leaves the 'std' column of @parameters empty for these
  # models, so the standardized solution is computed here from the
  # model-implied matrices:
  #
  #   lvm:    Var(eta) = (I - Beta)^-1 Sigma_zeta (I - Beta)^-T
  #           Sigma    = Lambda Var(eta) Lambda' + Sigma_epsilon
  #   varcov: Sigma is the implied covariance matrix itself.
  #
  # For a GGM-parameterized block psychonetrics forms
  # Sigma_* = Delta (I - Omega)^-1 Delta, so 'sigma_*' already accounts for the
  # network parameterization and needs no special-casing here.
  impSigma <- getBlockMatrix("sigma")
  impSigmaZeta <- if (framework == "lvm") getBlockMatrix("sigma_zeta") else NULL
  impBeta <- if (framework == "lvm") getBlockMatrix("beta") else NULL

  # A non-positive implied variance means the solution is improper (e.g. a
  # Heywood case); standardizing against it is undefined, so those parameters
  # keep std = NA and the affected variables are reported once below.
  improper <- character(0)
  sdFromVar <- function(v, nms, group)
  {
    bad <- !is.na(v) & v <= 0
    if (any(bad))
    {
      improper <<- union(improper, if (singleGroup) nms[bad] else paste0(nms[bad], " (", group, ")"))
    }
    sd <- rep(NA_real_, length(v))
    sd[!bad] <- sqrt(v[!bad])
    sd
  }

  sdOf <- list()
  for (g in grpLabels)
  {
    sdMan <- rep(NA_real_, length(manNames))
    if (!is.null(impSigma) && !is.null(impSigma[[g]]))
    {
      sdMan <- sdFromVar(diag(as.matrix(impSigma[[g]])), manNames, g)
    }
    sdLat <- rep(NA_real_, length(latNames))
    if (length(latNames) > 0 && !is.null(impSigmaZeta) && !is.null(impSigmaZeta[[g]]))
    {
      Sz <- as.matrix(impSigmaZeta[[g]])
      B <- if (!is.null(impBeta) && !is.null(impBeta[[g]])) as.matrix(impBeta[[g]]) else matrix(0, nrow(Sz), ncol(Sz))
      varEta <- tryCatch({
        IB <- solve(diag(nrow(B)) - B)
        IB %*% Sz %*% t(IB)
      }, error = function(e) NULL)
      if (!is.null(varEta)) sdLat <- sdFromVar(diag(varEta), latNames, g)
    }
    sdOf[[g]] <- list(man = sdMan, lat = sdLat)
  }

  if (length(improper) > 0)
  {
    warning("The model-implied variance is not positive for: ", paste(improper, collapse = ", "),
            ". This indicates an improper solution; standardized parameters involving ",
            "these variables are NA.")
  }

  # SDs belonging to a (co)variance block, for the derived rows added below:
  blockSD <- function(sfx, g)
  {
    if (sfx == "_zeta") sdOf[[g]]$lat else sdOf[[g]]$man
  }

  # Standardize one parameter from its matrix and row/col indices.
  stdValue <- function(mat, row, col, est, group)
  {
    sds <- sdOf[[group]]
    if (is.null(sds)) return(NA_real_)
    at <- function(x, i)
    {
      i <- suppressWarnings(as.integer(i))
      if (!is.na(i) && i >= 1 && i <= length(x)) x[i] else NA_real_
    }
    lat <- sds$lat
    man <- sds$man
    switch(as.character(mat),
      # Omega is by construction the partial correlation matrix and is thus
      # already standardized: with Kappa = Delta^-1 (I - Omega) Delta^-1,
      # -Kappa_ij / sqrt(Kappa_ii Kappa_jj) = Omega_ij for any Delta.
      "omega_zeta"    = est,
      "omega_epsilon" = est,
      "omega"         = est,
      "lambda"        = est * at(lat, col) / at(man, row),
      "beta"          = est * at(lat, col) / at(lat, row),
      "nu"            = est / at(man, row),
      "mu"            = est / at(man, row),
      "nu_eta"        = est / at(lat, row),
      "sigma_zeta"    = est / (at(lat, row) * at(lat, col)),
      "sigma_epsilon" = est / (at(man, row) * at(man, col)),
      "sigma"         = est / (at(man, row) * at(man, col)),
      # Delta is a scaling matrix; the corresponding scaling of the
      # standardized (correlation) metric is Delta / sd.
      "delta_zeta"    = est / at(lat, row),
      "delta_epsilon" = est / at(man, row),
      "delta"         = est / at(man, row),
      NA_real_)
  }

  # Rows that only encode a parameterization (not path-diagram parameters):
  parameterizationRows <- function(sfx)
  {
    pars$matrix %in% paste0(c("delta","lowertri","kappa"), sfx)
  }

  derivedPars <- list()
  addDerived <- function(nodes1, nodes2, ests, group, stds = NA_real_)
  {
    derivedPars[[length(derivedPars) + 1]] <<- data.frame(
      label = "",
      lhs = nodes1,
      edge = "<->",
      rhs = nodes2,
      est = ests,
      std = stds,
      group = group,
      fixed = FALSE,
      par = 0,
      stringsAsFactors = FALSE)
  }

  for (sfx in blockSuffix)
  {
    hasGGM <- paste0("omega", sfx) %in% allMatrices
    hasChol <- paste0("lowertri", sfx) %in% allMatrices
    hasPrec <- paste0("kappa", sfx) %in% allMatrices
    nodes <- blockNodes(sfx)

    if (hasGGM)
    {
      # Omega rows (op '--') stay as undirected edges. Delta rows:
      deltaRows <- pars$matrix == paste0("delta", sfx)
      if (delta == "show")
      {
        # Show raw scaling parameters as self-loops:
        pars$op[deltaRows] <- "~~"
        pars$var2[deltaRows] <- pars$var1[deltaRows]
      } else
      {
        pars <- pars[!deltaRows, ]
        if (delta == "variance")
        {
          sig <- getBlockMatrix(paste0("sigma", sfx))
          if (is.null(sig))
          {
            warning("Could not obtain the implied covariance matrix 'sigma", sfx,
                    "'; variances of this block are not shown.")
          } else
          {
            for (g in grpLabels)
            {
              vars <- diag(as.matrix(sig[[g]]))
              sdv <- blockSD(sfx, g)
              addDerived(nodes, nodes, vars, g, vars / (sdv * sdv))
            }
          }
        }
      }
    }

    if (hasChol || hasPrec)
    {
      # Display Cholesky / precision parameterized blocks as implied covariances:
      pars <- pars[!pars$matrix %in% paste0(c("lowertri","kappa"), sfx), ]
      sig <- getBlockMatrix(paste0("sigma", sfx))
      if (is.null(sig))
      {
        warning("Could not obtain the implied covariance matrix 'sigma", sfx,
                "'; this block is not shown.")
      } else
      {
        blockName <- switch(sfx, "_zeta" = "latent", "_epsilon" = "residual", "observed")
        message("The ", blockName, " block is parameterized as ",
                ifelse(hasChol, "'chol'", "'prec'"),
                "; displaying implied covariances.")
        for (g in grpLabels)
        {
          S <- as.matrix(sig[[g]])
          ind <- which(upper.tri(S, diag = TRUE), arr.ind = TRUE)
          ests <- S[ind]
          nonZero <- abs(ests) > sqrt(.Machine$double.eps)
          sdv <- blockSD(sfx, g)
          addDerived(nodes[ind[nonZero, 1]], nodes[ind[nonZero, 2]], ests[nonZero], g,
                     ests[nonZero] / (sdv[ind[nonZero, 1]] * sdv[ind[nonZero, 2]]))
        }
      }
    }
  }

  # Fill in the standardized solution, keeping any value psychonetrics itself
  # supplied:
  if (nrow(pars) > 0)
  {
    computedStd <- mapply(stdValue, pars$matrix, pars$row, pars$col, pars$est, pars$group,
                          SIMPLIFY = TRUE, USE.NAMES = FALSE)
    computedStd <- as.numeric(computedStd)
    computedStd[!is.finite(computedStd)] <- NA_real_
    pars$std <- ifelse(is.na(pars$std), computedStd, pars$std)
  }

  # Map psychonetrics operators to semPlot edge types:
  known <- pars$op %in% c("~1","~=","<-","~~","--")
  if (any(!known))
  {
    warning("Dropped ", sum(!known), " parameter(s) with unsupported operator(s): ",
            paste(unique(pars$op[!known]), collapse = ", "))
    pars <- pars[known, ]
  }

  Pars <- data.frame(
    label = "",
    lhs = ifelse(pars$op == "~1", "",
          ifelse(pars$op %in% c("~=","<-"), pars$var2, pars$var1)),
    edge = ifelse(pars$op == "~1", "int",
           ifelse(pars$op == "~=", "->",
           ifelse(pars$op == "<-", "~>",
           ifelse(pars$op == "~~", "<->", "--")))),
    rhs = ifelse(pars$op %in% c("~1","~=","<-"), pars$var1, pars$var2),
    est = pars$est,
    std = pars$std,
    group = pars$group,
    fixed = pars$fixed,
    par = pars$par,
    stringsAsFactors = FALSE)

  if (length(derivedPars) > 0)
  {
    Pars <- rbind(Pars, do.call(rbind, derivedPars))
  }

  # Single-group models get an empty group label (no per-group titles):
  if (singleGroup)
  {
    Pars$group <- ""
  }

  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- data.frame(
    name = c(manNames, latNames),
    manifest = c(rep(TRUE, length(manNames)), rep(FALSE, length(latNames))),
    exogenous = NA,
    stringsAsFactors = FALSE)
  semModel@Thresholds <- data.frame()
  semModel@Computed <- isTRUE(object@computed)

  ObsCovs <- lapply(object@sample@covs, as.matrix)
  ImpCovs <- getBlockMatrix("sigma")
  if (is.null(ImpCovs)) ImpCovs <- list()
  covNames <- if (singleGroup) "" else grpLabels
  if (length(ObsCovs) == length(covNames)) names(ObsCovs) <- covNames
  if (length(ImpCovs) == length(covNames)) names(ImpCovs) <- covNames
  semModel@ObsCovs <- ObsCovs
  semModel@ImpCovs <- ImpCovs

  semModel@Original <- list(object)

  return(semModel)
}
