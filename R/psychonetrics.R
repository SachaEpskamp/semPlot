# Importer for psychonetrics models (github.com/SachaEpskamp/psychonetrics).
#
# Supports the 'lvm', 'varcov' and 'cvm' frameworks, including every subtype
# of their (co)variance structures: "cov" (covariances), "ggm" (undirected
# network Omega + scaling Delta), "cor" (correlations Rho + standard
# deviations SD), "chol" (Cholesky) and "prec" (precision).
# GGM blocks are displayed as undirected ("--") edges and correlation blocks
# as (co)variance edges holding the correlations; the scaling parameters
# (Delta or SD) are, by default, replaced by the implied variances (e.g., diag
# of Sigma = Delta (I - Omega)^-1 Delta), which are directly interpretable.
# "chol"/"prec" blocks are displayed as implied covariances.
#
# The 'cvm' framework (cvm() and cnm()) contains composites next to latent
# variables. Composites are marked in Vars$composite (drawn as hexagons) and
# their weights are '~>' edges from the indicator into the composite, marked
# in Pars$composite. The '_c' block holds the (co)variances of the composite
# indicators, which is a network of these indicators for cnm(). The (residual)
# variances and intercepts of composites are computed by psychonetrics and are
# shown as fixed parameters.
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
  if (!framework %in% c("lvm","varcov","cvm"))
  {
    stop("psychonetrics framework '", framework, "' is not (yet) supported by semPlot. ",
         "Currently supported: 'lvm', 'varcov' (including the ggm/precision wrappers) and ",
         "'cvm' (composite models from cvm() and cnm()). ",
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
  # Variable indices: 'row' and 'col' index the matrix of a parameter, which in
  # the cvm framework covers only part of the variables for the '_epsilon' and
  # '_c' blocks. 'var1_id' and 'var2_id' always index all variables (or all
  # latent variables):
  if (is.null(pars$var1_id)) pars$var1_id <- pars$row
  if (is.null(pars$var2_id)) pars$var2_id <- pars$col
  pars <- pars[c(needed, "var1_id", "var2_id")]

  # Group labels in group_id order ('fullsample' for single-group models):
  grpLabels <- unique(pars$group[order(pars$group_id)])
  singleGroup <- length(grpLabels) == 1

  # Variable names (order matches psychonetrics' matrix rows). Latent names are
  # taken in 'lambda' column order, so that latNames[j] is the latent belonging
  # to column j of lambda/beta and can be used to index the standardizing SDs
  # below. This runs before the structural-zero rows are dropped, so every
  # column is still represented. In the cvm framework lambda has a column for
  # every construct, composites included.
  manNames <- object@sample@variables$label
  lambdaRows <- pars$matrix == "lambda"
  latNames <- unique(pars$var2[lambdaRows][order(pars$col[lambdaRows])])
  latNames <- latNames[!is.na(latNames)]

  # Composites (cvm framework):
  compNames <- character(0)
  if (framework == "cvm")
  {
    compvec <- object@extramatrices$compvec
    if (!is.null(compvec) && length(compvec) == length(latNames))
    {
      compNames <- latNames[as.vector(as.matrix(compvec)) == 1]
    } else
    {
      # Constructs with at least one weight that is not a structural zero:
      wRows <- pars$matrix == "weights" & !(pars$fixed & pars$est == 0)
      compNames <- unique(pars$var2[wRows])
    }
  }

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

  # (Co)variance blocks and the node names their matrices refer to:
  blockSuffix <- switch(framework,
                        "lvm" = c("_zeta","_epsilon"),
                        "cvm" = c("_zeta","_epsilon","_c"),
                        "")
  blockMatrices <- function(sfx)
  {
    paste0(c("sigma","kappa","lowertri","omega","delta","rho","SD"), sfx)
  }
  # The node of every row/column of a block, read from the full parameter
  # table of the first group (every parameterization has a matrix with the
  # diagonal in the table, so every node is covered):
  blockNodeList <- list()
  for (sfx in blockSuffix)
  {
    r <- pars[pars$matrix %in% blockMatrices(sfx) & pars$group_id == min(pars$group_id), , drop = FALSE]
    nodes <- character(0)
    if (nrow(r) > 0)
    {
      nodes <- rep(NA_character_, max(c(r$row, r$col)))
      nodes[r$row] <- r$var1
      nodes[r$col] <- ifelse(is.na(r$var2), r$var1, r$var2)
    }
    blockNodeList[[paste0("block", sfx)]] <- nodes
  }
  blockNodes <- function(sfx) blockNodeList[[paste0("block", sfx)]]

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

  # ---- Standardizing standard deviations --------------------------------
  # psychonetrics leaves the 'std' column of @parameters empty for these
  # models, so the standardized solution is computed here from the
  # model-implied matrices:
  #
  #   lvm/cvm: Var(eta) = (I - Beta)^-1 Sigma_zeta (I - Beta)^-T
  #            Sigma    = model-implied covariance matrix of the variables
  #   varcov:  Sigma is the implied covariance matrix itself.
  #
  # For a GGM-parameterized block psychonetrics forms
  # Sigma_* = Delta (I - Omega)^-1 Delta, so 'sigma_*' already accounts for the
  # network parameterization and needs no special-casing here. In the cvm
  # framework Sigma_zeta contains the computed (residual) variances of the
  # composites, so that Var(eta) contains the variances of the composites.
  impSigma <- getBlockMatrix("sigma")
  impSigmaZeta <- if (framework %in% c("lvm","cvm")) getBlockMatrix("sigma_zeta") else NULL
  impBeta <- if (framework %in% c("lvm","cvm")) getBlockMatrix("beta") else NULL

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

  # SDs belonging to the nodes of a (co)variance block, for the derived rows
  # added below:
  blockSD <- function(sfx, g)
  {
    nodes <- blockNodes(sfx)
    if (sfx == "_zeta") sdOf[[g]]$lat[match(nodes, latNames)] else sdOf[[g]]$man[match(nodes, manNames)]
  }

  # Standardize one parameter from its matrix and the (global) indices of its
  # variables.
  stdValue <- function(mat, v1, v2, est, group)
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
    mat <- as.character(mat)
    switch(mat,
      "lambda"  = est * at(lat, v2) / at(man, v1),
      "beta"    = est * at(lat, v2) / at(lat, v1),
      # Composite weights: from the indicator (var1) into the composite (var2):
      "weights" = est * at(man, v1) / at(lat, v2),
      "nu"      = est / at(man, v1),
      "mu"      = est / at(man, v1),
      "nu_eta"  = est / at(lat, v1),
      {
        # (Co)variance blocks: 'sigma_zeta', 'omega_epsilon', 'delta_c', ...
        type <- sub("_.*$", "", mat)
        sd <- if (grepl("_zeta$", mat)) lat else man
        switch(type,
          # Omega is by construction the partial correlation matrix and is thus
          # already standardized: with Kappa = Delta^-1 (I - Omega) Delta^-1,
          # -Kappa_ij / sqrt(Kappa_ii Kappa_jj) = Omega_ij for any Delta. The
          # same holds for the correlations in Rho.
          "omega" = est,
          "rho"   = est,
          "sigma" = est / (at(sd, v1) * at(sd, v2)),
          # Delta and SD are scaling matrices; the corresponding scaling of the
          # standardized (correlation) metric is Delta / sd.
          "delta" = est / at(sd, v1),
          "SD"    = est / at(sd, v1),
          NA_real_)
      })
  }

  derivedPars <- list()
  addDerived <- function(nodes1, nodes2, ests, group, stds = NA_real_)
  {
    # The (residual) variances of composites are computed, not estimated:
    computed <- nodes1 == nodes2 & nodes1 %in% compNames
    derivedPars[[length(derivedPars) + 1]] <<- data.frame(
      label = "",
      lhs = nodes1,
      edge = "<->",
      rhs = nodes2,
      est = ests,
      std = stds,
      group = group,
      fixed = computed,
      par = 0,
      stringsAsFactors = FALSE)
  }

  for (sfx in blockSuffix)
  {
    hasGGM <- paste0("omega", sfx) %in% allMatrices
    hasCor <- paste0("rho", sfx) %in% allMatrices
    hasChol <- paste0("lowertri", sfx) %in% allMatrices
    hasPrec <- paste0("kappa", sfx) %in% allMatrices
    nodes <- blockNodes(sfx)

    if (hasGGM || hasCor)
    {
      # Omega rows (op '--') stay as undirected edges and Rho rows (op '~~') as
      # (co)variance edges holding the correlations. The scaling rows (Delta or
      # SD):
      scaleRows <- pars$matrix == paste0(ifelse(hasGGM, "delta", "SD"), sfx)
      if (delta == "show")
      {
        # Show raw scaling parameters as self-loops:
        pars$op[scaleRows] <- "~~"
        pars$var2[scaleRows] <- pars$var1[scaleRows]
      } else
      {
        pars <- pars[!scaleRows, ]
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
        blockName <- switch(sfx, "_zeta" = "latent", "_epsilon" = "residual",
                            "_c" = "composite indicator", "observed")
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
    computedStd <- mapply(stdValue, pars$matrix, pars$var1_id, pars$var2_id, pars$est, pars$group,
                          SIMPLIFY = TRUE, USE.NAMES = FALSE)
    computedStd <- as.numeric(computedStd)
    computedStd[!is.finite(computedStd)] <- NA_real_
    pars$std <- ifelse(is.na(pars$std), computedStd, pars$std)
  }

  # Map psychonetrics operators to semPlot edge types. Factor loadings
  # (var2 =~ var1) and regressions (var1 <- var2) point from var2 to var1;
  # composite weights (var1 ~> var2) point from the indicator var1 into the
  # composite var2:
  known <- pars$op %in% c("~1","~=","<-","~>","~~","--")
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
           ifelse(pars$op %in% c("<-","~>"), "~>",
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

  # Composite weights:
  Pars$composite <- c(pars$op == "~>", rep(FALSE, nrow(Pars) - nrow(pars)))

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
    composite = c(rep(FALSE, length(manNames)), latNames %in% compNames),
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
