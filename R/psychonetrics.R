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

  # Variable names (order matches psychonetrics' matrix rows):
  manNames <- object@sample@variables$label
  latNames <- unique(pars$var2[pars$matrix == "lambda"])
  latNames <- latNames[!is.na(latNames)]

  # Thresholds not (yet) supported:
  if (any(pars$matrix == "tau"))
  {
    message("Thresholds ('tau') of psychonetrics models are not yet supported by semPlot and are not shown.")
    pars <- pars[pars$matrix != "tau", ]
  }

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

  # Rows that only encode a parameterization (not path-diagram parameters):
  parameterizationRows <- function(sfx)
  {
    pars$matrix %in% paste0(c("delta","lowertri","kappa"), sfx)
  }

  derivedPars <- list()
  addDerived <- function(nodes1, nodes2, ests, group)
  {
    derivedPars[[length(derivedPars) + 1]] <<- data.frame(
      label = "",
      lhs = nodes1,
      edge = "<->",
      rhs = nodes2,
      est = ests,
      std = NA,
      group = group,
      fixed = FALSE,
      par = 0,
      stringsAsFactors = FALSE)
  }

  for (sfx in blockSuffix)
  {
    matNames <- unique(pars$matrix)
    hasGGM <- paste0("omega", sfx) %in% matNames
    hasChol <- paste0("lowertri", sfx) %in% matNames
    hasPrec <- paste0("kappa", sfx) %in% matNames
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
              addDerived(nodes, nodes, diag(sig[[g]]), g)
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
          S <- sig[[g]]
          ind <- which(upper.tri(S, diag = TRUE), arr.ind = TRUE)
          ests <- S[ind]
          nonZero <- abs(ests) > sqrt(.Machine$double.eps)
          addDerived(nodes[ind[nonZero, 1]], nodes[ind[nonZero, 2]], ests[nonZero], g)
        }
      }
    }
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
