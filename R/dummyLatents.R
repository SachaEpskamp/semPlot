# Dummy latent variables (phantom single-indicator latents). In LISREL-style
# frameworks, such as psychonetrics' lvm() and cvm(), an observed variable can
# only enter the structural model through a latent variable that it measures
# perfectly: eta =~ 1*y with a residual variance of y fixed to zero, so y is
# eta. Such a dummy latent adds nothing to the path diagram, so semPaths()
# (with hideDummies = TRUE) draws the model with y in its place.
#
# A latent variable eta is a dummy latent if:
#   - it is not a composite;
#   - its only loading ('->' edge) goes to a single observed variable y and is
#     fixed to 1 (in every group);
#   - the residual variance of y is fixed to zero: a fixed '<->' row with
#     estimate 0, or no row at all (most importers, including psychonetrics,
#     omit parameters fixed to zero);
#   - y has no other incoming directed edges (loadings or regressions) and no
#     non-zero (or free) covariances with other variables;
#   - eta and y do not both have an intercept that is not fixed to zero.
#
# Only estimated models (and psychonetrics models, which omit all parameters
# fixed to zero) are considered: in a model given only as syntax, unspecified
# variances are listed as fixed to zero, and a latent drawn by the user stays.
#
# Returns the names of the dummy latents, named by their indicator.
.dummyLatents <- function(object)
{
  Pars <- object@Pars
  Vars <- object@Vars
  if (nrow(Pars) == 0 || nrow(Vars) == 0) return(character(0))
  if (!(isTRUE(object@Computed) ||
        any(vapply(object@Original, function(x) is(x, "psychonetrics"), logical(1))))) return(character(0))

  isZero <- function(rows) all(Pars$fixed[rows] %in% TRUE) && all(Pars$est[rows] %in% 0)
  isNone <- function(rows) length(rows) == 0 || isZero(rows)

  latNames <- Vars$name[!Vars$manifest & !Vars$composite]
  manNames <- Vars$name[Vars$manifest]

  dummies <- character(0)
  for (eta in latNames)
  {
    load <- which(Pars$lhs == eta & Pars$edge == "->")
    if (length(load) == 0) next
    y <- unique(Pars$rhs[load])
    if (length(y) != 1 || !y %in% manNames) next
    if (!all(Pars$fixed[load] %in% TRUE) || !all(Pars$est[load] %in% 1)) next
    # Do not chain dummies onto the same indicator:
    if (y %in% names(dummies)) next

    # Residual variance of y fixed to zero:
    if (!isNone(which(Pars$lhs == y & Pars$rhs == y & Pars$edge == "<->"))) next

    # No other incoming directed edges into y:
    inc <- which(Pars$rhs == y & Pars$edge %in% c("->", "~>") & Pars$lhs != eta)
    if (length(inc) > 0) next

    # No covariances / undirected edges between y and other variables:
    cov <- which(Pars$edge %in% c("<->", "--") & Pars$lhs != Pars$rhs &
                   (Pars$lhs == y | Pars$rhs == y))
    if (!isNone(cov)) next

    # At most one non-zero intercept:
    if (!isNone(which(Pars$edge == "int" & Pars$rhs == eta)) &&
        !isNone(which(Pars$edge == "int" & Pars$rhs == y))) next

    dummies[y] <- eta
  }

  dummies
}

# Collapse the dummy latents found by .dummyLatents() onto their indicators:
.collapseDummyLatents <- function(object, dummies)
{
  if (length(dummies) == 0) return(object)
  Pars <- object@Pars

  for (y in names(dummies))
  {
    eta <- dummies[[y]]

    # Remove the loading, the zero residual variance and zero covariances of y,
    # and the intercept that is fixed to zero:
    drop <- (Pars$lhs == eta & Pars$rhs == y & Pars$edge == "->") |
      ((Pars$lhs == y | Pars$rhs == y) & Pars$edge %in% c("<->", "--"))
    for (v in c(eta, y))
    {
      int <- Pars$edge == "int" & Pars$rhs == v
      if (any(int) && all(Pars$fixed[int] %in% TRUE) && all(Pars$est[int] %in% 0)) drop <- drop | int
    }
    Pars <- Pars[!drop, , drop = FALSE]

    # y takes the place of eta:
    Pars$lhs[Pars$lhs == eta] <- y
    Pars$rhs[Pars$rhs == eta] <- y
  }
  rownames(Pars) <- NULL
  object@Pars <- Pars
  object@Vars <- object@Vars[!object@Vars$name %in% dummies, , drop = FALSE]
  rownames(object@Vars) <- NULL

  object
}
