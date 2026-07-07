
semPlotModel.psem <- function(object, ...)
{
  if (!requireNamespace("piecewiseSEM", quietly = TRUE))
  {
    stop("The 'piecewiseSEM' package is required to import psem models. Install it with install.packages('piecewiseSEM').")
  }

  # coefs() returns a tidy data frame with (at least) the columns
  # Response, Predictor, Estimate, Std.Error, DF, Crit.Value, P.Value and
  # Std.Estimate. Warnings/messages (e.g. about standardization) are silenced.
  co <- suppressWarnings(suppressMessages(piecewiseSEM::coefs(object)))

  # Coerce defensively: older piecewiseSEM versions store these columns as
  # characters (with significance stars), in which case as.numeric() returns NA,
  # which is acceptable.
  est <- suppressWarnings(as.numeric(co$Estimate))
  std <- suppressWarnings(as.numeric(co$Std.Estimate))

  # Correlated errors (declared as x %~~% y) appear as rows whose Response and
  # Predictor both carry a "~~" prefix, with the correlation stored in Estimate.
  isCorr <- grepl("^~~", co$Response)

  lhs  <- ifelse(isCorr, sub("^~~", "", co$Predictor), co$Predictor)
  rhs  <- ifelse(isCorr, sub("^~~", "", co$Response),  co$Response)
  edge <- ifelse(isCorr, "<->", "~>")

  # For correlated errors Estimate is itself a correlation, so it doubles as the
  # standardized estimate.
  stdOut <- ifelse(isCorr, est, std)

  # Factor-level dummies (Predictor like "groupB") are imported as-is: the level
  # name is treated as an ordinary predictor variable.
  Pars <- data.frame(
    label = "",
    lhs = lhs,
    edge = edge,
    rhs = rhs,
    est = est,
    std = stdOut,
    group = "",
    fixed = FALSE,
    par = seq_len(nrow(co)),
    knot = 0,
    stringsAsFactors = FALSE)

  ## Split interactions:
  # coefs() reports products as "a:b" in Predictor (now lhs). Mirror the lm
  # importer: keep the first component in lhs with a fresh knot id and append
  # duplicate rows for the remaining components.
  if (any(grepl(":", Pars$lhs)))
  {
    colons <- grep(":", Pars$lhs)
    for (i in seq_along(colons))
    {
      labs <- strsplit(Pars$lhs[colons[i]], split = ":")[[1]]
      Pars$lhs[colons[i]] <- labs[1]
      Pars$knot[colons[i]] <- i
      for (j in 2:length(labs))
      {
        Pars <- rbind(Pars, Pars[colons[i], ])
        Pars$lhs[nrow(Pars)] <- labs[j]
      }
    }
  }

  # Variable dataframe:
  Vars <- data.frame(
    name = unique(c(Pars$lhs, Pars$rhs)),
    manifest = TRUE,
    exogenous = NA,
    stringsAsFactors = FALSE)
  Vars <- Vars[Vars$name != "", ]

  semModel <- new("semPlotModel")
  semModel@Pars <- Pars
  semModel@Vars <- Vars
  semModel@Computed <- TRUE
  semModel@Original <- list(object)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  semModel@Thresholds <- data.frame()

  return(semModel)
}
