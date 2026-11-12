# Add function:
'+.semPlotModel' <- function(x,y)
{
  stopifnot("semPlotModel"%in%class(x))
  stopifnot("semPlotModel"%in%class(y))
  
  # Update par in y:
  y@Pars$par[y@Pars$par>0] <- max(x@Pars$par) + y@Pars$par[y@Pars$par>0]
  
  # Composite columns, so that both models have them (and the same node is
  # not listed twice because one model lacks the column):
  x <- .normalizeComposites(x)
  y <- .normalizeComposites(y)

  # New model. Pars/Thresholds columns can differ between importers (e.g.
  # 'knot' from the lm importer, 'composite'/'BetweenWithin' from lavaan),
  # so bind with column union rather than plain rbind:
  semModel <- new("semPlotModel")
  semModel@Pars <- plyr::rbind.fill(x@Pars,y@Pars)
  semModel@Vars <- plyr::rbind.fill(x@Vars,y@Vars)
  semModel@Vars <- semModel@Vars[!duplicated(semModel@Vars),]
  semModel@Thresholds <- if (nrow(x@Thresholds) > 0 || nrow(y@Thresholds) > 0) plyr::rbind.fill(x@Thresholds,y@Thresholds) else data.frame()
  semModel@Computed <- x@Computed && y@Computed
  semModel@Original <- list(x@Original[[1]],y@Original[[1]])
  semModel@ObsCovs <- c(x@ObsCovs,y@ObsCovs)
  semModel@ImpCovs <- c(x@ImpCovs,y@ImpCovs)
  
  # Return:
  return(semModel)
}