# Importers for metaSEM (meta-analytic SEM). Stage-2 objects ('wls', from
# tssem2()/wls()) are converted through metaSEM's own bridge meta2semPlot(),
# which returns a semPlotModel. Stage-1 pooling objects carry no path
# structure and produce an informative error.

semPlotModel.wls <- function(object, ...)
{
  if (!requireNamespace("metaSEM", quietly = TRUE))
  {
    stop("The 'metaSEM' package is required to import metaSEM models. Install it with install.packages('metaSEM').")
  }
  # Drop semPaths' internal modelOpts arguments (e.g. mplusStd), which
  # meta2semPlot()'s downstream functions do not accept:
  dots <- list(...)
  dots$mplusStd <- NULL
  do.call(metaSEM::meta2semPlot, c(list(object), dots))
}

semPlotModel.meta <- function(object, ...)
{
  stop("This is a stage-1 metaSEM object (pooled correlations); it contains no path structure to draw. Fit the structural model with tssem2() or wls() and plot that instead.")
}

semPlotModel.tssem1FEM <- semPlotModel.meta
