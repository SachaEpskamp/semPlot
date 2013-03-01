# Add function:
'+.semPlotModel' <- function(x,y)
{
  stopifnot("semPlotModel"%in%class(x))
  stopifnot("semPlotModel"%in%class(y))
  
  # Update par in y:
  y@Pars$par[y@Pars$par>0] <- max(x@Pars$par) + y@Pars$par[y@Pars$par>0]
  
  # New model:
  semModel <- new("semPlotModel")
  semModel@Pars <- rbind(x@Pars,y@Pars)
  semModel@Vars <- rbind(x@Vars,y@Vars)
  semModel@Vars <- semModel@Vars[!duplicated(semModel@Vars),]
  semModel@Thresholds <- rbind(x@Thresholds,y@Thresholds)
  semModel@Computed <- x@Computed && y@Computed
  semModel@Original <- list(x@Original[[1]],y@Original[[1]])
  semModel@ObsCovs <- c(x@ObsCovs,y@ObsCovs)
  semModel@ImpCovs <- c(x@ImpCovs,y@ImpCovs)
  
  # Return:
  return(semModel)
}