# Add function:
'+.semPlotModel' <- function(x,y)
{
  stopifnot("semPlotModel"%in%class(x))
  stopifnot("semPlotModel"%in%class(y))
  
  # Update par in y:
  y@RAM$par[y@RAM$par>0] <- max(x@RAM$par) + y@RAM$par[y@RAM$par>0]
  
  # New model:
  semModel <- new("semPlotModel")
  semModel@RAM <- rbind(x@RAM,y@RAM)
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