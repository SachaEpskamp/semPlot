# semPaths.factanal <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 


### SINGLE GROUP MODEL ###
semPlotModel.factanal <- function(object)
{
  
  # Check if object is of class "sem":
  if (!"factanal"%in%class(object)) stop("Input must be a 'factanal' object")
  
  
  # Extract model:
  mod <- semPlotModel(loadings(object))
  manNames <- mod@Vars$name[mod@Vars$manifest]
  
  # Fix:
  mod@RAM$edge <- "->"
  
  # Add residuals:
  Uniqueness <- object$uniquenesses
  
  residRAM  <- data.frame(
    label = "", 
    lhs = manNames,
    edge = "<->",
    rhs = manNames,
    est = Uniqueness,
    std = Uniqueness,
    group = "",
    fixed = FALSE,
    par = 0,
    stringsAsFactors=FALSE)
  
  mod@RAM <- rbind(mod@RAM,residRAM)
  mod@RAM$par <- 1:nrow(mod@RAM)
  
  
  return(mod)
}

