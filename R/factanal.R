# semPaths.factanal <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 


### SINGLE GROUP MODEL ###
semPlotModel.factanal <- function(object, ...)
{
  
  # Check if object is of class "sem":
  if (!"factanal"%in%class(object)) stop("Input must be a 'factanal' object")
  
  
  # Extract model:
  mod <- semPlotModel(loadings(object))
  manNames <- mod@Vars$name[mod@Vars$manifest]
  
  # Fix:
  mod@Pars$edge <- "->"
  
  # Add residuals:
  Uniqueness <- object$uniquenesses
  
  residPars  <- data.frame(
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
  
  mod@Pars <- rbind(mod@Pars,residPars)
  mod@Pars$par <- 1:nrow(mod@Pars)
  
  
  return(mod)
}

