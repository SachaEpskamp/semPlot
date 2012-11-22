semPaths.princomp <- function(object,...) 
{
  invisible(semPaths(semPlotModel(object),...))
}



### SINGLE GROUP MODEL ###
semPlotModel.princomp <- function(object)
{
  
  # Check if object is of class "sem":
  if (!"princomp"%in%class(object)) stop("Input must be a 'princomp' object")
  
  
  # Extract model:
  mod <- semPlotModel(loadings(object))
  manNames <- mod@Vars$name[mod@Vars$manifest]
  
  # Fix:
  mod@RAM[c("lhs","rhs")] <- mod@RAM[c("rhs","lhs")]
  mod@RAM$edge <- "->"
    
  return(mod)
}

