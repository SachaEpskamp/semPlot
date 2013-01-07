# semPaths.principal <- function(object,...) 
# {
#   invisible(semPaths(semPlotModel(object),...))
# }
# 


### SINGLE GROUP MODEL ###
semPlotModel.principal <- function(object)
{
  
  # Check if object is of class "sem":
  if (!"principal"%in%class(object)) stop("Input must be a 'principal' object")
  
  # Extract model:
  mod <- semPlotModel(loadings(object))
  manNames <- mod@Vars$name[mod@Vars$manifest]
  
  # Fix:
  mod@RAM[c("lhs","rhs")] <- mod@RAM[c("rhs","lhs")]
  mod@RAM$edge <- "->"
    
  return(mod)
}

