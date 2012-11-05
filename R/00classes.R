
setClass( "semPlotModel", representation(
    RAM = "data.frame",
    Vars = "data.frame",
    Computed = "logical",
    ObsCovs = "list",
    ImpCovs = "list",
    Original = "list"))

setGeneric("semPlotModel.S4", function(object) {
  standardGeneric("semPlotModel.S4")
})

setGeneric("semPaths.S4", function(object,...) {
  standardGeneric("semPaths.S4")
})

semPaths <- function(object,...)
{
  if ("MxRAMModel"%in%class(object)) return(semPaths_MxRAMModel(object,...)) 
  if ("MxModel"%in%class(object)) return(semPaths_MxModel(object,...))
  if(isS4(object)) 
  {
    semPaths.S4(object, ...)
  } else
  {
    UseMethod("semPaths", object)
  }
}

semPlotModel <- function (object) {
  if ("MxRAMModel"%in%class(object)) return(semPlotModel_MxRAMModel(object)) 
  if ("MxModel"%in%class(object)) return(semPlotModel_MxModel(object))
  if(isS4(object)) 
  {
    semPlotModel.S4(object)
  } else
  {
    UseMethod("semPlotModel", object)
  }
}


semPaths.default <- function(object,...)
{
  if (is.character(object) && grepl("\\.out",object))
  {
    return(semPaths(readModels(object),...))
  }
}

semPlotModel.default <- function(object)
{
  if (is.character(object) && grepl("\\.out",object))
  {
    return(semPlotModel(readModels(object)))
  }
}
