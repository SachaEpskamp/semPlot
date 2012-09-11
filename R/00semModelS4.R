
setClass( "SEMmodel", representation(
    RAM = "data.frame",
    Vars = "data.frame",
    Computed = "logical",
    ObsCovs = "list",
    ImpCovs = "list",
    Original = "list"))

setGeneric("SEMmodel.S4", function(object) {
  standardGeneric("SEMmodel.S4")
})

setGeneric("SEMpaths.S4", function(object,...) {
  standardGeneric("SEMpaths.S4")
})

SEMpaths <- function(object,...)
{
  if ("MxRAMModel"%in%class(object)) return(SEMpaths_MxRAMModel(object,...)) 
  if ("MxModel"%in%class(object)) return(SEMpaths_MxModel(object,...))
  if(isS4(object)) 
  {
    SEMpaths.S4(object, ...)
  } else
  {
    UseMethod("SEMpaths", object)
  }
}

SEMmodel <- function (object) {
  if ("MxRAMModel"%in%class(object)) return(SEMmodel_MxRAMModel(object)) 
  if ("MxModel"%in%class(object)) return(SEMmodel_MxModel(object))
  if(isS4(object)) 
  {
    SEMmodel.S4(object)
  } else
  {
    UseMethod("SEMmodel", object)
  }
}

