## SemPlotModel
# Note on edge specification:
# '->' is factor loading
# '~>' is regression
# '<->' is (co)variance
# 'int' is an intercept

setClass( "semPlotModel", representation(
  Pars = "data.frame",
  Vars = "data.frame",
  Thresholds = "data.frame",
  Computed = "logical",
  ObsCovs = "list",
  ImpCovs = "list",
  Original = "list"))

setGeneric("semPlotModel.S4", function(object,...) {
  standardGeneric("semPlotModel.S4")
})
# 
# setGeneric("semPaths.S4", function(object,...) {
#   standardGeneric("semPaths.S4")
# })
# 
# semPaths <- function(object,...)
# {
#   if ("MxRAMModel"%in%class(object)) return(semPaths_MxRAMModel(object,...)) 
#   if ("MxModel"%in%class(object)) return(semPaths_MxModel(object,...))
#   if(isS4(object)) 
#   {
#     semPaths.S4(object, ...)
#   } else
#   {
#     UseMethod("semPaths", object)
#   }
# }

semPlotModel <- function (object, ...) {
  # Check if call contains a + operator, if so combine models:
  
  call <- paste(deparse(substitute(object)), collapse = "")
  if (grepl("\\+",call) & !grepl("\"",call) & !grepl("\'",call)) 
  {
    args <- unlist(strsplit(call,split="\\+"))
    obs <- lapply(args,function(x)semPlotModel(eval(parse(text=x))))
    Res <- obs[[1]]
    for (i in 2:length(obs)) Res <- Res + obs[[i]]
    return(Res)
  }
  
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

semPlotModel.semPlotModel <- function(object,...) object


# semPaths.default <- function(object,...)
# {
#   if (is.character(object) && grepl("\\.out",object))
#   {
#     return(semPaths(readModels(object),...))
#   }
# }

semPlotModel.default <- function(object,...)
{
  if (is(object,'data.frame'))
  {
    mod <- try(semPlotModel_lavaanModel(object,...),silent=TRUE)
    if (!"try-error"%in%class(mod)) return(mod)
  }
  
  if (is.character(object))
  {
    if (!file.exists(object))
    {
      mod <- try(semPlotModel_lavaanModel(object,...),silent=TRUE)
      if (!"try-error"%in%class(mod)) return(mod) else stop("Input string neither an existing file or Lavaan model.")
    }
    # Find file:
    if (grepl("\\.xml",object,ignore.case=TRUE))
    {
      return(semPlotModel_Onyx(object))
    }
    if (grepl("\\.AmosOutput",object,ignore.case=TRUE))
    {
      return(semPlotModel_Amos(object))
    }
    
    # Read first 100 lines:
    head <- readLines(object, 10)
    if (any(grepl("mplus",head,ignore.case=TRUE)))
    {
      return(semPlotModel.mplus.model(object))
    }
    
    if (any(grepl("l\\s*i\\s*s\\s*r\\s*e\\s*l",head,ignore.case=TRUE)))
    {
      return(semPlotModel(readLisrel(object)))
    }
    
    # If all else fais, just try everything and assume you get errors 
    # if it is wrong:
    mod <- try(semPlotModel_lavaanModel(object,...),silent=TRUE)
    if (!"try-error"%in%class(mod)) return(mod)
    
    mod <- try(semPlotModel.mplus.model(object),silent=TRUE)
    if (!"try-error"%in%class(mod)) return(mod)

    mod <- try(semPlotModel(readLisrel(object)),silent=TRUE)
    if (!"try-error"%in%class(mod)) return(mod)
    
    mod <- try(semPlotModel_Onyx(object),silent=TRUE)
    if (!"try-error"%in%class(mod)) return(mod)
    
    mod <- try(semPlotModel_Amos(object),silent=TRUE)
    if (!"try-error"%in%class(mod)) return(mod)
    
    # Well, we failed...
  }
  
  stop("Object not recognized as SEM model")
}
