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

setGeneric("semPlotModel_S4", function(object,...) {
  standardGeneric("semPlotModel_S4")
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

# Collect the terms of a (possibly nested) `a + b (+ c ...)` call,
# evaluated in 'env'. Only used for expressions whose top-level call is
# the binary `+` operator.
.collect_plus_terms <- function(expr, env)
{
  if (is.call(expr) && length(expr) == 3 && identical(expr[[1]], as.name("+")))
  {
    c(.collect_plus_terms(expr[[2]], env), .collect_plus_terms(expr[[3]], env))
  } else {
    list(eval(expr, env))
  }
}

# Interpret 'expr' (the substitute()d argument) as a combination of models
# written as `a + b`, for model classes that do not define `+` themselves
# (e.g. two lavaan fits). Returns a combined semPlotModel, or NULL when the
# expression is not a `+` call, a term does not evaluate, or R can evaluate
# the sum natively (in which case the caller proceeds normally). This
# replaces an old deparse-and-split-on-"+" heuristic that crashed on any
# inline call containing a `+`, such as semPaths(lm(y ~ x + z, data)).
.tryCombineModels <- function(expr, env)
{
  if (!(is.call(expr) && length(expr) == 3 && identical(expr[[1]], as.name("+")))) return(NULL)
  terms <- tryCatch(.collect_plus_terms(expr, env), error = function(e) NULL)
  if (is.null(terms)) return(NULL)
  # If R can evaluate the sum itself (semPlotModel objects, numerics, ...),
  # let the normal path handle it:
  canAdd <- tryCatch({ Reduce("+", terms); TRUE }, error = function(e) FALSE)
  if (canAdd) return(NULL)
  # Otherwise interpret `+` as model combination (errors propagate):
  Reduce("+", lapply(terms, semPlotModel))
}

semPlotModel <- function (object, ...) {
  # Check if the *unevaluated* argument is a `+` call, if so combine models:
  combined <- .tryCombineModels(substitute(object), parent.frame())
  if (!is.null(combined)) return(combined)

  # semPlotModel objects pass through unchanged (checked before the generic
  # S4 branch below, which would otherwise fail to dispatch on them):
  if (is(object, "semPlotModel")) return(object)

  if (is(object, "psychonetrics")) return(semPlotModel_psychonetrics(object, ...))
  if ("MxRAMModel"%in%class(object)) return(semPlotModel_MxRAMModel(object))
  if ("MxModel"%in%class(object)) return(semPlotModel_MxModel(object))
  if(isS4(object)) 
  {
    semPlotModel_S4(object)
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
      return(semPlotModel.mplus.model(object,...))
    }
    
    if (any(grepl("l\\s*i\\s*s\\s*r\\s*e\\s*l",head,ignore.case=TRUE)))
    {
      return(semPlotModel(readLisrel(object)))
    }
    
    # If all else fais, just try everything and assume you get errors 
    # if it is wrong:
    mod <- try(semPlotModel_lavaanModel(object,...),silent=TRUE)
    if (!"try-error"%in%class(mod)) return(mod)
    
    mod <- try(semPlotModel.mplus.model(object,...),silent=TRUE)
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
