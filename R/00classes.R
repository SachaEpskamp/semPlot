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

# Registry for third-party importers: packages (or users) can register a
# converter for their model class without modifying semPlot. See
# ?registerSemPlotImporter.
.importerRegistry <- new.env(parent = emptyenv())

registerSemPlotImporter <- function(class, fun)
{
  stopifnot(is.character(class), length(class) == 1, nzchar(class), is.function(fun))
  assign(class, fun, envir = .importerRegistry)
  invisible(NULL)
}

# Validate the invariants semPaths() relies on. Returns the object
# invisibly, or stops with a message listing all violations.
validateSemPlotModel <- function(object)
{
  problems <- character(0)
  if (!is(object, "semPlotModel"))
  {
    stop("Not a semPlotModel object.")
  }
  needed <- c("label","lhs","edge","rhs","est","std","group","fixed","par")
  missingCols <- needed[!needed %in% names(object@Pars)]
  if (length(missingCols) > 0)
  {
    problems <- c(problems, paste0("Pars is missing column(s): ", paste(missingCols, collapse = ", ")))
  }
  if ("edge" %in% names(object@Pars))
  {
    known <- c("->","~>","<->","--","int","|")
    bad <- unique(object@Pars$edge[!object@Pars$edge %in% known])
    if (length(bad) > 0)
    {
      problems <- c(problems, paste0("Unknown edge type(s): ", paste(bad, collapse = ", "),
                                     " (known: ", paste(known, collapse = ", "), ")"))
    }
    if (any(is.na(object@Pars$edge)))
    {
      problems <- c(problems, "Pars$edge contains NA")
    }
  }
  if (all(c("lhs","rhs","edge") %in% names(object@Pars)) && "name" %in% names(object@Vars))
  {
    nodes <- unique(c(object@Pars$lhs[object@Pars$edge != "int"], object@Pars$rhs))
    unknown <- nodes[!nodes %in% c(object@Vars$name, "")]
    if (length(unknown) > 0)
    {
      problems <- c(problems, paste0("Pars refer to variable(s) not in Vars: ", paste(unknown, collapse = ", ")))
    }
  }
  if (length(problems) > 0)
  {
    stop("Invalid semPlotModel:\n", paste0("  - ", problems, collapse = "\n"))
  }
  invisible(object)
}

semPlotModel <- function (object, ...) {
  # Check if the *unevaluated* argument is a `+` call, if so combine models:
  combined <- .tryCombineModels(substitute(object), parent.frame())
  if (!is.null(combined)) return(combined)

  # semPlotModel objects pass through unchanged (checked before the generic
  # S4 branch below, which would otherwise fail to dispatch on them):
  if (is(object, "semPlotModel")) return(object)

  # Registered third-party importers take precedence over built-ins:
  for (cl in class(object))
  {
    if (exists(cl, envir = .importerRegistry, inherits = FALSE))
    {
      return(get(cl, envir = .importerRegistry)(object, ...))
    }
  }

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
      if (!"try-error"%in%class(mod)) return(mod) else stop("Input string is neither an existing file nor a lavaan model. lavaan importer said:\n  ", attr(mod, "condition")$message, call. = FALSE)
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
    head <- readLines(object, n = 100, warn = FALSE)
    if (any(grepl("mplus",head,ignore.case=TRUE)))
    {
      return(semPlotModel.mplus.model(object,...))
    }
    
    if (any(grepl("l\\s*i\\s*s\\s*r\\s*e\\s*l",head,ignore.case=TRUE)))
    {
      if (!requireNamespace("lisrelToR", quietly = TRUE))
      {
        stop("This looks like LISREL output; the 'lisrelToR' package is required to import it. Install it with install.packages('lisrelToR').")
      }
      return(semPlotModel(lisrelToR::readLisrel(object)))
    }
    
    # If all else fais, just try everything and assume you get errors
    # if it is wrong. Accumulate each importer's error message so a
    # valid-but-failing input gets a useful diagnostic:
    attempts <- list(
      lavaan = function() semPlotModel_lavaanModel(object, ...),
      Mplus  = function() semPlotModel.mplus.model(object, ...),
      LISREL = function() {
        if (!requireNamespace("lisrelToR", quietly = TRUE)) stop("package 'lisrelToR' is not installed")
        semPlotModel(lisrelToR::readLisrel(object))
      },
      Onyx   = function() semPlotModel_Onyx(object),
      Amos   = function() semPlotModel_Amos(object)
    )
    errs <- character(0)
    for (nm in names(attempts))
    {
      mod <- try(attempts[[nm]](), silent = TRUE)
      if (!"try-error" %in% class(mod)) return(mod)
      errs[nm] <- attr(mod, "condition")$message
    }
    stop("Object not recognized as a SEM model. Importer attempts failed with:\n",
         paste0("  - ", names(errs), ": ", errs, collapse = "\n"), call. = FALSE)
  }

  stop("Object not recognized as SEM model")
}
