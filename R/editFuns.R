# Extract exogenous variables:
exo <- function(x)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$name[!is.na(x@Vars$exogenous)][x@Vars$exogenous[!is.na(x@Vars$exogenous)]]
}

# Set exogenous variables:
"exo<-" <- function(x,value)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$name[!is.na(x@Vars$exogenous)][x@Vars$exogenous[!is.na(x@Vars$exogenous)]] <- FALSE
  x@Vars$exogenous[x@Vars$name%in%value] <- TRUE
  return(x)
}

# Extract endogenous variables:
endo <- function(x)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$name[!is.na(x@Vars$exogenous)][!x@Vars$exogenous[!is.na(x@Vars$exogenous)]]
}

# Set endogenous variables:
"endo<-" <- function(x,value)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$name[!is.na(x@Vars$exogenous)][!x@Vars$exogenous[!is.na(x@Vars$exogenous)]] <- TRUE
  x@Vars$exogenous[x@Vars$name%in%value] <- FALSE
  return(x)
}

# Extract manifest variables:
man <- function(x)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$name[x@Vars$manifest]
}

# Set manifest variables:
"man<-" <- function(x,value)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$manifest[x@Vars$name%in%value] <- TRUE
  return(x)
}

# Extract latent variables:
lat <- function(x)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$name[!x@Vars$manifest]
}

# Set latent variables:
"lat<-" <- function(x,value)
{
  if (!"semPlotModel"%in%class(x)) stop("'semPlotModel' object is required")
  x@Vars$manifest[x@Vars$name%in%value] <- FALSE
  return(x)
}