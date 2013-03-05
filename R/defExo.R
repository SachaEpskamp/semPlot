defExo <- function(object,layout="tree")
{
  manNames <- object@Vars$name[object@Vars$manifest]
  latNames <- object@Vars$name[!object@Vars$manifest]

  # Define exogenous variables (only if any is NA):
  if (any(is.na(object@Vars$exogenous)))
  {
    if (any(!is.na(object@Vars$exogenous)))
    {
      exoOrig <- object@Vars$exogenous
      repExo <- TRUE
    } else repExo <- FALSE
    object@Vars$exogenous <- FALSE
    for (i in which(!object@Vars$manifest))
    {
      if (!any(object@Pars$edge[object@Pars$rhs==object@Vars$name[i]] %in% c("~>","->") & object@Pars$lhs[object@Pars$rhs==object@Vars$name[i]]%in%latNames))
      {
        object@Vars$exogenous[i] <- TRUE
      }
    }
    for (i in which(object@Vars$manifest))
    {
      if (all(object@Pars$lhs[object@Pars$rhs==object@Vars$name[i] & object@Pars$lhs%in%latNames]%in%object@Vars$name[object@Vars$exogenous]) &
            all(object@Pars$rhs[object@Pars$lhs==object@Vars$name[i] & object@Pars$rhs%in%latNames]%in%object@Vars$name[object@Vars$exogenous]) &
            !any(object@Pars$rhs==object@Vars$name[i] & object@Pars$edge=="~>"))
      {
        object@Vars$exogenous[i] <- TRUE
      }
    }
    
    # If all exo, treat all as endo:
    if (all(object@Vars$exogenous) | layout%in%c("circle","circle2","circle3"))
    {
      object@Vars$exogenous <- FALSE
    }
    # If al endo, treat formative manifest as exo (MIMIC mode), unless all manifest are formative.
    if (!any(object@Vars$exogenous))
    {
      if (any(object@Vars$manifest & (object@Vars$name%in%object@Pars$rhs[object@Pars$edge %in% c("~>","--","->")])))
        object@Vars$exogenous[object@Vars$manifest & !(object@Vars$name%in%object@Pars$rhs[object@Pars$edge %in% c("~>","--","->")])] <- TRUE
    }
    if (repExo)
    {
      object@Vars$exogenous[!is.na(exoOrig)] <- exoOrig[!is.na(exoOrig)]
    }
  }
  
  return(object)
}