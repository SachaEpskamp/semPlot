semSyntax <- function(object, syntax = "lavaan", allFixed = FALSE, file)
{
  if (!"semPlotModel" %in% class(object))
  {
    # Try to run semPlotModel on object, otherwise stop.
    object <- semPlotModel(object)
  }
  if (!syntax %in% c("lavaan","sem")) stop("Only 'lavaan' and 'sem' syntax is currently supported ")
 
  if (nrow(object@Thresholds) > 0) warning("Thresholds are not yet supported by semSyntax")
  
  # If all fixed, simply set all fixed = TRUE:
  if (allFixed)
  {
    object@Pars$fixed <- TRUE
  }
  
  ### LAVAAN ###
  if (syntax == "lavaan")
  {   
    Pars <- object@Pars
    
    # Reverse lhs and rhs:
    Pars[Pars$edge %in% c('~>','int'),c('lhs','rhs')] <- Pars[Pars$edge %in% c('~>','int'),c('rhs','lhs')]
    Pars[Pars$edge=='->'&!(Pars$lhs%in%object@Vars$name[!object@Vars$manifest] & Pars$rhs%in%object@Vars$name[object@Vars$manifest]),c('lhs','rhs')] <- Pars[Pars$edge=='->'&!(Pars$lhs%in%object@Vars$name[!object@Vars$manifest] & Pars$rhs%in%object@Vars$name[object@Vars$manifest]),c('rhs','lhs')]
    
    # Change operators:
    Pars$edge[Pars$edge=='->'&!(Pars$lhs%in%object@Vars$name[!object@Vars$manifest] & Pars$rhs%in%object@Vars$name[object@Vars$manifest])] <- "~"
    Pars$edge[Pars$edge=='->'&(Pars$lhs%in%object@Vars$name[!object@Vars$manifest] & Pars$rhs%in%object@Vars$name[object@Vars$manifest])] <- "=~"
    Pars$edge[Pars$edge == "~>"] <- "~"
    Pars$edge[Pars$edge == "<->"] <- "~~"
    Pars$rhs[Pars$edge == "int"] <- "1"
    Pars$edge[Pars$edge == "int"] <- "~"
    
    
    # Fixing parameters:
    Pars$rhs <- ifelse( Pars$fixed, paste0(Pars$est,"*",Pars$rhs), Pars$rhs)
    Pars$rhs <- ifelse( !Pars$fixed & Pars$par > 0 & (duplicated(Pars$par)|duplicated(Pars$par,fromLast=TRUE)), paste0("par",Pars$par,"*",Pars$rhs), Pars$rhs)
    
    # Combine and return:
    Mod <- paste(Pars$lhs,Pars$edge,Pars$rhs,collapse = "\n")
    
    # Print to console or file:
    if (missing(file))
    {
      cat("\nModel <- '\n",Mod,"\n'\n",sep="")
    } else 
    {
      write(paste0("\nModel <- '\n",Mod,"\n'\n"),file)
    }
    
    return(Mod)
  }
  
  ### SEM ###
  if (syntax == "sem")
  {
    browser()
    
    Pars <- object@Pars
    
    # Remove intercepts:
    if (any(Pars$edge == "int"))
    {
      warning("Intercepts removed from model for 'sem' syntax")
      Pars <- Pars[Pars$edge!="int",]
    }
    
    Pars$label[Pars$fixed] <- NA
    ## Fix parameter labels.
    if (max(Pars$par) > 0)
    {
      for (i in seq_len(max(Pars$par)))
      {
        # Check if unique to other par numbers:
        if (any(Pars$label[Pars$par!=i] %in% Pars$label[Pars$par==i] | any(Pars$label[Pars$par == i] == '')))
        {
          Pars$label[Pars$par==i] <- paste0("par",i)
        }
        
        # Check if labels are unique, else combine:
        if (length(unique(Pars$label[Pars$par == i])) > 1)
        {
          Pars$label[Pars$par==i] <- paste(Pars$label[Pars$par==i],collapse="_")
        }
      }
    }
    
    # Fix estimate:
    Pars$est[!Pars$fixed] <- NA
    
    # Fix edges:
    Pars$edge[Pars$edge == '~>'] <- '->'
    
    # Create model:
    Mod <- paste(paste(Pars$lhs, Pars$edge, Pars$rhs), Pars$label, Pars$est, sep = ",", collapse = "\n")
    
    # Print to console or file:
    if (missing(file))
    {
      cat("\nModel <- specifyModel()\n",Mod,"\n\n",sep="")
    } else 
    {
      write(paste0("\nModel <- specifyModel()\n",Mod,"\n\n",sep=""),file)
    }

    Mod <- specifyModel( textConnection( Mod ))
    
    return(Mod)
  }
  
  
  
}