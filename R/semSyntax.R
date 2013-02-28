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
    object@RAM$fixed <- TRUE
  }
  
  ### LAVAAN ###
  if (syntax == "lavaan")
  {   
    RAM <- object@RAM
    
    # Reverse lhs and rhs:
    RAM[RAM$edge %in% c('~>','int'),c('lhs','rhs')] <- RAM[RAM$edge %in% c('~>','int'),c('rhs','lhs')]
    RAM[RAM$edge=='->'&!(RAM$lhs%in%object@Vars$name[!object@Vars$manifest] & RAM$rhs%in%object@Vars$name[object@Vars$manifest]),c('lhs','rhs')] <- RAM[RAM$edge=='->'&!(RAM$lhs%in%object@Vars$name[!object@Vars$manifest] & RAM$rhs%in%object@Vars$name[object@Vars$manifest]),c('rhs','lhs')]
    
    # Change operators:
    RAM$edge[RAM$edge=='->'&!(RAM$lhs%in%object@Vars$name[!object@Vars$manifest] & RAM$rhs%in%object@Vars$name[object@Vars$manifest])] <- "~"
    RAM$edge[RAM$edge=='->'&(RAM$lhs%in%object@Vars$name[!object@Vars$manifest] & RAM$rhs%in%object@Vars$name[object@Vars$manifest])] <- "=~"
    RAM$edge[RAM$edge == "~>"] <- "~"
    RAM$edge[RAM$edge == "<->"] <- "~~"
    RAM$rhs[RAM$edge == "int"] <- "1"
    RAM$edge[RAM$edge == "int"] <- "~"
    
    
    # Fixing parameters:
    RAM$rhs <- ifelse( RAM$fixed, paste0(RAM$est,"*",RAM$rhs), RAM$rhs)
    RAM$rhs <- ifelse( !RAM$fixed & RAM$par > 0 & (duplicated(RAM$par)|duplicated(RAM$par,fromLast=TRUE)), paste0("par",RAM$par,"*",RAM$rhs), RAM$rhs)
    
    # Combine and return:
    Mod <- paste(RAM$lhs,RAM$edge,RAM$rhs,collapse = "\n")
    
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
    RAM <- object@RAM
    
    # Remove intercepts:
    if (any(RAM$edge == "int"))
    {
      warning("Intercepts removed from model for 'sem' syntax")
      RAM <- RAM[RAM$edge!="int",]
    }
    
    RAM$label[RAM$fixed] <- NA
    ## Fix parameter labels.
    if (max(RAM$par) > 0)
    {
      for (i in seq_len(max(RAM$par)))
      {
        # Check if unique to other par numbers:
        if (any(RAM$label[RAM$par!=i] %in% RAM$label[RAM$par==i] | any(RAM$label[RAM$par == i] == '')))
        {
          RAM$label[RAM$par==i] <- paste0("par",i)
        }
        
        # Check if labels are unique, else combine:
        if (length(unique(RAM$label[RAM$par == i])) > 1)
        {
          RAM$label[RAM$par==i] <- paste(RAM$label[RAM$par==i],collapse="_")
        }
      }
    }
    
    # Fix estimate:
    RAM$est[!RAM$fixed] <- NA
    
    # Fix edges:
    RAM$edge[RAM$edge == '~>'] <- '->'
    
    # Create model:
    Mod <- paste(paste(RAM$lhs, RAM$edge, RAM$rhs), RAM$label, RAM$est, sep = ",", collapse = "\n")
    
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