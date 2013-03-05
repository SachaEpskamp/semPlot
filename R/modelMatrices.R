# Function to extract parameters into model matrices:
# Object is semPlotModel or can be created:

# # Inner functions:
# getRAMmodel <- function(object)
# {
#   
#   # Parameters:
#   Nvar <- nrow(object@Vars)
#   Nman <- sum(object@Vars$manifest)
#   Names <- object@Vars$name
#   
#   # Empty matrices:
#   A <- S <- matrix(0,Nvar,Nvar)
#   F <- cbind(diag(1,Nman),matrix(0,Nman,Nvar-Nman))
#   F[,order(object@Vars$manifest,decreasing=TRUE)] <- F
#   rownames(A) <- colnames(A) <- rownames(S) <- colnames(S) <- colnames(F) <- Names
#   rownames(F) <- Names[object@Vars$manifest]  
#   
#   # Fill matrices:
#   for (i in seq_len(nrow(object@Pars)))
#   {
#     if (object@Pars$edge[i]=="<->")
#     {
#       S[which(Names==object@Pars$lhs[i])[1],which(Names==object@Pars$rhs[i])[1]] <- 
#         S[which(Names==object@Pars$rhs[i])[1],which(Names==object@Pars$lhs[i])[1]] <- object@Pars$est[i]
#     }
#     if (object@Pars$edge[i]%in%c("->","~>"))
#     {
#       A[which(Names==object@Pars$rhs[i])[1],which(Names==object@Pars$lhs[i])[1]] <- object@Pars$est[i]
#     }
#   }
#   
#   Res <- list(A=A,S=S,F=F)
#   class(Res) <- "RAM"
#   return(Res)  
# }

modelMatrices <- function(object,model="ram", endoOnly = FALSE)
{
  # Check if input is combination of models:
  call <- paste(deparse(substitute(object)), collapse = "")
  if (grepl("\\+",call)) 
  {
    args <- unlist(strsplit(call,split="\\+"))
    obs <- lapply(args,function(x)semPlotModel(eval(parse(text=x))))
    object <- obs[[1]]
    for (i in 2:length(obs)) object <- object + obs[[i]]
  }
  
  if (!"semPlotModel"%in%class(object)) object <- semPlotModel(object)
  stopifnot("semPlotModel"%in%class(object))
  
  ### SETUP ###
  Model <- list()
  class(Model) <- "semMatrixModel"
  
  # Define exogeneity:
  if (endoOnly)
  {
    object@Vars$exogenous <- FALSE  
  } else {
    if (any(is.na(object@Vars$exogenous)))
    {
      object <- defExo(object)
    }
  }
  
  
  ### RAM MODEL ###
  if (grepl("ram",model,ignore.case=TRUE))
  {
    # Extract names:
    man <- object@Vars$name[object@Vars$manifest]
    lat <- object@Vars$name[!object@Vars$manifest]
    all <- object@Vars$name
    
    # Extract matrices:
    Model[['A']] <- Pars2Matrix(object@Pars, c("->","~>"), all, all)
    Model[['S']] <- Pars2Matrix(object@Pars, "<->", all, all)
    Model[['F']] <- FilterMatrix(object@Pars, object@Vars)

    return(Model)
  }  
  
  
  ### LISREL MODEL ###:
  if (grepl("lis",model,ignore.case=TRUE))
  {
    # Extract names:
    manExo <- object@Vars$name[object@Vars$manifest & object@Vars$exogenous]
    manEndo <- object@Vars$name[object@Vars$manifest & !object@Vars$exogenous]
    latExo <- object@Vars$name[!object@Vars$manifest & object@Vars$exogenous]
    latEndo <- object@Vars$name[!object@Vars$manifest & !object@Vars$exogenous]
    
    # If any manifest var is used in regression, create dummy latents:
    if (any(object@Pars$lhs[object@Pars$edge%in%c("->","~>")] %in% c(manExo,manEndo)))
    {
      message("Latent dummy variables added to include manifest regressions")
      # Identify variables:
      manRegs <- c(manExo,manEndo)[c(manExo,manEndo)%in%object@Pars$lhs[object@Pars$edge%in%c("->","~>")]]
      newVars <- object@Vars[object@Vars$name %in% manRegs,]
      newVars$manifest <- FALSE
      newVars$name <- paste0(newVars$name,"@L@")
      object@Vars <- rbind(object@Vars,newVars)
      
      # Change regressions to latents:
      object@Pars$lhs[object@Pars$lhs %in% manRegs & object@Pars$edge%in%c("->","~>")] <- paste0(object@Pars$lhs[object@Pars$lhs %in% manRegs & object@Pars$edge%in%c("->","~>")],"@L@")
      
      manVarResids <- which(object@Pars$lhs %in% manRegs & object@Pars$rhs %in% manRegs & object@Pars$edge=="<->")
      
      object@Pars$lhs[manVarResids] <- paste0(object@Pars$lhs[manVarResids],"@L@")
      object@Pars$rhs[manVarResids] <- paste0(object@Pars$rhs[manVarResids],"@L@")
      
      
      # Add factor loadings:
      for (g in unique(object@Pars$group))
      {
        parLocs <- nrow(object@Pars)+seq_along(manRegs)
        object@Pars[parLocs,"lhs"] <- paste0(manRegs,"@L@")
        object@Pars[parLocs,"rhs"] <- manRegs
        object@Pars[parLocs,"label"] <- ""
        object@Pars[parLocs,"est"] <- 1
        object@Pars[parLocs,"std"] <- NA
        object@Pars[parLocs,"group"] <- g
        object@Pars[parLocs,"fixed"] <- TRUE
        object@Pars[parLocs,"par"] <- 0
      }
      
      # Extract names:
      manExo <- object@Vars$name[object@Vars$manifest & object@Vars$exogenous]
      manEndo <- object@Vars$name[object@Vars$manifest & !object@Vars$exogenous]
      latExo <- object@Vars$name[!object@Vars$manifest & object@Vars$exogenous]
      latEndo <- object@Vars$name[!object@Vars$manifest & !object@Vars$exogenous]
    }
    
    # Extract matrices:
    Model[['LY']] <- Pars2Matrix(object@Pars, c("->","~>"), manEndo, latEndo)
    Model[['TE']] <- Pars2Matrix(object@Pars, "<->", manEndo, manEndo)
    Model[['PS']] <- Pars2Matrix(object@Pars, "<->", latEndo, latEndo)
    Model[['BE']] <- Pars2Matrix(object@Pars, c("->","~>"), latEndo, latEndo)
    
    Model[['LX']] <- Pars2Matrix(object@Pars, c("->","~>"), manExo, latExo)
    Model[['TD']] <- Pars2Matrix(object@Pars, "<->", manExo, manExo)
    Model[['PH']] <- Pars2Matrix(object@Pars, "<->", latExo, latExo)
    Model[['GA']] <- Pars2Matrix(object@Pars, c("->","~>"), latEndo, latExo)
    
    Model[['TY']] <- Pars2Matrix(object@Pars, "int", manEndo, "1")
    Model[['TX']] <- Pars2Matrix(object@Pars, "int", manExo, "1")
    Model[['AL']] <- Pars2Matrix(object@Pars, "int", latEndo, "1")
    Model[['KA']] <- Pars2Matrix(object@Pars, "int", latExo, "1")
    
    return(Model)
  }
  
  
  
  ### Mplus MODEL ###:
  if (grepl("mplus",model,ignore.case=TRUE))
  {
    # Extract names (exo only if manifest has outgoing cons. error if in and outgoing):
    man <- object@Vars$name[object@Vars$manifest]
    lat <- object@Vars$name[!object@Vars$manifest]
    
    # Control input:
    if (any(sapply(man, function(m)  any((object@Pars$lhs==m & object@Pars$edge %in% c("->","~>")) & (object@Pars$rhs==m & object@Pars$edge %in% c("->","~>")))))) stop("Manifest variable found with both incoming and outgoing edge. This is not yet supported in modelMatrices.")
    if (any(object@Pars$rhs %in% man & object@Pars$lhs %in% lat & object@Pars$edge == "~>"))
    {
      warning("Can not place regression (ON) from latent to manifest in a model matrix. Interpreted as factor loading (BY).")
      object@Pars$edge[object@Pars$rhs %in% man & object@Pars$lhs %in% lat & object@Pars$edge == "~>"] <- "->"
    }
    
    trueExo <- sapply(man, function(m)  any((object@Pars$lhs==m & object@Pars$edge %in% c("->","~>")) & !(object@Pars$rhs==m & object@Pars$edge %in% c("->","~>"))))
    manEndo <- man[!trueExo]
    manExo <- man[trueExo]
    
    
    ## Extract matrices:
    # BY matrices:
    Model[['Nu']] <- Pars2Matrix(object@Pars, "int", manEndo, "1")
    Model[['Lambda']] <- Pars2Matrix(object@Pars, c("->","~>"), manEndo, lat)
    Model[['Theta']] <- Pars2Matrix(object@Pars, "<->", manEndo, manEndo)
    
    # ON matrices:
    Model[['Kappa']] <- Pars2Matrix(object@Pars, c("->","~>"), manEndo, manExo)
    Model[['Alpha']] <- Pars2Matrix(object@Pars, "int", lat, "1")
    Model[['Beta']] <- Pars2Matrix(object@Pars, c("->","~>"), lat, lat)    
    Model[['Gamma']] <- Pars2Matrix(object@Pars, c("->","~>"), lat, manExo)
    Model[['Psi']] <- Pars2Matrix(object@Pars, "<->", lat, lat)

    return(Model)
  }
  
  
  else stop(paste("Model",model,"is not supported."))
  
}