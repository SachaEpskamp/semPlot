semPlotModel.regsem <- function(object){

    
    ## First we divide the strings of the model's operations
 object1 <- object$lav.model@ParTable
 object2 <- object$lav.model@Model@dimNames
 varnames <- object2[[1]][2]
 mannames <- object2[[1]][1]
 names(varnames) <- 'name'
 names(mannames) <- 'manifest'
 
 
 #$lav.model@ParTable$
    
    ## Create a S4 list
    semModel <- new("semPlotModel")
    
    ## Create a Pars data frame
    semModel@Pars <- data.frame(
      label = rep("", length(object1$id)),
      lhs = object1$lhs,
      edge = object1$op,
      rhs = object1$rhs,
      est = object1$est,  # check if we should take estimates from other model, if estimates are same as in regsem
      std = NA,
      group = object1$group,
      fixed = object1$free == 0,
      par = object1$free,
      stringsAsFactors=FALSE)
    row.names(semModel@Pars) <- 1:length(object1$id)
    
    ## translate operators
    semModel@Pars$edge[object1$op=="~~"] <- "<->"  
    semModel@Pars$edge[object1$op=="~*~"] <- "<->"  
    semModel@Pars$edge[object1$op=="~"] <- "~>"
    semModel@Pars$edge[object1$op=="=~"] <- "->"
    semModel@Pars$edge[object1$op=="~1"] <- "int"
    semModel@Pars$edge[grepl("\\|",object1$op)] <- "|"
    
    semModel@Pars <- semModel@Pars[!object$op%in%c(':=','<','>','==','|','<', '>'),]
    
    ## Create a vars data frame
    '%!in%' <- function(x,y)!('%in%'(x,y))
    
    semModel@Vars <- data.frame(
      name = varnames$name[1:length(varnames$name)],
      manifest = varnames$name[1:length(varnames$name)] %in% mannames$manifest[1:length(mannames$manifest)],
      exogenous = NA,
      stringsAsFactors = FALSE
    )
    
    ## Miscellaneous data frames
    semModel@Thresholds <- data.frame()
    semModel@ObsCovs <- list()  
    semModel@ImpCovs <- list()
    semModel@Computed <- FALSE
    semModel@Original <- list(object)
    
    return(semModel)
  }