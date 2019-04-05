

semPlotModel.regsemplot <- function(object,...){
  
  ## Save parts of the output in objects 
  object1 <- object$lav.model@ParTable  # parameters 
  object2 <- object$lav.model@Model@dimNames  # variable names 
  varnames <- unique(c(object1$lhs, object1$rhs))  # all names
  mannames <- object2[[1]][1]  # manifest variables
  names(mannames) <- 'manifest'
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  ## Add the fixed relations to the parameter estimates of regsem 
  namelist <- strsplit(names(object$out$pars)," ")  # split names and operators
  inout <- data.frame(1,2)
  for(i in 1:length(namelist)){
    inout[i,1] <- namelist[[i]][1]
    inout[i,2] <- namelist[[i]][3]
  }  # create data frame of regsem variables
  
  int <- data.frame(1,2)
  for(i in 1:length(object1$lhs)){
    int[i,1] <- ifelse(object1$op[i]=="~"|object1$op[i]=="~1",object1$rhs[i],object1$lhs[i])   
    int[i,2] <- ifelse(object1$op[i]=="~"|object1$op[i]=="~1",object1$lhs[i],object1$rhs[i])  
  }  # create data frame of lavaan variables
  
  ## Paste together 
  pinout <- with(inout, paste0(X1, X2))
  pint <- with(int, paste0(X1, X2))
  counter <- 0
  
  for(i in 1:length(pint)){ 
    if(pint[i] %!in% pinout){
      object1$regest[i] <- 1
      counter <- counter + 1 
    } else {
      object1$regest[i] <- object$out$pars[i - counter]
    }
  }  # match regsem estimates with lavaan variables, set fixed to 1
  
  
  ## Create a S4 list
  semModel <- new("semPlotModel")
  
  ## Create a Pars data frame
  semModel@Pars <- data.frame(
    label = rep("", length(object1$id)),
    lhs = ifelse(object1$op=="~"|object1$op=="~1",object1$rhs,object1$lhs),  # first went from left to right without checking relationship 
    edge = "--",
    rhs = ifelse(object1$op=="~"|object1$op=="~1",object1$lhs,object1$rhs),
    est = object1$regest,  # check if we should take estimates from other model, if estimates are same as in regsem
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
  semModel@Vars <- data.frame(
    name = varnames,
    manifest = varnames[1:length(varnames)] %in% mannames$manifest[1:length(mannames$manifest)],
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
