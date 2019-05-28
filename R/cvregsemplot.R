

semPlotModel.cvregsem <- function(object,model,...){
  if (missing(model)){
    stop("Please supply lavaan model with 'model' argument!")
  }
  ## Save parts of the output in objects 
  object1 <- object  # parameters
  object2 <- model@ParTable  # lavaan parameters
  varnames <- unique(c(object2$lhs, object2$rhs))  # all names
  mannames <- model@Model@dimNames[[1]][1]  # manifest variables
  names(varnames) <- 'name'
  names(mannames) <- 'manifest'
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  ## Add the fixed relations to the parameter estimates of regsem 
  namelist <- strsplit(names(object1$final_pars)," ")   # split names and operators
  inout <- data.frame(1,2)
  for(i in 1:length(namelist)){
    inout[i,1] <- namelist[[i]][1]
    inout[i,2] <- namelist[[i]][3]
  }  # create data frame of regsem variables
  
  int <- data.frame(1,2)
  for(i in 1:length(object2$lhs)){
    int[i,1] <- ifelse(object2$op[i]=="~"|object2$op[i]=="~1",object2$rhs[i],object2$lhs[i])   
    int[i,2] <- ifelse(object2$op[i]=="~"|object2$op[i]=="~1",object2$lhs[i],object2$rhs[i])  
  }  # create data frame of lavaan variables
 
  ## paste together 
  pinout <- with(inout, paste0(X1, X2))
  pint <- with(int, paste0(X1, X2))
  
  counter <- 0
  for(i in 1:length(object2$free)){  # if free before, 
    if(object2$free[i] == 0){
      object1$regest[i] <- 1
      counter = counter + 1
    } else{
      object1$regest[i] <- object1$final_pars[i - counter]
    }
  }  # match regsem estimates with lavaan variables, set fixed to 1
  
  
  ## Create a S4 list
  semModel <- new("semPlotModel")
  
  ## Create a Pars data frame
  semModel@Pars <- data.frame(
    label = rep("", length(object2$id)),
    lhs = ifelse(object2$op=="~"|object2$op=="~1",object2$rhs,object2$lhs),  # first went from left to right without checking relationship 
    edge = "--",
    rhs = ifelse(object2$op=="~"|object2$op=="~1",object2$lhs,object2$rhs),
    est = object1$regest,  # check if we should take estimates from other model, if estimates are same as in regsem
    std = NA,
    group = object2$group,
    fixed = object2$free == 0,
    par = object2$free,
    stringsAsFactors=FALSE)
  row.names(semModel@Pars) <- 1:length(object2$id)
  
  ## translate operators
  semModel@Pars$edge[object2$op=="~~"] <- "<->"  
  semModel@Pars$edge[object2$op=="~*~"] <- "<->"  
  semModel@Pars$edge[object2$op=="~"] <- "~>"
  semModel@Pars$edge[object2$op=="=~"] <- "->"
  semModel@Pars$edge[object2$op=="~1"] <- "int"
  semModel@Pars$edge[grepl("\\|",object2$op)] <- "|"
  
  semModel@Pars <- semModel@Pars[!object2$op%in%c(':=','<','>','==','|','<', '>'),]
  
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
