# Function to extract parameters into model matrices:
# Object is semPlotModel or can be created:

# Inner functions:
getRAMmodel <- function(object)
{
  
  # Parameters:
  Nvar <- nrow(object@Vars)
  Nman <- sum(object@Vars$manifest)
  Names <- object@Vars$name
  
  # Empty matrices:
  A <- S <- matrix(0,Nvar,Nvar)
  F <- cbind(diag(1,Nman),matrix(0,Nman,Nvar-Nman))
  F[,order(object@Vars$manifest,decreasing=TRUE)] <- F
  rownames(A) <- colnames(A) <- rownames(S) <- colnames(S) <- colnames(F) <- Names
  rownames(F) <- Names[object@Vars$manifest]  
  
  # Fill matrices:
  for (i in seq_len(nrow(object@Pars)))
  {
    if (object@Pars$edge[i]=="<->")
    {
      S[which(Names==object@Pars$lhs[i])[1],which(Names==object@Pars$rhs[i])[1]] <- 
        S[which(Names==object@Pars$rhs[i])[1],which(Names==object@Pars$lhs[i])[1]] <- object@Pars$est[i]
    }
    if (object@Pars$edge[i]%in%c("->","~>"))
    {
      A[which(Names==object@Pars$rhs[i])[1],which(Names==object@Pars$lhs[i])[1]] <- object@Pars$est[i]
    }
  }
  
  Res <- list(A=A,S=S,F=F)
  class(Res) <- "RAM"
  return(Res)  
}

modelMatrices <- function(object,model="ram")
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
  
  ### RAM MODEL ###
  
}