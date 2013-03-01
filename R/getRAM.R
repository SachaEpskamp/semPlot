getRAM <- function(object)
{
  
  # Parameters:
  Nvar <- nrow(object@Vars)
  Nman <- sum(object@Vars$manifest)
  Names <- object@Vars$name
  
  # Empty matrices:
  A <- S <- matrix(0,Nvar,Nvar)
  F <- matrix(0,Nman,Nvar)
  rownames(A) <- colnames(A) <- rownames(S) <- colnames(S) <- colnames(F) <- Names
  rownames(F) <- Names[object@Vars$manifest]  
  
  
  
}