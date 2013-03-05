Imin <- function(x,inverse=FALSE)
{
  if (any(dim(x)==0)) 
  {
    return(array(0,dim=dim(x))) 
  } else {
    x <- diag(1,nrow(x),ncol(x)) - x
    if (inverse) 
    {
      res <- tryCatch(solve(x), error = function(e) FALSE)
    } else {
      res <- x
    }
    if (is.matrix(res)) return(res) else 
    {
      warning("Uninvertable matrix found. Standardized solutions are not proper.")
      return(array(0, dim=dim(x)))
    }
  }
}
