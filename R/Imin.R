Imin <- function(x,inverse=FALSE)
{
  if (any(dim(x)==0)) 
  {
    return(array(0,dim=dim(x))) 
  } else {
    x <- diag(1,nrow(x),ncol(x)) - x
    if (inverse) 
    {
      res <- tryCatch(solve(x), error = function(e) FALSE, silent = TRUE)
      if (is.matrix(res)) return(res) else 
      {
        res <- tryCatch(pseudoinverse(x), error = function(e) FALSE, silent = TRUE)
        if (is.matrix(res))
        {
          warning("Psuedoinverse used for singular matrix. Standardized solution might not be proper.")
          return(res) 
        } else 
        {
          warning("Uninvertable matrix found and psuedoinverse could not be computed. Standardized solutions probably not proper.")
          return(array(0, dim=dim(x)))
        }
      }
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
