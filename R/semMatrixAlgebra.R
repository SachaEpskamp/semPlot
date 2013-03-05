semMatrixAlgebra <- function(object, algebra, group, simplify = TRUE)
{
  if ("lisrel"%in%class(object)) object <- object$matrices
  stopifnot("semMatrixModel"%in%class(object))
  
  if (missing(group)) group <- seq_len(max(sapply(object,length)))
  
  Mats <- lapply(object,lapply,'[[','est')
  Res <- list()
  for (i in seq_along(group))
  {
    GroupMats <- lapply(Mats,'[[',i) 
    Res[[i]] <- eval(substitute(algebra), GroupMats)
  }
  
  if (simplify) if (length(Res)==1) Res <- Res[[1]]
  return(Res)
}