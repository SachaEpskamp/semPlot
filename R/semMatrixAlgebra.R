semMatrixAlgebra <- function(object, algebra, group, simplify = TRUE, model)
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
  
  if ("lisrel"%in%class(object)) object <- object$matrices
  if (!"semMatrixModel"%in%class(object)) 
  {
    if (missing(model)) stop("'model' must be assigned if input is not semMatrixModel")
    object <- modelMatrices(object,model)
  }
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