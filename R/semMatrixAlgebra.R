semMatrixAlgebra <- function(object, algebra, group, simplify = TRUE, model, endoOnly = FALSE)
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
    if (missing(model))
    {
      if (any(grepl("(LY)|(TE)|(PS)|(BE)|(LX)|(TD)|(PH)|(GA)|(TY)|(TX)|(AL)|(KA)",deparse(substitute(algebra)))))      {
        model <- "lisrel"
        message("model set to 'lisrel'")
      } else if (any(grepl("(Lambda)|(Nu)|(Theta)|(Kappa)|(Alpha)|(Beta)|(Gamma)|(Psi)",deparse(substitute(algebra)))))
      {
        model <- "mplus"
        message("model set to 'mplus'")
      } else if (any(grepl("A|S|F",deparse(substitute(algebra)))))
      {
        model <- "ram"
        message("model set to 'ram'")
      } else stop("'model' could not be detected")
    } 
    object <- modelMatrices(object,model,endoOnly = endoOnly)
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