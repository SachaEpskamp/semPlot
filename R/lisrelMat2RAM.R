modMat2Pars <- function(x,edge,exprname,symmetric=FALSE,vec=FALSE,cols,rows,group="",exprsup="")
{
  # Define x Pars:
  if (length(x)>0)
  {
    if (symmetric)
    {
      if (!isSymmetric(x$est)) stop(paste0("'",deparse(substitute(x)),"' matrix must be symmetrical."))
      x$est[upper.tri(x$est)] <- 0  
    }
    
    Pars <- data.frame(
      label = "", 
      lhs = rep(cols,each=length(rows)),
      edge = edge,
      rhs = rep(rows,times=length(cols)),
      est = c(x$est),
      std = NA,
      group = group,
      fixed = FALSE,
      par = 0,
      stringsAsFactors=FALSE)
    
    if (!vec)
    {
      Pars$label <- c(outer(1:nrow(x$est),1:ncol(x$est),function(x,y)paste0(exprname,"[",x,y,"]",exprsup)))
    } else {
      Pars$label <- paste0(exprname,"[",1:length(x$est),"]",exprsup)
    }
    if (!is.null(x[['std']]))
    {
      Pars[['std']] <- c(x[['std']])
    }
    if (!is.null(x[['par']]))
    {
      Pars[['par']] <- c(x[['par']])
    }
    if (!is.null(x[['fixed']]))
    {
      Pars[['fixed']] <- c(x[['fixed']])
    }
    
  } else Pars <- data.frame(
    label = character(0), 
    lhs = character(0),
    edge = character(0),
    rhs = character(0),
    est = numeric(0),
    std = numeric(0),
    group = character(0),
    fixed = logical(0),
    par = numeric(0),
    stringsAsFactors=FALSE)
  
  return(Pars)
}

