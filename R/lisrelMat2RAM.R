lisrelMat2RAM <- function(x,edge,exprname,symmetric=FALSE,vec=FALSE,cols,rows,group="",exprsup="")
{
  # Define x RAM:
  if (length(x)>0)
  {
    if (symmetric)
    {
      if (!isSymmetric(x$est)) stop(paste0("'",deparse(substitute(x)),"' matrix must be symmetrical."))
      x$est[upper.tri(x$est)] <- 0  
    }
    
    RAM <- data.frame(
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
      RAM$label <- c(outer(1:nrow(x$est),1:ncol(x$est),function(x,y)paste0(exprname,"[",x,y,"]",exprsup)))
    } else {
      RAM$label <- paste0(exprname,"[",1:length(x$est),"]",exprsup)
    }
    if (!is.null(x[['std']]))
    {
      RAM[['std']] <- c(x[['std']])
    }
    if (!is.null(x[['par']]))
    {
      RAM[['par']] <- c(x[['par']])
    }
    if (!is.null(x[['fixed']]))
    {
      RAM[['fixed']] <- c(x[['fixed']])
    }
    
  } else RAM <- data.frame(
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
  
  return(RAM)
}

