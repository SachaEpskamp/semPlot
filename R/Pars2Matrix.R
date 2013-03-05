# Inner function, computes matrix from Pars subsection:
# Pars: Sub of Pars
# rows: Rownames
# cols: Colnames
# lhsisrow: lhs variable is interpreted as row (default to FALSE)
Pars2Matrix <- function(Pars, edges, rows, cols, symmetrical, lhsisrow = FALSE)
{
  if (missing(symmetrical))
  {
    symmetrical <- any(grepl("<->",edges))
  }
  if (lhsisrow) Pars[c('lhs','rhs')] <- Pars[c('rhs','lhs')]
  Groups <- unique(Pars$group)
  Pars$lhs[Pars$edge=="int"] <- "1"
  
  Pars <- Pars[Pars$edge %in% edges & Pars$lhs %in% cols & Pars$rhs %in% rows,]
  
  ResMatrix <- list()
  empMatrix <- matrix(0, length(rows), length(cols))
  rownames(empMatrix) <- gsub("@L@","",rows)
  colnames(empMatrix) <- gsub("@L@","",cols)
  for (i in seq_along(Groups))
  {
    ResMatrix[[i]] <- list()
    ResMatrix[[i]]$est <- empMatrix
    ResMatrix[[i]]$std <- empMatrix
    ResMatrix[[i]]$par <- empMatrix
    ResMatrix[[i]]$fixed <- empMatrix
    mode(ResMatrix[[i]]$fixed) <- "logical"
    for (j in seq_len(nrow(Pars)))
    {
      ResMatrix[[i]]$est[match(Pars$rhs[j],rows),match(Pars$lhs[j],cols)] <- Pars$est[j]
      ResMatrix[[i]]$std[match(Pars$rhs[j],rows),match(Pars$lhs[j],cols)] <- Pars$std[j]
      ResMatrix[[i]]$fixed[match(Pars$rhs[j],rows),match(Pars$lhs[j],cols)] <- Pars$fixed[j]
      ResMatrix[[i]]$par[match(Pars$rhs[j],rows),match(Pars$lhs[j],cols)] <- Pars$par[j]
      if (symmetrical)
      {
        ResMatrix[[i]]$est[match(Pars$lhs[j],rows),match(Pars$rhs[j],cols)] <- Pars$est[j]
        ResMatrix[[i]]$std[match(Pars$lhs[j],rows),match(Pars$rhs[j],cols)] <- Pars$std[j]
        ResMatrix[[i]]$fixed[match(Pars$lhs[j],rows),match(Pars$rhs[j],cols)] <- Pars$fixed[j]
        ResMatrix[[i]]$par[match(Pars$lhs[j],rows),match(Pars$rhs[j],cols)] <- Pars$par[j]        
      }
    }
  }
  names(ResMatrix) <- Groups
  return(ResMatrix)
}


FilterMatrix <- function(Pars, Vars)
{
  Groups <- unique(Pars$group)

  ResMatrix <- list()
  Nvar <- nrow(Vars)
  Nman <- sum(Vars$manifest)
  
  for (i in seq_along(Groups))
  {
    ResMatrix[[i]] <- list()
    ResMatrix[[i]]$est <- cbind(diag(1,Nman),matrix(0,Nman,Nvar-Nman))
    ResMatrix[[i]]$est[,order(Vars$manifest,decreasing=TRUE)] <- ResMatrix[[i]]$est
    rownames(ResMatrix[[i]]$est) <- Vars$name[Vars$manifest]
    colnames(ResMatrix[[i]]$est) <- Vars$name
  }
  
  names(ResMatrix) <- Groups
  return(ResMatrix)
}