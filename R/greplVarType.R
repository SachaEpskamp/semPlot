# grepl on varnames with special keywords:
# - MAN
# - LAT
# - ENDO
# - EXO
# - INT

matchVar <- function(x, Vars, manIntsExo, manIntsEndo, latIntsExo, latIntsEndo)
{
  
  n <- nrow(Vars) + nrow(manIntsEndo)  + nrow(manIntsExo)  + nrow(latIntsEndo)  + nrow(latIntsExo) 
  
  Man <- c(Vars$manifest, rep(FALSE,n-nrow(Vars)))
  Man[c(manIntsEndo[,1],manIntsExo[,1])] <- TRUE
  
  Exo <- c(Vars$exogenous, rep(FALSE,n-nrow(Vars)))
  Exo[c(manIntsExo[,1],latIntsExo[,1])] <- TRUE
  
  isInt <- c(rep(FALSE,nrow(Vars)), rep(TRUE, n-nrow(Vars)))
  
  # match:
  matchRes <- match(x,Vars$name)
  matchRes <- matchRes[!is.na(matchRes)]
  
  # keywords:
  select <- rep(grepl("(EXO)|(ENDO)|(MAN)|(LAT)|(INT)|(VAR)",x),n)
  
  if (any(select))
  {
    
    if (grepl("(ENDO)|(EXO)",x))
    {      
      # First node first / endo:
      select <- select & ((grepl("ENDO",x) & !Exo) | 
                            (grepl("EXO",x) & Exo )
      )
    }
    
    if (grepl("(LAT)|(MAN)",x))
    {
      
      # Any node man / latent
      select <- select & ((grepl("LAT",x) & !Man) | 
                            (grepl("MAN",x) & Man )
      )
    }
    
    if (grepl("(INT)|(VAR)",x))
    {
      
      # Any node man / latent
      select <- select & ((grepl("VAR",x) & !isInt) | 
                            (grepl("INT",x) & isInt )
      )
    }
     
  }
  
  return(c(matchRes,which(select)))
}