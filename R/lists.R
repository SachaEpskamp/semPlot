
semPlotModel.list <- function(object)
{
  isModel <- sapply(object,function(x)"semPlotModel"%in%class(x))
  object[!isModel] <- lapply(object[!isModel],semPlotModel)
  if (length(object)>1) 
  {
    Res <- object[[1]]
    for (i in 2:length(object)) Res <- Res + object[[i]]
    return(Res) 
  } else return(object)
}