
semPlotModel.list <- function(object, ...)
{
  if ("mplus.model"%in%class(object)) return(semPlotModel.mplus.model(object))
  
  mod <- try(semPlotModel_lavaanModel(object,...),silent=TRUE)
  if (!"try-error"%in%class(mod)) return(mod)
  
  isModel <- sapply(object,function(x)"semPlotModel"%in%class(x))
  object[!isModel] <- lapply(object[!isModel],semPlotModel)
  if (length(object)>1) 
  {
    Res <- object[[1]]
    for (i in 2:length(object)) Res <- Res + object[[i]]
    return(Res) 
  } else return(object)
}