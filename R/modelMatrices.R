# Function to extract parameters into model matrices:
# Object is semPlotModel or can be created:

modelMatrices <- function(object,model="ram")
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
  
  if (!"semPlotModel"%in%class(object)) object <- semPlotModel(object)
  stopifnot("semPlotModel"%in%class(object))
  
  ### RAM MODEL ###
  
}