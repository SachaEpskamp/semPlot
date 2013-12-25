
## Mode function:
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Function to scale and rotate layouts:
LayoutScaler <- function(x, xrange=1, yrange=1)
{
  if ((max(x[,1]) - min(x[,1])) == 0) x[,1] <- mean(xrange) else x[,1] <- (x[,1] - min(x[,1])) / (max(x[,1]) - min(x[,1])) * 2 - 1
  if ((max(x[,2]) - min(x[,2])) == 0) x[,2] <- mean(yrange) else x[,2] <- (x[,2] - min(x[,2])) / (max(x[,2]) - min(x[,2])) * 2 - 1
  
  x[,1] <- x[,1] * xrange
  x[,2] <- x[,2] * yrange
  
  return(x)
}

# Rotation function:
RotMat <- function(d,w2hrat=1) 
{
  matrix(c(cos(-d),sin(-d),-sin(-d),cos(-d)),2,2)
}


## Function to compute reingold-tilford layout using igraph:
rtLayout <- function(roots,GroupPars,Edgelist,layout,exoMan)
{
  # Reverse intercepts in graph:
  #   revNodes <- which((GroupPars$edge == "int" | Edgelist[,2] %in% exoMan) & !Edgelist[,1] %in% roots )
  #   revNodes <- which((GroupPars$edge == "int" & !Edgelist[,1] %in% roots) | Edgelist[,2] %in% exoMan )
  #   Edgelist[revNodes,1:2] <- Edgelist[revNodes,2:1]
  # Remove double headed arrows:
  Edgelist <- Edgelist[GroupPars$edge != "<->",]
  
  # Make igraph object:
  Graph <- graph.edgelist(Edgelist, FALSE)
  # Compute layout:
  Layout <- layout.reingold.tilford(Graph,root=roots,circular = FALSE) 
  
  return(Layout)
}

## Function to mix color vector x with weight w
mixColfun <- function(x,w)
{
  # x = vector of colors
  # w = weights
  if (missing(w)) w <- rep(1,length(x))
  if (length(w)==1) w <- rep(w,length(x))
  
  RGB <- col2rgb(x)
  wMeans <- apply(RGB,1,weighted.mean,w=w)
  return(rgb(wMeans[1],wMeans[2],wMeans[3],maxColorValue=255))
}

loopOptim <- function(x,Degrees)
{
  NotinRange <- sum(sapply(Degrees,function(d)!any(c(d,d-2*pi,d+2*pi)>(x-pi/4) & c(d,d-2*pi,d+2*pi)<(x+pi/4))))
  Dist2Edges <- sapply(Degrees,function(d)min(abs(x - c(d,d-2*pi,d+2*pi))))
  return(NotinRange * 2 * pi * 2 + sum(sort(Dist2Edges)[1:2]))
}

# RotMat <- function(d) matrix(c(cos(-d),sin(-d),-sin(-d),cos(-d)),2,2)

mixInts <- function(vars,intMap,Layout,trim=FALSE,intAtSide=TRUE)
{
  n <- length(vars)
  
  if (intAtSide)
  {
    if (!trim)
    {
      if (n+nrow(intMap)==1)
      {
        sq <- 0
      }
      if (n+nrow(intMap) == 2)
      {
        sq <- c(0,0.5) 
      } else {
        sq <- seq(-1,1,length=n+nrow(intMap))
      }
    } else {
      if (n+nrow(intMap) == 2)
      {
        sq <- c(0,0.5) 
      } else {
        sq <- seq(-1,1,length=n+nrow(intMap)+2)[-c(1,n+nrow(intMap)+2)]
      }
    }
    cent <- median(1:n)
    c <- 1
    for (i in seq_along(vars))
    {
      if (vars[i]%in%intMap[,2])
      {
        if (i < cent)
        {
          Layout[intMap[intMap[,2]==vars[i],1],1] <- sq[c]
          Layout[vars[i],1] <- sq[c+1]
          c <- c+2
        } else
        {
          Layout[intMap[intMap[,2]==vars[i],1],1] <- sq[c+1]
          Layout[vars[i],1] <- sq[c]
          c <- c+2                   
        }
      } else
      {
        Layout[vars[i],1] <- sq[c]
        c <- c+1
      }
    }
  } else {
    if (!trim)
    {
      if (n==1)
      {
        sq <- 0
      } else if (n == 2)
      {
        sq <- c(-1,1) 
      } else {
        sq <- seq(-1,1,length=n)
      }
    } else {
      if (n == 1)
      {
        sq <- 0
      } else if (n == 2)
      {
        sq <- c(-0.5,0.5) 
      } else {
        sq <- seq(-1,1,length=n+2)[-c(1,n+2)]
      }
    }
    c <- 1
    for (i in seq_along(vars))
    {
      if (vars[i]%in%intMap[,2])
      {
        Layout[intMap[intMap[,2]==vars[i],1],1] <- sq[c]
        Layout[vars[i],1] <- sq[c]
        c <- c + 1 
      } else
      {
        Layout[vars[i],1] <- sq[c]
        c <- c+1
      }
    }    
  }
  return(Layout)
}
