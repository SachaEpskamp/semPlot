# Arguments:
# rotation: 1 = normal (endo manifests under), 2 3 and 4 flip counterclockwise.
# 2 = endo manifests righy
# 3 = endo man up
# 4 = endo man left
# allVars: TRUE includes variables that are not in model (e.g. with between-within group models)

# Layout modes:
# "tree"
# "circle"
# "spring"
# igraph function
# "tree2" and  "circle2" for layout.reingold.tilford

# manifests: vector of manifest labels ordered
# latents: vector of latents ordered
# fixedStyle: if coercible to numeric lty is assigned this value, else a color for color representation. If this argument is not a number or color representation the edge is not displayed differently.


## Mode function:
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## Function to compute reingold-tilford layout using igraph:
rtLayout <- function(roots,GroupRAM,Edgelist,layout,exoMan)
{
  # Reverse intercepts in graph:
#   revNodes <- which((GroupRAM$edge == "int" | Edgelist[,2] %in% exoMan) & !Edgelist[,1] %in% roots )
  revNodes <- which((GroupRAM$edge == "int" & !Edgelist[,1] %in% roots) | Edgelist[,2] %in% exoMan )
  Edgelist[revNodes,1:2] <- Edgelist[revNodes,2:1]
  # Remove double headed arrows:
  Edgelist <- Edgelist[GroupRAM$edge != "<->",]
  
  # Make igraph object:
  Graph <- graph.edgelist(Edgelist, TRUE)
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

RotMat <- function(d) matrix(c(cos(-d),sin(-d),-sin(-d),cos(-d)),2,2)

mixInts <- function(vars,intMap,Layout,trim=FALSE,residuals=TRUE)
{
  n <- length(vars)
  
  if (residuals)
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
      }
      if (n == 2)
      {
        sq <- c(-1,1) 
      } else {
        sq <- seq(-1,1,length=n)
      }
    } else {
      if (n == 1)
      {
        sq <- 0
      }
      if (n == 2)
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

semPaths <- function(object,what="paths",whatLabels,style,layout="tree",intercepts=TRUE,residuals=TRUE,thresholds=TRUE,intStyle="multi",rotation=1,curve,nCharNodes=3,nCharEdges=3,sizeMan = 5,sizeLat = 8,sizeInt = 2,ask,mar,title,title.color="black",include,combineGroups=FALSE,manifests,latents,groups,color,residScale,gui=FALSE,allVars=FALSE,edge.color,reorder=TRUE,structural=FALSE,ThreshAtSide=FALSE,threshold.color,fixedStyle=2,freeStyle=1,as.expression=character(0),optimizeLatRes=FALSE,mixCols=TRUE,curvePivot,...){
  
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
  
  # if (gui) return(do.call(semPathsGUI,as.list(match.call())[-1]))
  
  # Check:
  if (!rotation%in%1:4)
  {
    stop("Rotation must be 1, 2 3 or 4.")
  }
  if (any(object@RAM$edge=="int")) 
  {
    object@Vars$name[object@Vars$name=="1"] <- "_1"
    object@RAM$lhs[object@RAM$lhs=="1"] <- "_1"
    object@RAM$rhs[object@RAM$rhs=="1"] <- "_1"
  }
  
  # Check if layout is function. If so, set layout <- "spring" as dummy:
  if (is.function(layout))
  {
    layoutFun <- layout
    layout <- "spring"
  } else layoutFun <- NULL
  
  
  if (missing(curve))
  {
    if (layout == "tree")
    {
      curve <- 1
    } else {
      curve <- 0
    }
  }
  curveDefault <- curve
  
  if (missing(curvePivot))
  {
    curvePivot <- grepl("tree",layout)
  }
  
  if (missing(whatLabels))
  {
    edge.labels <- TRUE
  } else
  {
    edge.labels <- FALSE    
  }
  
  if (missing(as.expression)) as.expression <- ""
  
  # Set and check style: 
  if (missing(style)) style <- "OpenMx"
  if (!grepl("mx|lisrel",style,ignore.case=TRUE)) stop("Only OpenMx or LISREL style is currently supported.")
  #   if (grepl("mx",style,ignore.case=TRUE) & !missing(residScale)) warning("'residScale' ingored in OpenMx style")
  if (missing(residScale)) residScale <- 2*sizeMan
  
  #   residScale <- residScale * 1.75
  # Remove means if means==FALSE
  if (intercepts==FALSE)
  {
    object@RAM <- object@RAM[object@RAM$edge!="int",]
  }
  # Remove residuals if residuals=FALSE
  if (residuals==FALSE)
  {
    object@RAM <- object@RAM[!(object@RAM$edge=="<->"&object@RAM$lhs==object@RAM$rhs),]
  }  
  # Combine groups if combineGroups=TRUE:
  if (combineGroups)
  {
    object@RAM$group <- ""
  }
  # Set title:
  if (missing(title))
  {
    # Check titles:
    title <- length(unique(object@RAM$group))>1
  }
  # If structural, remove all manifest from RAM:
  if (structural)
  {
    #     browser()
    object@RAM <- object@RAM[!(object@RAM$lhs %in% object@Vars$name[object@Vars$manifest] | object@RAM$rhs %in% object@Vars$name[object@Vars$manifest]),]
    object@Vars <- object@Vars[!object@Vars$manifest,]
    object@Thresholds <- data.frame()
  }  
  
  # Add rows for bidirectional edges:
  if (any(object@RAM$edge=="<->" & object@RAM$lhs != object@RAM$rhs))
  {
    bidirs <- object@RAM[object@RAM$edge=="<->" & object@RAM$lhs != object@RAM$rhs,]
    bidirs[c("lhs","rhs")] <- bidirs[c("rhs","lhs")]
    bidirs$par <- -1
    object@RAM <- rbind(object@RAM,bidirs)
  }
  object@RAM <- object@RAM[!duplicated(object@RAM),]
  
  # Extract names:
  manNames <- object@Vars$name[object@Vars$manifest]
  if (!missing(manifests))
  {
    if (!(all(manNames%in%manifests) & length(manifests) == length(manNames)))
    {
      stop(paste("Argument 'manifests' should be a vector containing reordered elements of the vector",dput(manNames)))
    }
    manNames <- manifests
  }
  latNames <- object@Vars$name[!object@Vars$manifest]
  if (!missing(latents))
  {
    if (!(all(latNames%in%latents) & length(latents) == length(latNames)))
    {
      stop(paste("Argument 'latents' should be a vector containing reordered elements of the vector",dput(latNames)))
    }
    latNames <- latents
  }
  Labels <- c(manNames,latNames)
  nM <- length(manNames)
  nL <- length(latNames)
  nN <- length(Labels)
  object@Vars <- object@Vars[match(Labels,object@Vars$name),]
  
  # Define groups and colors setup:
  DefaultColor <- FALSE
  if (!missing(groups))
  {
    if (is.character(groups))
    {
      if (any(grepl("man",groups,ignore.case=TRUE)) & any(grepl("lat",groups,ignore.case=TRUE)))
      {
        groups <- as.list(c(manNames,latNames))
      } else if (any(grepl("man",groups,ignore.case=TRUE)) & !any(grepl("lat",groups,ignore.case=TRUE)))
      {
        groups <- as.list(manNames)
      } else if (!any(grepl("man",groups,ignore.case=TRUE)) & any(grepl("lat",groups,ignore.case=TRUE)))
      {
        groups <- as.list(latNames)
      } else stop("Character specification of 'groups' must contain 'man','lat' or both")
    }
    if (is.factor(groups) | is.character(groups)) groups <- tapply(1:length(groups),groups,identity)
    
    if (!is.list(groups)) stop("'groups' argument is not a factor or list")
    
    if (missing(color)) 
    {
      color <- rainbow(length(groups))
    }
  } else {
    if (missing(color)) 
    {
      color <- "background"
    } 
  }
  
  # Define exogenous variables (only if any is NA):
  if (any(is.na(object@Vars$exogenous)))
  {
    if (any(!is.na(object@Vars$exogenous)))
    {
      exoOrig <- object@Vars$exogenous
      repExo <- TRUE
    } else repExo <- FALSE
    object@Vars$exogenous <- FALSE
    for (i in which(!object@Vars$manifest))
    {
      if (!any(object@RAM$edge[object@RAM$rhs==object@Vars$name[i]] %in% c("~>","->") & object@RAM$lhs[object@RAM$rhs==object@Vars$name[i]]%in%latNames))
      {
        object@Vars$exogenous[i] <- TRUE
      }
    }
    for (i in which(object@Vars$manifest))
    {
      if (all(object@RAM$lhs[object@RAM$rhs==object@Vars$name[i] & object@RAM$lhs%in%latNames]%in%object@Vars$name[object@Vars$exogenous]) &
            all(object@RAM$rhs[object@RAM$lhs==object@Vars$name[i] & object@RAM$rhs%in%latNames]%in%object@Vars$name[object@Vars$exogenous]) &
            !any(object@RAM$rhs==object@Vars$name[i] & object@RAM$edge=="~>"))
      {
        object@Vars$exogenous[i] <- TRUE
      }
    }
    
    # If all exo, treat all as endo:
    if (all(object@Vars$exogenous) | layout%in%c("circle","circle2"))
    {
      object@Vars$exogenous <- FALSE
    }
    # If al endo, treat formative manifest as exo (MIMIC mode), unless all manifest are formative.
    if (!any(object@Vars$exogenous))
    {
      if (any(object@Vars$manifest & (object@Vars$name%in%object@RAM$rhs[object@RAM$edge %in% c("~>","--","->")])))
        object@Vars$exogenous[object@Vars$manifest & !(object@Vars$name%in%object@RAM$rhs[object@RAM$edge %in% c("~>","--","->")])] <- TRUE
    }
    if (repExo)
    {
      object@Vars$exogenous[!is.na(exoOrig)] <- exoOrig[!is.na(exoOrig)]
    }
  }
  
  Groups <- unique(object@RAM$group)
  qgraphRes <- list()
  if (missing(ask))
  {
    if (length(Groups)>1) ask <- TRUE else ask <- FALSE
  }
  askOrig <- par("ask")
  
  if (missing(include)) include <- 1:length(Groups)
  
  # Reassign labels (temporary solution for excluding vars in multi group)
  AllLabs <- Labels
  AllMan <- manNames
  AllLat <- latNames
  
  par(ask=ask)
  for (gr in Groups[(1:length(Groups))%in%include])
  {
    GroupRAM <- object@RAM[object@RAM$group==gr,]
    GroupVars <- object@Vars
    GroupThresh <- object@Thresholds[object@Thresholds$group==gr,]
    
    # Restore Labels, manNames and latNames:
    Labels <- AllLabs
    manNames <- AllMan
    latNames <- AllLat
    nM <- length(AllMan)
    nL <- length(AllLat)
    
    ### Reorder nodes in order of factors
    if (reorder)
    {
      ConOrd <- function(nodes)
      {
        E <- GroupRAM[c("lhs","rhs")]
        subE <- rbind(as.matrix(E[E[,1]%in%nodes,1:2]),as.matrix(E[E[,2]%in%nodes,2:1]))
        subE <- subE[subE[,2]%in%GroupVars$name[!GroupVars$manifest],,drop=FALSE]
        ranks <- rank(match(subE[,2],GroupVars$name))
        #       ranks <- match(subE[,2],object@Vars$name)
        avgCon <- sapply(nodes,function(x)mean(ranks[subE[,1]==x]))
        return(order(avgCon))
      }
      
      # Endo:    
      EnM <- which(GroupVars$manifest & !GroupVars$exogenous)
      if (length(EnM) > 0)
      {
        GroupVars[EnM,] <- GroupVars[EnM,][ConOrd(GroupVars$name[EnM]),]
      }
      rm(EnM)
      
      ExM <- which(GroupVars$manifest & GroupVars$exogenous)
      # Exo:    
      if (length(ExM) > 0)
      {
        GroupVars[ExM,] <- GroupVars[ExM,][ConOrd(GroupVars$name[ExM]),]
      }
      rm(ExM)
      
      manNames <- GroupVars$name[GroupVars$manifest]
      Labels <- c(manNames,latNames)
    }
    
    Ni <- sum(GroupRAM$edge=="int")
    # Add intercept:
    if (any(object@RAM$edge=="int")) 
    {
      Labels[Labels=="1"] <- "_1"
      if (intStyle == "single") 
      {
        Labels <- c(Labels,"1")
      } else if (intStyle == "multi")
      {
        Labels <- c(Labels,rep("1",Ni))
      } 
    }
    
    nN <- length(Labels)
    
    # Extract edgelist:
    Edgelist <- GroupRAM[c("lhs","rhs")]
    Edgelist$lhs <- match(Edgelist$lhs,Labels)
    Edgelist$lhs[GroupRAM$edge=="int"] <- (nM+nL+1):nN
    Edgelist$rhs <- match(Edgelist$rhs,Labels)
    
    # Coerce to numeric matrix:
    Edgelist$lhs <- as.numeric(Edgelist$lhs)
    Edgelist$rhs <- as.numeric(Edgelist$rhs)
    Edgelist <- as.matrix(Edgelist)
    
    if (!allVars)
    {
      NodesInGroup <- sort(unique(c(Edgelist[,1],Edgelist[,2])))
      incl <- 1:nN %in% NodesInGroup
      
      Edgelist[,1] <- match(Edgelist[,1],NodesInGroup)
      Edgelist[,2] <- match(Edgelist[,2],NodesInGroup)
      
    } else incl <- 1:nN
    
    Labels <- Labels[incl]
    nN <- length(Labels)
    nM <- sum(manNames%in%Labels)
    nL <- sum(latNames%in%Labels)
    
    GroupVars <- GroupVars[GroupVars$name%in%Labels,]
    
    manInts <- Edgelist[GroupRAM$edge=="int" & GroupRAM$rhs%in%manNames,,drop=FALSE]
    latInts <- Edgelist[GroupRAM$edge=="int" & GroupRAM$rhs%in%latNames,,drop=FALSE]
    
    manIntsEndo <- manInts[!GroupVars$exogenous[manInts[,2]],,drop=FALSE]
    manIntsExo <- manInts[GroupVars$exogenous[manInts[,2]],,drop=FALSE]
    latIntsEndo <- latInts[!GroupVars$exogenous[latInts[,2]],,drop=FALSE]
    latIntsExo <- latInts[GroupVars$exogenous[latInts[,2]],,drop=FALSE]
    
    endoMan <- which(Labels%in%manNames&Labels%in%GroupVars$name[!GroupVars$exogenous])
    exoMan <- which(Labels%in%manNames&Labels%in%GroupVars$name[GroupVars$exogenous])
    endoLat <- which(Labels%in%latNames&Labels%in%GroupVars$name[!GroupVars$exogenous])
    exoLat <- which(Labels%in%latNames&Labels%in%GroupVars$name[GroupVars$exogenous])
    
    # Bidirectional:
    Bidir <- GroupRAM$edge == "<->"
    if (!grepl("mx",style,ignore.case=TRUE))
    {
      Bidir[GroupRAM$lhs==GroupRAM$rhs] <- FALSE
    }
    
    # Shape:
    Shape <- c(rep("square",nM),rep("circle",nL))
    if (any(GroupRAM$edge=="int")) Shape <- c(Shape,rep("triangle",Ni))
    
    Curve <- curve
    
    # Layout:
    if (layout=="tree" | layout=="circle" | layout=="circular")
    {
      #       if (all(!object@Vars$exogenous))
      #       {
      if (intStyle=="single")
      {
        # Curves:
        Curve <- ifelse(GroupRAM$lhs != GroupRAM$rhs & ((GroupRAM$lhs%in%manNames & GroupRAM$rhs%in%manNames) | (GroupRAM$lhs%in%latNames & GroupRAM$rhs%in%latNames)),curve,NA)
        Curve <- ifelse(GroupRAM$lhs%in%manNames,Curve,-1*Curve)
        Curve <- ifelse(GroupRAM$edge=="int" & GroupRAM$rhs%in%latNames,curve,-1*Curve)
        
        # Empty layout:
        Layout <- matrix(,length(Labels),2)
        
        # Add vertical levels:
        Layout[,2] <- ifelse(Labels%in%manNames,1,2)
        
        # Add vertical levels:
        Layout[Labels%in%manNames,1] <- seq(-1,1,length=nM)
        if (any(GroupRAM$edge=="int"))
        {
          sq <- seq(-1,1,length=nL+1)
          cent <- floor(median(1:(nL+1)))
          Layout[!Labels%in%manNames,1] <- sq[c(which(1:(nL+1) < cent),which(1:(nL+1) > cent),cent)]
        } else
        {
          Layout[Labels%in%latNames,1] <- seq(-1,1,length=nL)
        }
        
      } else if (intStyle=="multi")
      {          
        
        # Empty layout:
        Layout <- matrix(,length(Labels),2)
        
        # Add vertical levels:
        Layout[endoMan,2] <- 1
        Layout[endoLat,2] <- 2
        Layout[exoLat,2] <- 3
        Layout[exoMan,2] <- 4
        Layout[latIntsEndo[,1],2] <- 2
        Layout[latIntsExo[,1],2] <- 3
        
        if (residuals)
        {
          Layout[manIntsExo[,1],2] <- 4
          Layout[manIntsEndo[,1],2] <- 1
        } else {
          Layout[manIntsExo[,1],2] <- 5
          Layout[manIntsEndo[,1],2] <- 0            
        }
        
        # Add horizontal levels:
        if (nrow(manIntsEndo)>0)
        {
          Layout <- mixInts(endoMan,manIntsEndo,Layout,residuals=residuals)
        } else
        {
          if (length(endoMan)==1) Layout[endoMan,1] <- 0 else Layout[endoMan,1] <- seq(-1,1,length=length(endoMan))
        }
        if (nrow(manIntsExo)>0)
        {
          Layout <- mixInts(exoMan,manIntsExo,Layout,residuals=residuals)
        } else
        {
          if (length(exoMan)==1) Layout[exoMan,1] <- 0 else Layout[exoMan,1] <- seq(-1,1,length=length(exoMan))
        }
        
        if (nrow(latIntsEndo)>0)
        {
          Layout <- mixInts(endoLat,latIntsEndo,Layout,trim=TRUE)
        } else
        {
          Layout[endoLat,1] <- seq(-1,1,length=length(endoLat)+2)[-c(1,length(endoLat)+2)]
        }
        if (nrow(latIntsExo)>0)
        {
          Layout <- mixInts(exoLat,latIntsExo,Layout,trim=TRUE)
        } else
        {
          Layout[exoLat,1] <- seq(-1,1,length=length(exoLat)+2)[-c(1,length(exoLat)+2)]
        }
      } else stop("MeanStyle not supported")
      
      # Optimize layout if reorder = TRUE
      
    } else if (layout %in% c("tree2","circle2"))
    {
      # reingold-tilford layout
#       Layout <- layout.reingold.tilford
      
      # Set roots, in order of precedence:
        # exo man intercepts (incomming edges on exo man reversed)
        # Any exo man with outward edges (incomming edges on exo man reversed)
        # All exo man
        # Exo latent intercept
        # Exo latentes
        # Endo latent intercept
        # Endo latents
        # Manifest intercepts
        # Endo man with most outgoing edges
        # Endo man with least incoming edges
      # For all roots, base graph on graph with:
        # Double headed arrows removed
        # Arrow on other intercepts reversed
#       
#       if (any(GroupRAM$lhs %in% GroupVars$name[exoMan] & GroupRAM$edge %in% c("->","~>")))
#       {
#         roots <- sort(unique(Edgelist[,1][which(GroupRAM$lhs %in% GroupVars$name[exoMan] & GroupRAM$edge %in% c("->","~>"))]))
#         if (any(roots %in% manIntsExo[,2]))
#         {
#           roots <- manIntsExo[match(roots,manIntsExo[,2]),1]
#         } 
#       } else 
#         
      if (nrow(manIntsExo) > 0)
      {
        roots <- manIntsExo[,1]
      } else if (length(exoMan) > 0)
      {
        roots <- exoMan
      } else if (nrow(latIntsExo) > 0)
      {
        roots <- latIntsExo[,1]
      } else if (length(exoLat) > 0)
      {
        roots <- exoLat
      } else if (nrow(latIntsEndo) > 0)
      {
        roots <- latIntsEndo[,1]
      } else if (length(endoLat) > 0)
      {
        roots <- endoLat
      } else if (nrow(rbind(manIntsExo,manIntsEndo)) > 0)
      {
        roots <- rbind(manIntsExo,manIntsEndo)[,1]
      } else if (any(GroupRAM$edge %in% c("->","~>")))
      {
        roots <- Mode(Edgelist[,1][GroupRAM$edge %in% c("->","~>")])
      } else {
        roots <- Mode(c(Edgelist[,1],Edgelist[,2]))
      }
      
      Layout <- rtLayout(roots,GroupRAM,Edgelist,layout,exoMan)
      # Fix top level to use entire range:
      # Layout[Layout[,2]==max(Layout[,2]),1] <- seq(min(Layout[,1]),max(Layout[,1]),length.out=sum(Layout[,2]==max(Layout[,2])))
      # Center all horizontal levels:
      if (length(roots)>1) Layout[,1] <- ave(Layout[,1],Layout[,2],FUN = function(x) scale(x,TRUE,FALSE))
      
    } else Layout <- layout
    
    
    # loopRotation:
    if (layout%in%c("tree","tree2"))
    {
      loopRotation <- rep(0,nN)
      loopRotation[endoMan] <- pi
      loopRotation[exoMan] <- 0
      
      loopRotation[endoLat] <- 0
      noCons <- sapply(endoLat,function(x)nrow(Edgelist[(Edgelist[,1]==x|Edgelist[,2]==x) & (Edgelist[,1]%in%endoMan|Edgelist[,2]%in%endoMan),])==0)
      if (length(noCons)==0) noCons <- logical(0)
      loopRotation[endoLat][noCons] <- pi
      if (length(endoLat) > 1 & !(length(exoLat)==0&length(exoMan)==0))
      {
        if (length(exoLat) > 0 | any(endoLat %in% latIntsEndo[,2]))
        {
          loopRotation[endoLat[which.min(Layout[endoLat,1])]] <- ifelse(noCons[which.min(Layout[endoLat,1])],5/4*pi,7/4*pi)
          loopRotation[endoLat[which.max(Layout[endoLat,1])]] <- ifelse(noCons[which.min(Layout[endoLat,1])],3/4*pi,1/4*pi)          
        } else {
          loopRotation[endoLat[which.min(Layout[endoLat,1])]] <- 6/4 * pi
          loopRotation[endoLat[which.max(Layout[endoLat,1])]] <- 2/4 * pi
        }

      } else if (length(endoLat) == 1 && endoLat %in% latIntsEndo[,2])
      {
        loopRotation[endoLat] <- 6/4 * pi
      }
      
      loopRotation[exoLat] <- pi
      noCons <- sapply(exoLat,function(x)nrow(Edgelist[(Edgelist[,1]==x|Edgelist[,2]==x) & (Edgelist[,1]%in%exoMan|Edgelist[,2]%in%exoMan),])==0)
      if (length(noCons)==0) noCons <- logical(0)
      loopRotation[exoLat][noCons] <- 0
      if (length(exoLat) > 1 & length(exoMan)>0)
      {
        if (length(endoLat) > 0 | any(exoLat %in% latIntsExo[,2]))
        {
          loopRotation[exoLat[which.min(Layout[exoLat,1])]] <- ifelse(noCons[which.min(Layout[exoLat,1])],7/4*pi,5/4*pi)
          loopRotation[exoLat[which.max(Layout[exoLat,1])]] <- ifelse(noCons[which.min(Layout[exoLat,1])],1/4*pi,3/4*pi)
        } else {
          loopRotation[exoLat[which.min(Layout[exoLat,1])]] <- 6/4*pi
          loopRotation[exoLat[which.max(Layout[exoLat,1])]] <- 2/4*pi          
        }
      } else if (length(exoLat) == 1 && exoLat %in% latIntsExo[,2])
      {
        loopRotation[exoLat] <- 6/4 * pi
      }
      
      loopRotation <- loopRotation - 0.5 * (rotation-1) *pi
      
      if (any(GroupVars$exogenous) & optimizeLatRes)
      {
        ### For latents that have loops, find opposite of mean angle
        for (i in which(Labels%in%latNames & Labels%in%GroupRAM$lhs[GroupRAM$lhs==GroupRAM$rhs]))
        {
          # Layout subset of all connected:
          subEdgelist <- Edgelist[(Edgelist[,1]==i|Edgelist[,2]==i)&(Edgelist[,1]!=Edgelist[,2]),,drop=FALSE]              
          conNodes <- c(subEdgelist[subEdgelist[,1]==i,2],subEdgelist[subEdgelist[,2]==i,1])
          
          # Test for empty:
          if (nrow(subEdgelist)==0) conNodes <- sort(unique(c(Edgelist[,1:2])))
          
          subLayout <- Layout[conNodes,,drop=FALSE]
          Degrees <- apply(subLayout,1,function(x)atan2(x[1]-Layout[i,1],x[2]-Layout[i,2]))
          if (!grepl("lisrel",style,ignore.case=TRUE) | !any((Edgelist[,1]==i|Edgelist[,2]==i)&(Edgelist[,1]!=Edgelist[,2])&GroupRAM$edge=="<->"))
          {
            loopRotation[i] <- optimize(loopOptim,c(0,2*pi),Degrees=Degrees,maximum=TRUE)$maximum
          } else {
            # Layout subset of all connected:
            subEdgelist <- Edgelist[(Edgelist[,1]==i|Edgelist[,2]==i)&(Edgelist[,1]!=Edgelist[,2])&GroupRAM$edge=="<->",]
            conNodes <- c(subEdgelist[subEdgelist[,1]==i,2],subEdgelist[subEdgelist[,2]==i,1])
            
            # Test for empty:
            if (nrow(subEdgelist)==0) conNodes <- sort(unique(c(Edgelist[,1:2])))
            
            subLayout <- Layout[conNodes,]
            goodDegrees <- apply(subLayout,1,function(x)atan2(x[1]-Layout[i,1],x[2]-Layout[i,2]))
            loopRotation[i] <- optimize(loopOptim,c(min(goodDegrees-pi/4),max(goodDegrees+pi/4)),Degrees=Degrees,maximum=TRUE)$maximum
          }
        }
      }
    } else loopRotation <- NULL
    
    # Set curves and rotate:    
    if (layout %in% c("tree","tree2"))
    {
      
      inBetween <- function(x)
      {
        if (Layout[x[1],2]!=Layout[x[2],2]) return(0) else return(sum(Layout[Layout[,2]==Layout[x[1],2],1] > min(Layout[x,1]) & Layout[Layout[,2]==Layout[x[1],2],1] < max(Layout[x,1])))
      }
      # Curves:
      inBet <- apply(Edgelist,1,inBetween)
      inBet[inBet>0] <- as.numeric(as.factor(inBet[inBet>0]))
      if (!grepl("lisrel",style,ignore.case=TRUE) | all(!GroupVars$exogenous) | !residuals)
      {
        Curve <- ifelse(Layout[Edgelist[,1],2]==Layout[Edgelist[,2],2]&Edgelist[,1]!=Edgelist[,2]&GroupRAM$edge!="int",ifelse(inBet<1,0,curve+inBet/max(inBet)*curve),NA)
      } else {
        Curve <- ifelse(Layout[Edgelist[,1],2]==Layout[Edgelist[,2],2]&Edgelist[,1]!=Edgelist[,2]&GroupRAM$edge!="int" & Labels[Edgelist[,1]]%in%manNames & Labels[Edgelist[,2]]%in%manNames,ifelse(inBet<1,0,curve+inBet/max(inBet)*curve),NA)
      }
      ### If origin node is "right" of  destination node, flip curve:
      Curve[Layout[Edgelist[,1],1] > Layout[Edgelist[,2],1]] <- -1 * Curve[Layout[Edgelist[,1],1] > Layout[Edgelist[,2],1]]
      ## If endo man, flip again:
      Curve <- ifelse(Edgelist[,1]%in%endoMan | Edgelist[,2]%in%endoMan, -1 *  Curve, Curve)
      
      ### ORDINALIZE LAYOUT ###
      if (layout=="tree")
      {
        Layout[Layout[,2]>0&Layout[,2]<5,2] <- as.numeric(as.factor(Layout[Layout[,2]>0&Layout[,2]<5,2]))
        Layout[Layout[,2]==0,2] <- 0.5
        Layout[Layout[,2]==5,2] <- max(Layout[Layout[,2]<5,2]) + 0.5        
      }

      #         for (i in unique(Layout[,2]))
      #         {
      #           Layout[Layout[,2]==i,1] <- (as.numeric(as.factor(Layout[Layout[,2]==i,1])) - 1) / (sum(Layout[,2]==i) - 1)
      #         }
      # FLIP LAYOUT ###
      if (rotation==2) 
      {
        Layout <- Layout[,2:1]
        Layout[,1] <- -1 * Layout[,1]
      }
      if (rotation==3) 
      {
        Layout[,1] <- -1 * Layout[,1]
        Layout[,2] <- -1 * Layout[,2]
      }
      if (rotation==4) 
      {
        Layout <- Layout[,2:1]
        Layout[,2] <- -1 * Layout[,2]
      }
      Layout[,2] <- Layout[,2]-max(Layout[,2]) + 0.5
    }

    
    # Edge labels:
    if (edge.labels)
    {
      eLabels <- GroupRAM$label
    } else eLabels <- rep("",nrow(Edgelist))
    
    # vsize:
    vSize <- numeric(nN)
    vSize[Labels%in%manNames] <- sizeMan
    vSize[Labels%in%latNames] <- sizeLat
    vSize[Labels=="1"] <- sizeInt
    
    eColor <- rep(NA,nrow(Edgelist))
    #     tColor <- rep(rgb(0.5,0.5,0.5),nrow(GroupThresh))
    if (missing(threshold.color)) 
    {
      tColor <- rep("black",nrow(GroupThresh)) 
    } else {
      tColor <- rep(threshold.color,nrow(GroupThresh))
    }
    
    ### WHAT TO PLOT? ###
    if (grepl("path|diagram|mod",what,ignore.case=TRUE))
    {
      
    } else if (grepl("stand|std",what,ignore.case=TRUE))
    {
      Edgelist <- cbind(Edgelist,GroupRAM$std)
      if (edge.labels) eLabels <- as.character(round(GroupRAM$std,2))
    } else if (grepl("est|par",what,ignore.case=TRUE))
    {
      Edgelist <- cbind(Edgelist,GroupRAM$est)
      if (edge.labels) eLabels <- as.character(round(GroupRAM$est,2))
    } else if (grepl("eq|cons",what,ignore.case=TRUE))
    {
      #       eColor <- rep(rgb(0.5,0.5,0.5),nrow(Edgelist))
      unPar <- unique(object@RAM$par[object@RAM$par>0 & duplicated(object@RAM$par)])
      cols <- rainbow(max(c(object@RAM$par,GroupThresh$par)))
      for (i in unPar)
      {
        eColor[GroupRAM$par==i] <- cols[i]
      }
      if (nrow(GroupThresh) > 0)
      {
        for (i in 1:nrow(GroupThresh))
        {
          if (GroupThresh$par[i]>0 & sum(GroupThresh$par[i] == object@Thresholds$par) > 1 )
          {
            tColor[i] <- cols[GroupThresh$par[i]]
          }
        }
      }
    } else if (!grepl("col",what,ignore.case=TRUE)) stop("Could not detect use of 'what' argument")
    
    ### VERTEX COLOR ###
    # Set default if list:
    if (is.list(color))
    {
      colList <- color
      color <- rep("",nN)
      if (!is.null(colList$man))
      {
        color[Labels%in%manNames] <- colList$man
      }
      if (!is.null(colList$lat))
      {
        color[Labels%in%latNames] <- colList$lat
      }
      if (!is.null(colList$int))
      {
        color[Labels=="1"] <- colList$int
      }
    }
    
    if (!missing(groups))
    {
      NodeGroups <- groups
      
      Ng <- length(NodeGroups)
      
      if (length(color)==1)
      {
        Vcolors <- rep(color,nN)
      } else if (length(color)==nM)
      {
        Vcolors <- c(color,rep("",nN-nM))
      } else if (length(color)==nN)
      {
        Vcolors <- color  
      } else if (length(color)!=Ng)
      {
        stop("'color' vector not of appropriate length")
      }
      
      if (missing(manifests) & any(sapply(NodeGroups,mode)!="character")) warning("Groups specified numerically and 'manifests' not supplied. Results might be unexpected.")
      
      if (length(color)==Ng)
      {
        Vcolors <- rep("",nN)
        for (g in 1:Ng)
        {
          if (mode(NodeGroups[[g]])=="character") NodeGroups[[g]] <- match(NodeGroups[[g]],Labels)
          Vcolors[NodeGroups[[g]]] <- color[g]
        }
        
        #         Vcolors[Vcolors=="" & Labels%in%manNames] <- "white"
      }
      
    } else 
    {
      NodeGroups <- NULL
      
      if (length(color)==1)
      {
        Vcolors <- rep(color,nN)
      } else if (length(color)==nM)
      {
        Vcolors <- c(color,rep(NA,nN-nM))
      } else if (length(color)==nN)
      {
        Vcolors <- color  
      } else stop("'color' vector not of appropriate length")
    }
    
    # If missing color, obtain weighted mix of connected colors:
    if (any(Vcolors==""))
    {
      if (mixCols)
      {
        VcolorsBU <- Vcolors
        W <- 1
  #       for (i in 1:(nM+nL))
        for (i in 1:nN)
        {
          if (Vcolors[i]=="")
          {
            cons <- c(Edgelist[Edgelist[,1]==i,2],Edgelist[Edgelist[,2]==i,1])
            if (ncol(Edgelist) == 3)
            {
              W <- abs(c(Edgelist[Edgelist[,1]==i,3],Edgelist[Edgelist[,2]==i,3]))
              W <- W[VcolorsBU[cons]!=""]
            }
            cons <- cons[VcolorsBU[cons]!=""]
            if (length(cons)>0)
            {
              Vcolors[i] <- mixColfun(VcolorsBU[cons],W)
            } else Vcolors[i] <- NA
          }
        }
      }
      Vcolors[Vcolors==""] <- NA
    }
    
    
    if (grepl("col",what,ignore.case=TRUE))
    {
      #       eColor <- character(nrow(Edgelist))
      for (i in 1:nrow(Edgelist))
      {
        cols <- Vcolors[Edgelist[i,]]
        if (!all(cols=="background")) eColor[i] <- mixColfun(cols[cols!="background"])
      }
    }
    
    if (!missing(whatLabels))
    {
      if (grepl("path|diagram|model|name|label",whatLabels,ignore.case=TRUE))
      {
        eLabels <- GroupRAM$label
      } else if (grepl("stand|std",whatLabels,ignore.case=TRUE))
      {
        eLabels <- as.character(round(GroupRAM$std,2))
      } else if (grepl("est|par",whatLabels,ignore.case=TRUE))
      {
        eLabels <- as.character(round(GroupRAM$est,2))
      } else if (grepl("eq|cons",whatLabels,ignore.case=TRUE))
      {
        eLabels <- GroupRAM$par
      } else if (grepl("no|omit|hide|invisible",whatLabels,ignore.case=TRUE))
      {
        eLabels <- rep("",nrow(Edgelist))
      } else stop("Could not detect use of 'whatLabels' argument")
    }
    
    # Abbreviate:
    if (nCharEdges>0 & !"edges"%in%as.expression)
    {
      eLabels <- abbreviate(eLabels,nCharEdges)
    }
    if (nCharNodes>0 & !"nodes"%in%as.expression)
    {
      Labels <- abbreviate(Labels,nCharNodes)
    }
    
    #     ### CONVERT TO LISREL STYLE ###
    if (grepl("lisrel",style,ignore.case=TRUE) & residuals)
    {
      isResid <- GroupRAM$edge == "<->" & GroupRAM$lhs != GroupRAM$rhs
    } else isResid <- rep(FALSE,nrow(Edgelist))
    
    
    #       nResid <- length(whichResid)
    #       Edgelist[whichResid,1] <- (nN+1):(nN+nResid)
    #       rots <- loopRotation[Edgelist[whichResid,2]]
    #       Lresid <- matrix(,nResid,2)
    #       hLength <- diff(range(Layout[,1]))
    #       vLength <- diff(range(Layout[,2]))
    #       for (i in 1:nResid)
    #       {
    #         Lresid[i,1] <- Layout[Edgelist[whichResid[i],2],1] + sin(rots[i]) * residScale * 0.25 * hLength/vLength
    #         Lresid[i,2] <- Layout[Edgelist[whichResid[i],2],2] + cos(rots[i]) * residScale * 0.25
    #       }
    #       
    #       # Add nodes:
    #       Layout <- rbind(Layout,Lresid)
    #       Labels <- c(Labels,rep("",nResid))
    #       Shape <- c(Shape,rep("circle",nResid))
    #       loopRotation <- NULL
    #       vSize <- c(vSize,rep(0,nResid))
    #       Vcolors <- c(Vcolors,rep(rgb(0,0,0,0),nResid))
    #     }
    
    if (grepl("mx",style,ignore.case=TRUE)) LoopAsResid <- FALSE else LoopAsResid <- TRUE
    
    ### ROTATE IF CIRCLE:
    if (layout%in%c("circle","circle2"))
    {
      if (rotation%in%c(2,4)) stop("Circle layout only supported if rotation is 1 or 3")
      underMean <- Layout[,2] < mean(Layout[,2])
      Layout[,2] <- -1*Layout[,2] + max(Layout[,2]) + 0.5
      Ltemp <- Layout
      unVert <- sort(unique(Layout[,2]))
      for (i in unVert)
      {
        l <- sum(Layout[,2]==i)
        sq <- seq(0,2*pi,length=l+1)[-(l+1)] + pi/l
        c <- 1
        for (j in order(Layout[Layout[,2]==i,1]))
        {
          Ltemp[which(Layout[,2]==i)[j],] <- c(RotMat(sq[c])%*%c(0,i))
          c <- c + 1
        }
      }
      Layout <- Ltemp
      
      # loopRotation:
      loopRotation <- apply(Layout,1,function(x)atan2(x[1],x[2]))
      loopRotation <- ifelse(underMean,loopRotation,(loopRotation+pi)%%(2*pi))
    }
    
    if (layout == "spring") loopRotation <- NULL
    
    
    if (!allVars)
    {
      NodeGroups2 <- NodeGroups
      if (!is.null(NodeGroups2))
      {
        newNodes <- match(1:length(AllLabs),(1:length(AllLabs))[incl])
        for (g in 1:length(NodeGroups2))
        {
          NodeGroups2[[g]] <- newNodes[NodeGroups2[[g]]]
          NodeGroups2[[g]] <- NodeGroups2[[g]][!is.na(NodeGroups2[[g]])]
        }
      }
    }
    
    ### Compute margins ###
    if (missing(mar))
    {
      Mar <- c(5,5,5,5)
      
      # Add 3 to top and bottom for residuals if lisrel style is used:
      if (grepl("lisrel",style,ignore.case=TRUE)) Mar[c(1,3)] <- Mar[c(1,3)] + 3
      
#       # Add 4 to bottom if there are endo man residual correlations:
#       if (any(Edgelist[,1]%in%endoMan & Edgelist[,2]%in%endoMan & Edgelist[,1]!=Edgelist[,2]))
#       {
#         Mar[1] <- Mar[1] + 4
#       }
#       
#       # Add 4 to top if there are endo man residual correlations:
#       if (any(Edgelist[,1]%in%exoMan & Edgelist[,2]%in%exoMan & Edgelist[,1]!=Edgelist[,2]))
#       {
#         Mar[3] <- Mar[3] + 4
#       }
#       
      # Add 3 to top if top consist of latent residuals:
      if (length(exoMan)==0)
      {
        Mar[3] <- Mar[3] + 3
      }
      
      # Rotate:
      Mar <- Mar[(0:3 + (rotation-1)) %% 4 + 1]
      
      # Add 2 to top for title:
      if (title) Mar[3] <- Mar[3] + 2
    } else Mar <- mar
    
    # Overwrite edge colors if appropriate:
    if (!missing(edge.color))
    {
      eColor <- edge.color
      if (thresholds & missing(threshold.color))
      {
        tColor <- rep(edge.color,length(tColor))
      }
    }
    
    # Fixed and free edges:
    if (length(freeStyle) > 2 | length(fixedStyle) > 2) warning("'freeStyle' and 'fixedStyle' are assumed to be vectors of at most length 2. Unexpected results will probably occur.")
    # lty:
    lty <- rep(1,nrow(GroupRAM))
    
    # fixedStyle
    if (any(is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle))) lty <- ifelse(GroupRAM$fixed,as.numeric(fixedStyle[is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle)]),lty) 
    
    if (any(qgraph:::isColor(fixedStyle) & !(is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle)))) eColor[GroupRAM$fixed] <- fixedStyle[qgraph:::isColor(fixedStyle) & !(is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle))]
    
    
    # freeStyle:
    if (any(is.numeric(freeStyle) | grepl("\\d+",freeStyle))) lty <- ifelse(GroupRAM$fixed,lty,as.numeric(freeStyle[is.numeric(freeStyle) | grepl("\\d+",freeStyle)]))
    
    if (any(qgraph:::isColor(freeStyle) & !(is.numeric(freeStyle) | grepl("\\d+",freeStyle)))) eColor[!GroupRAM$fixed] <- freeStyle[qgraph:::isColor(freeStyle) & !(is.numeric(freeStyle) | grepl("\\d+",freeStyle))]
    
    # Directed settings:
    
    Directed <- GroupRAM$edge!="--"
    
    # Convert labels to expressions:
    if ("edges"%in%as.expression)
    {
      eLabels <- as.expression(parse(text=eLabels))
    }    # Convert labels to expressions:
    if ("nodes"%in%as.expression)
    {
      Labels <- as.expression(parse(text=Labels))
    }
    # Restore layout function:
    if (!is.null(layoutFun)) Layout <- layoutFun
    
    qgraphRes[[which(Groups==gr)]] <- qgraph(Edgelist,
                                             labels=Labels,
                                             bidirectional=Bidir,
                                             directed=Directed,
                                             shape=Shape,
                                             layout=Layout,
                                             lty=lty,
                                             loopRotation=loopRotation,
                                             curve=Curve,
                                             edge.labels=eLabels,
                                             mar=Mar,
                                             vsize = vSize,
                                             edge.color=eColor,
                                             groups=NodeGroups2,
                                             color=Vcolors,
                                             residuals=LoopAsResid,
                                             residScale = residScale,
                                             residEdge = isResid,
                                             edgelist = TRUE,
                                             curveDefault = curveDefault,
                                             knots = GroupRAM$knot,
                                             curvePivot = curvePivot,
                                             ...)
    
    if (thresholds)
    {
      # Overwrite color to white if bg is dark (temporary solution)
      if (missing(threshold.color) & missing(edge.color)) 
      {
        if (mean(col2rgb(qgraphRes[[which(Groups==gr)]]$background)/255) <= 0.5) tColor <- rep("white",length(tColor))
      }
      if (nrow(GroupThresh) > 0)
      {
        for (i in 1:nrow(GroupThresh))
        {
          node <- which(Labels==GroupThresh$lhs[i])
          # Compute side:
          IntSide <- 1
          if (layout=="tree")
          {
            if (rotation%in%c(1,3))
            {
              IntSide <- ifelse(Layout[node,2]>mean(Layout[,2]),3,1)
            } else {
              IntSide <- ifelse(Layout[node,1]>mean(Layout[,1]),4,2)
            }
          } else {
            IntSide <- sum((atan2(qgraphRes[[which(Groups==gr)]]$layout[node,1],qgraphRes[[which(Groups==gr)]]$layout[node,2])+pi)%%(2*pi) > c(0,pi/2,pi,1.5*pi))
          }
          IntInNode(qgraphRes[[which(Groups==gr)]]$layout[node,,drop=FALSE],vSize[node],Shape[node],pnorm(GroupThresh$est[i]),width=0.5,triangles=FALSE,col=tColor[i],IntSide,!ThreshAtSide)
        }
      }
    }
    
    if (title)
    {
      #       if (length(Groups)==1) title("Path Diagram",line=3) else title(paste0("Path Diagram for group '",gr,"'"),line=3)
      title(gr,line=3,col.main=title.color)
    }
  }
  par(ask=askOrig)
  if (length(qgraphRes)==1) qgraphRes <- qgraphRes[[1]]
  invisible(qgraphRes)
}
