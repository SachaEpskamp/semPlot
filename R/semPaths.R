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
# Boker

# manifests: vector of manifest labels ordered
# latents: vector of latents ordered
# fixedStyle: if coercible to numeric lty is assigned this value, else a color for color representation. If this argument is not a number or color representation the edge is not displayed differently.


semPaths <- function(object,what="paths",whatLabels,style,layout="tree",intercepts=TRUE,residuals=TRUE,thresholds=TRUE,
                     intStyle="multi",rotation=1,curve, curvature = 1, nCharNodes=3,nCharEdges=3,sizeMan = 5,sizeLat = 8,
                     sizeInt = 2,  sizeMan2 ,sizeLat2 ,sizeInt2, shapeMan, shapeLat, shapeInt = "triangle", ask,mar,title,title.color="black",
                     title.adj = 0.1, title.line = -1, title.cex = 0.8,
                     include,combineGroups=FALSE,manifests,latents,groups,color, residScale,gui=FALSE,allVars=FALSE,edge.color,
                     reorder=TRUE,structural=FALSE,ThreshAtSide=FALSE,thresholdColor,thresholdSize = 0.5, fixedStyle=2,freeStyle=1,
                     as.expression=character(0),optimizeLatRes=FALSE,inheritColor=TRUE,levels,nodeLabels,edgeLabels,
                     pastel=FALSE,rainbowStart=0,intAtSide,springLevels=FALSE,nDigits=2,exoVar,exoCov=TRUE,centerLevels=TRUE,
                     panelGroups=FALSE,layoutSplit = FALSE, measurementLayout = "tree", subScale, subScale2, subRes = 4, 
                     subLinks, modelOpts = list(), curveAdjacent = "<->", edge.label.cex = 0.6, cardinal =  "none", 
                     equalizeManifests = FALSE,  covAtResiduals = TRUE, bifactor, optimPoints = 1:8 * (pi/4),
                     ...){
  
#   c("exo cov","load dest","endo man cov")
  
  # Check if input is combination of models:
  call <- paste(deparse(substitute(object)), collapse = "")
  if (grepl("\\+",call)) 
  {
    args <- unlist(strsplit(call,split="\\+"))
    obs <- lapply(args,function(x)semPlotModel(eval(parse(text=x))))
    object <- obs[[1]]
    for (i in 2:length(obs)) object <- object + obs[[i]]
  }
  
  if (!"semPlotModel"%in%class(object)) object <- do.call(semPlotModel,c(list(object),modelOpts))
  stopifnot("semPlotModel"%in%class(object))
  
  # if (gui) return(do.call(semPathsGUI,as.list(match.call())[-1]))
  
  ### edgeConnectPoints dummy:
  ECP <- NULL
  
  # Set defaults size and shape:
  if (missing(sizeMan2)) sizeMan2 <- sizeMan
  if (missing(sizeLat2)) sizeLat2 <- sizeLat
  if (missing(sizeInt2)) sizeInt2 <- sizeInt
  
  if (missing(shapeMan))
  {
    if (sizeMan == sizeMan2) shapeMan <- "square" else shapeMan <- "rectangle"
  }
  
  if (missing(shapeLat))
  {
    if (sizeLat == sizeLat2) shapeLat <- "circle" else shapeLat <- "ellipse"
  }
  
  # Check:
  if (missing(intAtSide)) intAtSide <- !residuals
  
  if (!rotation%in%1:4)
  {
    stop("Rotation must be 1, 2 3 or 4.")
  }
  if (any(object@Pars$edge=="int")) 
  {
    object@Vars$name[object@Vars$name=="1"] <- "_1"
    object@Pars$lhs[object@Pars$lhs=="1"] <- "_1"
    object@Pars$rhs[object@Pars$rhs=="1"] <- "_1"
  }
  
  # Check if layout is not character of length 1. If so, set layout <- "spring" as dummy:
  if (!is.character(layout) || length(layout) > 1)
  {
    layoutFun <- layout
    layout <- "spring"
  } else layoutFun <- NULL
  
  if (!missing(bifactor) & !layout %in% c("tree2","tree3","circle2","circle3"))
  {
    warning("'bifactor' argument only supported in layouts 'tree2', 'tree3', 'circle2' and 'circle3'")
  }
  
  
  if (missing(curve))
  {
    if (layout %in% c("tree","tree2","tree3"))
    {
      curve <- 1
    } else {
      curve <- 0
    }
  }
  curveDefault <- curve
  
#   if (missing(curvePivot))
#   {
#     curvePivot <- grepl("tree",layout)
#   }
  
  if (missing(whatLabels))
  {
    edge.labels <- TRUE
  } else
  {
    edge.labels <- FALSE    
  }
  
  if (missing(as.expression)) 
  {
    if ("lisrel" %in% unlist(sapply(object@Original,class)))
    {
      as.expression <- "edges"
    } else as.expression <- ""
  }
  
  # Set and check style: 
  if (missing(style))
  {
    if ("lisrel" %in% unlist(sapply(object@Original,class)))
    {
      style <- "lisrel"
    } else style <- "OpenMx"
  }
  if (grepl("ram",style,ignore.case=TRUE)) style <- "OpenMx"
  if (!grepl("mx|lisrel",style,ignore.case=TRUE)) stop("Only OpenMx (ram) or LISREL style is currently supported.")
  #   if (grepl("mx",style,ignore.case=TRUE) & !missing(residScale)) warning("'residScale' ingored in OpenMx style")
  if (missing(residScale)) residScale <- sizeMan
  
  # Set exoVar default:
  if (missing(exoVar)) exoVar <- !grepl("lis",style,ignore.case=TRUE)
  
  #   residScale <- residScale * 1.75
  # Remove means if means==FALSE
  if (intercepts==FALSE)
  {
    object@Pars <- object@Pars[object@Pars$edge!="int",]
  }
  # Set true exogenous:
  object@Vars$trueExo <- !object@Vars$name %in% object@Pars$rhs[object@Pars$edge %in% c("->","~>")]
  # Remove true exo variances:
  if (!exoVar)
  {
    object@Pars <- object@Pars[!((object@Pars$lhs %in% object@Vars$name[object@Vars$trueExo] | object@Pars$rhs %in% object@Vars$name[object@Vars$trueExo]) & object@Pars$edge == "<->" & (object@Pars$rhs == object@Pars$lhs)), ]
  }
  if (!exoCov)
  {
    object@Pars <- object@Pars[!((object@Pars$lhs %in% object@Vars$name[object@Vars$trueExo] | object@Pars$rhs %in% object@Vars$name[object@Vars$trueExo]) & object@Pars$edge == "<->" & (object@Pars$rhs != object@Pars$lhs)), ]
  }
  
  # Remove residuals if residuals=FALSE
  if (residuals==FALSE)
  {
    object@Pars <- object@Pars[!(object@Pars$edge=="<->"&object@Pars$lhs==object@Pars$rhs),]
  }  
  
  # Combine groups if combineGroups=TRUE:
  if (combineGroups)
  {
    object@Pars$group <- ""
  }
  # Within - Between framework:
  if (is.null(object@Pars$BetweenWithin))
  {
    object@Pars$BetweenWithin <- ''
    if (nrow(object@Thresholds) > 0)
    {
      object@Thresholds$BetweenWithin <- ''
    }
  }
  
  if ((length(unique(object@Pars$BetweenWithin)) > 1 && !all(unique(object@Pars$BetweenWithin) %in% c('Within','Between'))) | length(unique(object@Pars$BetweenWithin)) > 2) stop("BetweenWithin must be labeled 'Between' and 'Within' only")
  
  if (length(unique(object@Pars$BetweenWithin)) == 2)
  {
    object@Pars$group <- paste(object@Pars$group,'-',object@Pars$BetweenWithin)
    object@Pars$group <- gsub('\\s+\\-\\s+(?=Within$)','',object@Pars$group,perl=TRUE)
    object@Pars$group <- gsub('\\s+\\-\\s+(?=Between$)','',object@Pars$group,perl=TRUE)
    
    if (nrow(object@Thresholds) > 0)
    {
      object@Thresholds$group <- paste(object@Thresholds$group,'-',object@Thresholds$BetweenWithin)
      object@Thresholds$group <- gsub('\\s+\\-\\s+(?=Within$)','',object@Thresholds$BetweenWithin,perl=TRUE)
      object@Thresholds$group <- gsub('\\s+\\-\\s+(?=Between$)','',object@Thresholds$BetweenWithin,perl=TRUE)
    }
    
  }
  # Set title:
  if (missing(title))
  {
    # Check titles:
    title <- length(unique(object@Pars$group))>1
  }
  # If structural, remove all manifest from Pars:
  if (structural)
  {
    object@Pars <- object@Pars[!(object@Pars$lhs %in% object@Vars$name[object@Vars$manifest] | object@Pars$rhs %in% object@Vars$name[object@Vars$manifest]),]
    object@Vars <- object@Vars[!object@Vars$manifest,]
    object@Thresholds <- data.frame()
  }  
  
  # Add rows for bidirectional edges:
  if (any(object@Pars$edge=="<->" & object@Pars$lhs != object@Pars$rhs))
  {
    bidirs <- object@Pars[object@Pars$edge=="<->" & object@Pars$lhs != object@Pars$rhs,]
    bidirs[c("lhs","rhs")] <- bidirs[c("rhs","lhs")]
    bidirs$par <- -1
    object@Pars <- rbind(object@Pars,bidirs)
  }
  object@Pars <- object@Pars[!duplicated(object@Pars),]
  
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
      if (pastel)
      {
        if (length(groups) == 1) color <- "white" else 
        color <- rainbow_hcl(length(groups), start = rainbowStart * 360, end = (360 * rainbowStart + 360*(length(groups)-1)/length(groups)))
      } else {
        if (length(groups) == 1) color <- "white" else 
        color <- rainbow(length(groups), start = rainbowStart, end = (rainbowStart + (max(1,length(groups)-1))/length(groups)) %% 1)   
      }
    }
  } else {
    if (missing(color)) 
    {
      color <- "background"
    } 
  }

  
  #   # Define exogenous variables (only if any is NA):
  #   if (any(is.na(object@Vars$exogenous)))
  #   {
  #     if (any(!is.na(object@Vars$exogenous)))
  #     {
  #       exoOrig <- object@Vars$exogenous
  #       repExo <- TRUE
  #     } else repExo <- FALSE
  #     object@Vars$exogenous <- FALSE
  #     for (i in which(!object@Vars$manifest))
  #     {
  #       if (!any(object@Pars$edge[object@Pars$rhs==object@Vars$name[i]] %in% c("~>","->") & object@Pars$lhs[object@Pars$rhs==object@Vars$name[i]]%in%latNames))
  #       {
  #         object@Vars$exogenous[i] <- TRUE
  #       }
  #     }
  #     for (i in which(object@Vars$manifest))
  #     {
  #       if (all(object@Pars$lhs[object@Pars$rhs==object@Vars$name[i] & object@Pars$lhs%in%latNames]%in%object@Vars$name[object@Vars$exogenous]) &
  #             all(object@Pars$rhs[object@Pars$lhs==object@Vars$name[i] & object@Pars$rhs%in%latNames]%in%object@Vars$name[object@Vars$exogenous]) &
  #             !any(object@Pars$rhs==object@Vars$name[i] & object@Pars$edge=="~>"))
  #       {
  #         object@Vars$exogenous[i] <- TRUE
  #       }
  #     }
  #     
  #     # If all exo, treat all as endo:
  #     if (all(object@Vars$exogenous) | layout%in%c("circle","circle2","circle3"))
  #     {
  #       object@Vars$exogenous <- FALSE
  #     }
  #     # If al endo, treat formative manifest as exo (MIMIC mode), unless all manifest are formative.
  #     if (!any(object@Vars$exogenous))
  #     {
  #       if (any(object@Vars$manifest & (object@Vars$name%in%object@Pars$rhs[object@Pars$edge %in% c("~>","--","->")])))
  #         object@Vars$exogenous[object@Vars$manifest & !(object@Vars$name%in%object@Pars$rhs[object@Pars$edge %in% c("~>","--","->")])] <- TRUE
  #     }
  #     if (repExo)
  #     {
  #       object@Vars$exogenous[!is.na(exoOrig)] <- exoOrig[!is.na(exoOrig)]
  #     }
  #   }
  object <- defExo(object, layout)
  
  Groups <- unique(object@Pars$group)
  qgraphRes <- list()
  if (missing(ask))
  {
    if (length(Groups)>1) ask <- TRUE else ask <- FALSE
  }
  askOrig <- par("ask")
  
  if (missing(include)) include <- 1:length(Groups)
  
  if (panelGroups)
  {
    layout(t(1:length(include)))
  }
  
  # Reassign labels (temporary solution for excluding vars in multi group)
  AllLabs <- Labels
  AllMan <- manNames
  AllLat <- latNames
  
  par(ask=ask)

  ### If no sub, set sub to 0 (root sub)
  if (is.null(object@Pars$sub)) 
  {    
    if (!layoutSplit)
    {
      object@Pars$sub <- 0
    } else {
      object@Pars$sub <- 1
      
      ### Detect manifest children of each latent:
      for (i in seq_along(latNames))
      {
        # Connected manifests:
        object@Pars$sub[object@Pars$lhs == latNames[i] & object@Pars$rhs%in%manNames & object@Pars$edge%in%c('->','~>')] <- i+1
        
        # Intercepts and variances connected to these manifests:
        object@Pars$sub[object@Pars$rhs%in%object@Pars$rhs[object@Pars$rhs%in%manNames & object@Pars$sub==(i+1)] & object@Pars$edge == "int"] <- i+1
        
        object@Pars$sub[object@Pars$rhs%in%object@Pars$rhs[object@Pars$rhs%in%manNames & object@Pars$sub==(i+1)] & object@Pars$lhs%in%object@Pars$rhs[object@Pars$rhs%in%manNames & object@Pars$sub==(i+1)] & object@Pars$edge == "<->"] <- i+1
      }
      
      # Remove manifests already in a submodel from model 1:
      ManInSub <- manNames[ manNames %in% object@Pars$lhs[object@Pars$sub > 1] | manNames %in% object@Pars$rhs[object@Pars$sub > 1]]
      object@Pars$sub[ (object@Pars$lhs %in% ManInSub | object@Pars$rhs %in% ManInSub) & object@Pars$sub == 1] <- 0
    }
  }
  if (missing( subLinks))
  {
    subLinks <- latNames
  }
  if (missing(subScale))
  {
    subScale <- 0.1 + 0.3 * (1 / max(1,max(object@Pars$sub)))
  }
  if (missing(subScale2))
  {
    subScale2 <- subScale * 1.5
  }
  
  layoutMain <- layout
  rotationMain <- rotation


  for (gr in Groups[(1:length(Groups))%in%include])
  {  
    grSub <- object@Pars$sub[object@Pars$group==gr]
    if (length(unique(grSub)) == 1) grSub[] <- 0
    
    # List to store results (layout and curve of submodel 1 - n, 1 being root)
    subModList <- list()
    
    # Start sub loop (one loop per sub and once for whole graph, if number of subs is 1 ignore)
    # -1 is rerun for main graph:
    for (Sub in (max(grSub):0)[max(grSub):0%in%c(grSub,0)])
    {
      
      if (Sub > 0)
      {
        GroupPars <- object@Pars[object@Pars$group==gr & object@Pars$sub==Sub,]
        GroupVars <- object@Vars
        GroupThresh <- object@Thresholds[object@Thresholds$group==gr & object@Pars$sub==Sub,]        
        
      } else 
      {
        GroupPars <- object@Pars[object@Pars$group==gr,]
        GroupVars <- object@Vars
        GroupThresh <- object@Thresholds[object@Thresholds$group==gr,]
      }
      
      if (Sub > 1)
      {
        GroupVars$exogenous <- FALSE
        rotation <- 1
        layout <- measurementLayout
      } else {
        rotation <- rotationMain
        layout <- layoutMain
      }

      
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
          E <- GroupPars[c("lhs","rhs")]
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
      
      Ni <- sum(GroupPars$edge=="int")
      # Add intercept:
      if (any(object@Pars$edge=="int")) 
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
      Edgelist <- GroupPars[c("lhs","rhs")]
      Edgelist$lhs <- match(Edgelist$lhs,Labels)
      Edgelist$lhs[GroupPars$edge=="int"] <- (nM+nL+1):nN
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
      
      manInts <- Edgelist[GroupPars$edge=="int" & GroupPars$rhs%in%manNames,,drop=FALSE]
      latInts <- Edgelist[GroupPars$edge=="int" & GroupPars$rhs%in%latNames,,drop=FALSE]
      
      manIntsEndo <- manInts[!GroupVars$exogenous[manInts[,2]],,drop=FALSE]
      manIntsExo <- manInts[GroupVars$exogenous[manInts[,2]],,drop=FALSE]
      latIntsEndo <- latInts[!GroupVars$exogenous[latInts[,2]],,drop=FALSE]
      latIntsExo <- latInts[GroupVars$exogenous[latInts[,2]],,drop=FALSE]
      
      endoMan <- which(Labels%in%manNames&Labels%in%GroupVars$name[!GroupVars$exogenous])
      exoMan <- which(Labels%in%manNames&Labels%in%GroupVars$name[GroupVars$exogenous])
      endoLat <- which(Labels%in%latNames&Labels%in%GroupVars$name[!GroupVars$exogenous])
      exoLat <- which(Labels%in%latNames&Labels%in%GroupVars$name[GroupVars$exogenous])
      
      # Bidirectional:
      Bidir <- GroupPars$edge == "<->"
      if (!grepl("mx",style,ignore.case=TRUE))
      {
        Bidir[GroupPars$lhs==GroupPars$rhs] <- FALSE
      }
      
      # Shape:
      Shape <- c(rep(shapeMan,nM),rep(shapeLat,nL))
      if (any(GroupPars$edge=="int")) Shape <- c(Shape,rep(shapeInt,Ni))
      
      Curve <- curve
      
      # Layout:
      if (layout=="tree" | layout=="circle" | layout=="circular")
      {
        #       if (all(!object@Vars$exogenous))
        #       {
        if (intStyle=="single")
        {
          # Curves:
          Curve <- ifelse(GroupPars$lhs != GroupPars$rhs & ((GroupPars$lhs%in%manNames & GroupPars$rhs%in%manNames) | (GroupPars$lhs%in%latNames & GroupPars$rhs%in%latNames)),curve,NA)
          Curve <- ifelse(GroupPars$lhs%in%manNames,Curve,-1*Curve)
          Curve <- ifelse(GroupPars$edge=="int" & GroupPars$rhs%in%latNames,curve,-1*Curve)
          
          # Empty layout:
          Layout <- matrix(,length(Labels),2)
          
          # Add vertical levels:
          Layout[,2] <- ifelse(Labels%in%manNames,1,2)
          
          # Add vertical levels:
          Layout[Labels%in%manNames,1] <- seq(-1,1,length=nM)
          if (any(GroupPars$edge=="int"))
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
          
          if (intAtSide)
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
            Layout <- mixInts(endoMan,manIntsEndo,Layout,intAtSide=intAtSide)
          } else
          {
            if (length(endoMan)==1) Layout[endoMan,1] <- 0 else Layout[endoMan,1] <- seq(-1,1,length=length(endoMan))
          }
          if (nrow(manIntsExo)>0)
          {
            Layout <- mixInts(exoMan,manIntsExo,Layout,intAtSide=intAtSide)
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

          
          if (equalizeManifests)
          {
            # Max of number of nodes in lvls 1 and 4:
            EndoHorRange <- max(sapply(c(0,1,4,5), function(x) sum(Layout[,2] == x)))
            for (lvl in c(0,1,4,5)) Layout[Layout[,2]==lvl,1] <- sum(Layout[,2]==lvl) * Layout[Layout[,2]==lvl,1] / EndoHorRange
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
        #       if (any(GroupPars$lhs %in% GroupVars$name[exoMan] & GroupPars$edge %in% c("->","~>")))
        #       {
        #         roots <- sort(unique(Edgelist[,1][which(GroupPars$lhs %in% GroupVars$name[exoMan] & GroupPars$edge %in% c("->","~>"))]))
        #         if (any(roots %in% manIntsExo[,2]))
        #         {
        #           roots <- manIntsExo[match(roots,manIntsExo[,2]),1]
        #         } 
        #       } else 
        #         
        
        # If bifactor is assigned and exists, bifactor becomes root and all edges not connected to bifactor are reversed:    
        if (!missing(bifactor) && any(bifactor %in% Labels))
        {
          roots <- which(Labels %in% bifactor)
        } else {
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
          } else if (any(GroupPars$edge %in% c("->","~>")))
          {
            roots <- Mode(Edgelist[,1][GroupPars$edge %in% c("->","~>")])
          } else {
            roots <- Mode(c(Edgelist[,1],Edgelist[,2]))
          }
        }
        
        Layout <- rtLayout(roots,GroupPars,Edgelist,layout,exoMan)
        # Fix top level to use entire range:
        # Layout[Layout[,2]==max(Layout[,2]),1] <- seq(min(Layout[,1]),max(Layout[,1]),length.out=sum(Layout[,2]==max(Layout[,2])))
        # Center all horizontal levels:
#         if (centerLevels) if (length(roots)>1) Layout[,1] <- ave(Layout[,1],Layout[,2],FUN = function(x) scale(x,TRUE,FALSE))
        if (centerLevels) Layout[,1] <- ave(Layout[,1],Layout[,2],FUN = function(x) scale(x,TRUE,FALSE))
        
      } else if (layout%in%c("tree3","circle3"))
      {
        
        # Igraph:
        # Select only directed edges:
        Edgelist2 <- Edgelist[GroupPars$edge%in%c("->","~>"),]
        
        # Flip edges connected to manifest indicators of exogenous latents:
        Edgelist2[Edgelist2[,2]%in%which(GroupVars$manifest&GroupVars$exogenous),] <- Edgelist2[Edgelist2[,2]%in%which(GroupVars$manifest&GroupVars$exogenous),2:1]
        
        # Flip all edges that are not connected to the bifactor:
        if (!missing(bifactor) && any(bifactor %in% Labels))
        {
          Edgelist2[!Labels[Edgelist2[,1]] %in% bifactor & !Labels[Edgelist2[,2]] %in% bifactor,1:2] <- Edgelist2[!Labels[Edgelist2[,1]] %in% bifactor & !Labels[Edgelist2[,2]] %in% bifactor,2:1]
        }
        
        iG <- graph.edgelist(Edgelist2)
        sp <- shortest.paths(iG,mode="out")
        sp[!is.finite(sp)] <- 0
        maxPaths <- apply(sp,1,max)
        # Mix in intercepts:
        if (any(GroupPars$edge=="int"))
        {
          maxPathsInts <- maxPaths[Edgelist[GroupPars$edge=="int",2]]
          if (!intAtSide)
          {
            maxPathsInts[maxPathsInts==min(maxPaths)] <- min(maxPaths) - 1
            maxPathsInts[maxPathsInts==max(maxPaths)] <- max(maxPaths) + 1
          }
          maxPaths <- c(maxPaths,maxPathsInts)
        }
        if (springLevels)
        {
          Cons <- cbind(NA,maxPaths)
          Layout <- qgraph.layout.fruchtermanreingold(Edgelist,vcount=length(maxPaths),constraints=Cons*sqrt(length(maxPaths)))
        } else {
          Layout <- cbind(NA,maxPaths)
          Layout[,1] <- ave(Layout[,2],Layout[,2],FUN=function(x)seq(-1,1,length=length(x)+2)[-c(1,length(x)+2)])
          
          # Mix intercepts:
          if (any(GroupPars$edge=="int"))
          {
            intMap <- rbind(manInts,latInts)
            for (i in sort(unique(Layout[,2])))
            {
              if (any(which(Layout[,2]==i)%in%intMap[,1]))
              {
                conInts <- which(Layout[,2]==i)
                conInts <- conInts[conInts%in%intMap[,1]]
                Layout <- mixInts(intMap[intMap[,1]%in%conInts,2],intMap,Layout,trim=TRUE,intAtSide=intAtSide)  
              }
              
            }
          }
        }
        
      } else Layout <- layout
      
      # loopRotation:
      if (layout%in%c("tree","tree2","tree3"))
      {
        loopRotation <- rep(0,nN)
        loopRotation[endoMan] <- pi
        loopRotation[exoMan] <- 0
        
        loopRotation[endoLat] <- 0
        noCons <- sapply(endoLat,function(x)nrow(Edgelist[(Edgelist[,1]==x|Edgelist[,2]==x) & (Edgelist[,1]%in%endoMan|Edgelist[,2]%in%endoMan),,drop=FALSE])==0)
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
        noCons <- sapply(exoLat,function(x)nrow(Edgelist[(Edgelist[,1]==x|Edgelist[,2]==x) & (Edgelist[,1]%in%exoMan|Edgelist[,2]%in%exoMan),,drop=FALSE])==0)
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
        
        
        if (any(GroupVars$exogenous) & optimizeLatRes)
        {        
          ### For latents that have loops, find a nice angle:
          for (i in which(Labels%in%latNames & Labels%in%GroupPars$lhs[GroupPars$lhs==GroupPars$rhs]))
          {
            # Layout subset of all connected:
            subEdgelist <- Edgelist[(Edgelist[,1]==i|Edgelist[,2]==i)&(Edgelist[,1]!=Edgelist[,2]),,drop=FALSE]              
            conNodes <- c(subEdgelist[subEdgelist[,1]==i,2],subEdgelist[subEdgelist[,2]==i,1])
            
            # Test for empty:
            if (nrow(subEdgelist)==0) conNodes <- sort(unique(c(Edgelist[,1:2])))
            
            subLayout <- Layout[conNodes,,drop=FALSE]
            
            # Add degree of edges passing node:
            lower <- which(Layout[,2] < Layout[i,2])
            higher <- which( Layout[,2] > Layout[i,2])
            passNode <- which((Edgelist[,1] %in% lower & Edgelist[,2] %in% higher) | (Edgelist[,2] %in% lower & Edgelist[,1] %in% higher))
            
            if (length(passNode) > 0)
            {
              passLayout <- do.call(rbind,lapply(passNode, function(ii)c(mean(Layout[Edgelist[ii,],1]), mean(Layout[Edgelist[ii,],2]))))
              subLayout <- rbind(subLayout, passLayout) 
            }
            
            Degrees <- apply(subLayout,1,function(x)atan2(x[1]-Layout[i,1],x[2]-Layout[i,2]))
            
            loopRotation[i] <- optimPoints[which.max(sapply(optimPoints,loopOptim,Degrees=Degrees))]
            
            # Completely forgot point of this whole thing here:
#             if (!grepl("lisrel",style,ignore.case=TRUE) | !any((Edgelist[,1]==i|Edgelist[,2]==i)&(Edgelist[,1]!=Edgelist[,2])&GroupPars$edge=="<->"))
#             {
# #               loopRotation[i] <- optimize(loopOptim,c(0,2*pi),Degrees=Degrees,maximum=TRUE)$maximum
#               loopRotation[i] <- optimPoints[which.max(sapply(optimPoints,loopOptim,c(0,2*pi),Degrees=Degrees))]
#             } else {
#               # Layout subset of all connected:
#               subEdgelist <- Edgelist[(Edgelist[,1]==i|Edgelist[,2]==i)&(Edgelist[,1]!=Edgelist[,2])&GroupPars$edge=="<->",]
#               conNodes <- c(subEdgelist[subEdgelist[,1]==i,2],subEdgelist[subEdgelist[,2]==i,1])
#               
#               # Test for empty:
#               if (nrow(subEdgelist)==0) conNodes <- sort(unique(c(Edgelist[,1:2])))
#               
#               subLayout <- Layout[conNodes,]
#               goodDegrees <- apply(subLayout,1,function(x)atan2(x[1]-Layout[i,1],x[2]-Layout[i,2]))
#               loopRotation[i] <- optimize(loopOptim,c(min(goodDegrees-pi/4),max(goodDegrees+pi/4)),Degrees=Degrees,maximum=TRUE)$maximum
#             }
          }
        }
        #     } else if (layout=="tree3"|layout=="circle3")
        #     {
        #       loopRotation <- rep(NA,nN)
        #       loopRotation[endoMan] <- pi
        #       loopRotation[exoMan] <- 0
      } else loopRotation <- rep(NA, length(Labels))
      
      ### ORDINALIZE LAYOUT ###
      if (layout=="tree")
      {
        Layout[Layout[,2]>0&Layout[,2]<5,2] <- as.numeric(as.factor(Layout[Layout[,2]>0&Layout[,2]<5,2]))
        Layout[Layout[,2]==0,2] <- (1*!residuals) * 0.25
        Layout[Layout[,2]==5,2] <- max(Layout[Layout[,2]<5,2]) + (1 - (1*!residuals)*0.25)
      }
      
      # Level layout:
      if (!missing(levels)&layout%in%c("tree","tree2","tree3","circle","circle2","circle3"))
      {
        if (length(levels)<length(unique(Layout[,2]))) stop(paste("'levels' argument must have at least",length(unique(Layout[,2])),"elements"))
        Layout[,2] <- levels[as.numeric(as.factor(Layout[,2]))]
      }
      
      ECP <- matrix(NA,nrow(Edgelist),2)
      
      # Set curves, edgeConnectPoints and rotate:    
      if (layout %in% c("tree","tree2","tree3"))
      {
        inBetween <- function(x)
        {
          if (Layout[x[1],2]!=Layout[x[2],2]) return(0) else return(sum(Layout[Layout[,2]==Layout[x[1],2],1] > min(Layout[x,1]) & Layout[Layout[,2]==Layout[x[1],2],1] < max(Layout[x,1])))
        }
        # Curves:
        inBet <- apply(Edgelist,1,inBetween)
        inBet[inBet>0] <- as.numeric(as.factor(inBet[inBet>0]))

#         if (!grepl("lisrel",style,ignore.case=TRUE) | all(!GroupVars$exogenous) | !residuals)
#         {
        if (isTRUE(curveAdjacent))
        {
          percurveAdjacent <- rep(TRUE,nrow(Edgelist))  
        } else {
          percurveAdjacent <- rep(FALSE,nrow(Edgelist))  
        
          curveAdjacent <- gsub("<->","cov",curveAdjacent)
          curveAdjacent <- gsub("(->)|(~>)","reg",curveAdjacent)
        
          if (is.character(curveAdjacent))
          {
            percurveAdjacent[(any(grepl("reg",curveAdjacent,ignore.case=TRUE))&GroupPars$edge%in%c("->","~>",ignore.case=TRUE))|(any(grepl("cov",curveAdjacent))&GroupPars$edge%in%c("<->"))] <- TRUE
          }
        }

        # Original curve:
          Curve <- ifelse(Layout[Edgelist[,1],2]==Layout[Edgelist[,2],2]&Edgelist[,1]!=Edgelist[,2]&GroupPars$edge!="int",ifelse(inBet<(1-percurveAdjacent),0,curve+curvature*(inBet)/max(1,max(inBet))*curve),NA)
#           Curve <- ifelse(Layout[Edgelist[,1],2]==Layout[Edgelist[,2],2]&Edgelist[,1]!=Edgelist[,2]&GroupPars$edge!="int",ifelse(inBet<(1-percurveAdjacent),0,curve),NA)
        
        
        #         } else {
#           Curve <- ifelse(Layout[Edgelist[,1],2]==Layout[Edgelist[,2],2]&Edgelist[,1]!=Edgelist[,2]&GroupPars$edge!="int" & Labels[Edgelist[,1]]%in%manNames & Labels[Edgelist[,2]]%in%manNames,ifelse(inBet<1,0,curve+inBet/max(inBet)*curve),NA)
#         }
        ### If origin node is "right" of  destination node, flip curve:
        Curve[Layout[Edgelist[,1],1] > Layout[Edgelist[,2],1]] <- -1 * Curve[Layout[Edgelist[,1],1] > Layout[Edgelist[,2],1]]
        ## If endo man, flip again:
        Curve <- ifelse(Edgelist[,1]%in%endoMan | Edgelist[,2]%in%endoMan, -1 *  Curve, Curve)
        
        
        ### Cardinal options:
        # Fuzzy matching:
        # exo/endo
        # man/lat
        # cov/reg/load
        # start/end
        
        if (any(grepl("all",cardinal)) || isTRUE(cardinal)) cardinal <- "exo endo man lat cov reg load source dest"
        
        ### Edge connect points:
        if (length(cardinal) > 0 && !identical(cardinal,FALSE) && !all(grepl("none",cardinal)))
        {
          
          if (packageDescription("qgraph")$Version == "1.2.3") warning("'cardinal' argument requires qgraph version 1.2.4")
          
          ECP <- matrix(NA,nrow(Edgelist),2)
          
          lvlDiff <- Layout[Edgelist[,1],2] - Layout[Edgelist[,2],2]
          
          ECP[lvlDiff>0,1] <- pi
          ECP[lvlDiff>0,2] <- 0
          
          ECP[lvlDiff<0,1] <- 0
          ECP[lvlDiff<0,2] <- pi
          
          
          ECP[lvlDiff==0,1] <- ifelse(Curve!=0, ifelse((loopRotation[Edgelist[,1]]+0.5*pi)%%2*pi <= pi, pi, 0), NA)[lvlDiff==0]
          ECP[lvlDiff==0,2] <- ifelse(Curve!=0, ifelse((loopRotation[Edgelist[,1]]+0.5*pi)%%2*pi <= pi, pi, 0), NA)[lvlDiff==0]
          
          
          ECP <- (ECP - 0.5 * (rotation-1) * pi ) %% (2*pi)
          
          # All nonspecified to NA:
          allSelect <- matrix(FALSE,nrow(ECP),2)
          for (cardGroup in cardinal)
          {
            select <- matrix(grepl("(exo)|(endo)|(man)|(lat)|(cov)|(reg)|(load)|(src)|(source)|(dest)",cardGroup),nrow(ECP),2)
            
            if (grepl("(endo)|(exo)",cardGroup))
            {
              
              # First node first / endo:
              select <- select & ((grepl("endo",cardGroup,ignore.case=TRUE) & !GroupVars$trueExo[Edgelist[,1]]) | 
                                    (grepl("exo",cardGroup,ignore.case=TRUE) & GroupVars$trueExo[Edgelist[,1]] )
              )
            }
            
            if (grepl("(lat)|(man)",cardGroup))
            {
              
              # Any node man / latent
              select <- select & ((grepl("lat",cardGroup,ignore.case=TRUE) & (!GroupVars$manifest[Edgelist[,1]] |  !GroupVars$manifest[Edgelist[,2]])) | 
                                    (grepl("man",cardGroup,ignore.case=TRUE) & (GroupVars$manifest[Edgelist[,1]] |  GroupVars$manifest[Edgelist[,2]]) )
              )
            }
            
            if (grepl("(cov)|(reg)|(load)",cardGroup))
            {        
              # Edge is cov/reg/loading:
              select <- select & ((grepl("cov",cardGroup,ignore.case=TRUE) & (GroupPars$edge=="<->")) | 
                                    (grepl("reg",cardGroup,ignore.case=TRUE) & (GroupPars$edge%in%c("->","~>") & !(!GroupVars$manifest[Edgelist[,1]] & GroupVars$manifest[Edgelist[,2]]))) | 
                                    (grepl("load",cardGroup,ignore.case=TRUE) & (GroupPars$edge%in%c("->","~>") & (!GroupVars$manifest[Edgelist[,1]] & GroupVars$manifest[Edgelist[,2]]))                         )
              )
            }
            
            if (grepl("(src)|(source)|(dest)",cardGroup))
            {                
              # Start/end:
              select[,1] <- select[,1] & grepl("(src)|(source)",cardGroup,ignore.case=TRUE)
              select[,2] <- select[,2] & grepl("dest",cardGroup,ignore.case=TRUE)
            }
            
            allSelect[select] <- TRUE
          }
          
          ECP[!allSelect] <- NA
          
        }

        
        ### Flip ECP, loopRotation and curve for any non bifactor variables, if specified:
        if (!missing(bifactor) && any(bifactor %in% Labels) && layout %in% c("tree2", "tree3", "circle2", "circle3"))
        {
          loopRotation[!Labels %in% bifactor] <- loopRotation[!Labels %in% bifactor] + pi
          ECP[!GroupPars$lhs%in%bifactor,1] <- ECP[!GroupPars$lhs%in%bifactor,1] + pi
          ECP[!GroupPars$lhs%in%bifactor,2] <- ECP[!GroupPars$lhs%in%bifactor,2] + pi
          Curve[!GroupPars$lhs%in%bifactor & !GroupPars$lhs%in%bifactor] <- -1 * Curve[!GroupPars$lhs%in%bifactor & !GroupPars$lhs%in%bifactor]
        }
        
        ### Rotate loopRotation:
        loopRotation <- loopRotation - 0.5 * (rotation-1) *pi
        
        
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
      
      
      ### ROTATE IF CIRCLE:
      if (layout%in%c("circle","circle2","circle3"))
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
      
      if (layout == "spring") loopRotation <- rep(NA, length(Labels))
      
 
      if (layoutSplit & length(unique(grSub)) > 1 & Sub > 0)
      {
        if (is.character(Layout))
        {
          Layout <- qgraph(Edgelist, layout = Layout, DoNotPlot = TRUE, edgelist=TRUE)$layout
        }
      ## Store in submodel list (could well be moved earlier but whatever)
      subModList[[Sub]] <- list(
        Layout = Layout,
        loopRotation = loopRotation,
        Curve = Curve,
        Labels = Labels,
        ECP = ECP
        )
        
      }
    }
    
    ### COMBINE SUB MODELS ###
    if (layoutSplit & length(unique(grSub)) > 1 & length(subModList) > 1)
    {
    
      ### Rescale subScale to height in width relative to diameter of device in inches ###
      din <- par("din")
      diamet <- sqrt(sum(din^2))
      subDim <- diamet * c(subScale, subScale2)
          
      if (is.character(Layout))
      {
        Layout <- qgraph(Edgelist, layout = Layout, DoNotPlot = TRUE)$layout
      }
      # Rescale main layout:
      Layout <- LayoutScaler(Layout,  din[1]/2, din[2]/2)
      subModList[[1]]$Layout <- LayoutScaler(subModList[[1]]$Layout)
      # Angle from center:
      centAngles <- atan2(subModList[[1]]$Layout[,1],subModList[[1]]$Layout[,2]) + pi
      
      subModList[[1]]$Layout <- LayoutScaler(subModList[[1]]$Layout,  din[1]/2, din[2]/2)
      
      centAngles[subModList[[1]]$Layout[,1]==0&subModList[[1]]$Layout[,2]==0] <- mean(centAngles[!(subModList[[1]]$Layout[,1]==0&subModList[[1]]$Layout[,2]==0)]) + pi
      if (layout %in% c('tree','tree2','tree3'))
      {
        err <- 1.1
      } else 
      {
        err <- 1.1
      }
      srot <- ifelse(rotation%in%c(1,3),1/err,err)

      centAngles <- atan2(srot*sin(centAngles),cos(centAngles))
      if (subRes != 0)
      {
        centAngles <- round_any(centAngles%%(2*pi), (2*pi)/subRes)  
      }      
      # Rescale and rotate sub layouts and enter in main layout:
      for (g in rev(seq_along(subModList)))
      {
        if (g > 1 && !is.null(subModList[[g]]))
        {
          link <- c(which(subModList[[1]]$Labels == subLinks[g-1]), which(subModList[[g]]$Labels == subLinks[g-1])  )
          
          # Scale:
#           subDim2 <-  abs(c(RotMat(centAngles[link[1]]) %*% subDim))
#           subDim2 <- subDim2 / abs(c(RotMat(centAngles[link[1]]) %*% rev(din)))
                    
#           subModList[[g]]$Layout <- LayoutScaler(subModList[[g]]$Layout, c(-1,1) * subDim[1]/2, c(-1,1) * subDim[2]/2)
          
          # Map to inch coordinates:
          subModList[[g]]$Layout <- LayoutScaler(subModList[[g]]$Layout, subDim[1]/2, subDim[2]/2)
          
          # Center to link:
          subModList[[g]]$Layout[,1] <- subModList[[g]]$Layout[,1] - subModList[[g]]$Layout[link[2],1]
          subModList[[g]]$Layout[,2] <- subModList[[g]]$Layout[,2] - subModList[[g]]$Layout[link[2],2]
                    
          # Rotate:              
          subModList[[g]]$Layout <-  t(RotMat(centAngles[link[1]]) %*% t(subModList[[g]]$Layout))
        
          # Map back to usr coordinates:
          #           subModList[[g]]$Layout[,1] <- subModList[[g]]$Layout[,1] / (din[1]/2)
          #           subModList[[g]]$Layout[,2] <- subModList[[g]]$Layout[,2] / (din[2]/2)
          
          # Center to Layout:
          subModList[[g]]$Layout[,1] <- subModList[[g]]$Layout[,1] + subModList[[1]]$Layout[link[1],1]
          subModList[[g]]$Layout[,2] <- subModList[[g]]$Layout[,2] + subModList[[1]]$Layout[link[1],2]
        } 
          # Enter in general model:
          subLabnums <- match(subModList[[g]]$Labels[subModList[[g]]$Labels!='1'],Labels)
          subLabnums <- c(subLabnums,manInts[match(match(GroupPars$rhs[GroupPars$sub==g & GroupPars$edge == "int" & GroupPars$rhs %in% manNames],Labels),manInts[,2]),1])
          subLabnums <- c(subLabnums,latInts[match(match(GroupPars$rhs[GroupPars$sub==g & GroupPars$edge == "int" & GroupPars$rhs %in% latNames],Labels),latInts[,2]),1])
                  
          Layout[subLabnums,] <- subModList[[g]]$Layout
          Curve[GroupPars$sub == g] <- subModList[[g]]$Curve
        
        if (g == 1)
        {
          loopRotation[subLabnums] <- subModList[[g]]$loopRotation
          ECP[object@Pars$sub == g,]  <- subModList[[g]]$ECP
          
          for (g2 in length(subModList):2)
          {
            if (!is.null(subModList[[g2]]))
            {
              loopRotation[Labels==latNames[g2-1]] <- centAngles[g2-1]
#               ECP[object@Pars$sub == g,]  <- centAngles[g2-1]
            }
          }
          
        } else 
        {
          loopRotation[subLabnums] <- (subModList[[g]]$loopRotation + centAngles[link[1]]) %% ( 2*pi)
          ECP[object@Pars$sub == g,] <- (subModList[[g]]$ECP + centAngles[link[1]]) %% ( 2*pi)
        }
        
      }
      
    }

    # Edge labels:
    if (edge.labels)
    {
      eLabels <- GroupPars$label
    } else eLabels <- rep("",nrow(Edgelist))
    
    # vsize:
    vSize <- numeric(nN)
    vSize[Labels%in%manNames] <- sizeMan
    vSize[Labels%in%latNames] <- sizeLat
    vSize[Labels=="1"] <- sizeInt
    
    vSize2 <- numeric(nN)
    vSize2[Labels%in%manNames] <- sizeMan2
    vSize2[Labels%in%latNames] <- sizeLat2
    vSize2[Labels=="1"] <- sizeInt2
    
    eColor <- rep(NA,nrow(Edgelist))
    #     tColor <- rep(rgb(0.5,0.5,0.5),nrow(GroupThresh))
    if (missing(thresholdColor)) 
    {
      tColor <- rep("border", nN) 
    } else {
      tColor <- rep(thresholdColor, nN)
    }
    
    ### WHAT TO PLOT? ###
    if (grepl("path|diagram|mod",what,ignore.case=TRUE))
    {
      
    } else if (grepl("stand|std",what,ignore.case=TRUE))
    {
      Edgelist <- cbind(Edgelist,GroupPars$std)
      if (edge.labels) eLabels <- GroupPars$std
    } else if (grepl("est|par",what,ignore.case=TRUE))
    {
      Edgelist <- cbind(Edgelist,GroupPars$est)
      if (edge.labels) eLabels <- GroupPars$est
    } else if (grepl("eq|cons",what,ignore.case=TRUE))
    {
      #       eColor <- rep(rgb(0.5,0.5,0.5),nrow(Edgelist))
      unPar <- unique(object@Pars$par[object@Pars$par>0 & duplicated(object@Pars$par)])
      if (pastel)
      {
        cols <- rainbow_hcl(max(c(object@Pars$par,GroupThresh$par)), c = 35, l = 85)
      } else {
        cols <- rainbow(max(c(object@Pars$par,GroupThresh$par)))
      }
      for (i in unPar)
      {
        eColor[GroupPars$par==i] <- cols[i]
      }
      if (nrow(GroupThresh) > 0)
      {
        warning("Equality constraints of Thresholds currently not supported")
#         for (i in 1:nrow(GroupThresh))
#         {
#           if (GroupThresh$par[i]>0 & sum(GroupThresh$par[i] == object@Thresholds$par) > 1 )
#           {
#             tColor[i] <- cols[GroupThresh$par[i]]
#           }
#         }
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
          if (mode(NodeGroups[[g]])=="character") 
          {
            ### hier grepl!
            NodeGroups[[g]] <- matchVar(NodeGroups[[g]], GroupVars, manIntsExo, manIntsEndo, latIntsExo, latIntsEndo)
              #             NodeGroups[[g]] <- match(NodeGroups[[g]],Labels)
          }
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
      if (inheritColor)
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
        eLabels <- GroupPars$label
      } else if (grepl("stand|std",whatLabels,ignore.case=TRUE))
      {
        eLabels <- GroupPars$std
      } else if (grepl("est|par",whatLabels,ignore.case=TRUE))
      {
        eLabels <- GroupPars$est
      } else if (grepl("eq|cons",whatLabels,ignore.case=TRUE))
      {
        eLabels <- GroupPars$par
      } else if (grepl("no|omit|hide|invisible",whatLabels,ignore.case=TRUE))
      {
        eLabels <- rep("",nrow(Edgelist))
      } else stop("Could not detect use of 'whatLabels' argument")
    }
    
    # Abbreviate:
    if (!"edges"%in%as.expression)
    {
      if (is.numeric(eLabels))
      {
        eLabels <- ifelse(is.na(eLabels),"",formatC(eLabels, format=ifelse(all(eLabels%%1==0),'d','f'), digits=nDigits))
      } else 
      {
        if (nCharEdges>0)  eLabels <- abbreviate(eLabels,nCharEdges)
      }
    }
    if (!"nodes"%in%as.expression)
    {
      if (is.numeric(Labels))
      {
        Labels <- ifelse(is.na(Labels),"",formatC(Labels, format=ifelse(all(Labels%%1==0),'d','f'), digits=nDigits))
      } else 
      {
        if (nCharNodes>0 )  Labels <- abbreviate(Labels,nCharNodes)
      }
    }
    
    #     ### CONVERT TO LISREL STYLE ###
    if (grepl("lisrel",style,ignore.case=TRUE) & residuals & covAtResiduals)
    {
      isResid <- GroupPars$edge == "<->" & GroupPars$lhs != GroupPars$rhs & (GroupPars$lhs %in% GroupPars$lhs[GroupPars$edge == "<->" & GroupPars$lhs == GroupPars$rhs] & GroupPars$rhs %in% GroupPars$rhs[GroupPars$edge == "<->" & GroupPars$lhs == GroupPars$rhs])
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

      if (!layoutSplit)
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
#         if (length(exoMan)==0)
#         {
#           Mar[3] <- Mar[3] + 3
#         }
        
        # Rotate:
        Mar <- Mar[(0:3 + (rotation-1)) %% 4 + 1]
        
        # Add 2 to top for title:
        if (title) Mar[3] <- Mar[3] + 2
      } else 
      {
        Mar <- c(3,3,3,3)
      }
    } else Mar <- mar
    
    # Overwrite edge colors if appropriate:
    if (!missing(edge.color))
    {
      eColor <- edge.color
      if (thresholds & missing(thresholdColor))
      {
        tColor <- rep(edge.color,length(tColor))
      }
    }
    
    # Fixed and free edges:
    if (length(freeStyle) > 2 | length(fixedStyle) > 2) warning("'freeStyle' and 'fixedStyle' are assumed to be vectors of at most length 2. Unexpected results will probably occur.")
    # lty:
    lty <- rep(1,nrow(GroupPars))
    
    # fixedStyle
    if (any(is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle))) lty <- ifelse(GroupPars$fixed,as.numeric(fixedStyle[is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle)]),lty) 
    
    if (any(isColor(fixedStyle) & !(is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle)))) eColor[GroupPars$fixed] <- fixedStyle[isColor(fixedStyle) & !(is.numeric(fixedStyle) | grepl("^\\d+$",fixedStyle))]
    
    
    # freeStyle:
    if (any(is.numeric(freeStyle) | grepl("\\d+",freeStyle))) lty <- ifelse(GroupPars$fixed,lty,as.numeric(freeStyle[is.numeric(freeStyle) | grepl("\\d+",freeStyle)]))
    
    if (any(isColor(freeStyle) & !(is.numeric(freeStyle) | grepl("\\d+",freeStyle)))) eColor[!GroupPars$fixed] <- freeStyle[isColor(freeStyle) & !(is.numeric(freeStyle) | grepl("\\d+",freeStyle))]
    
    # Directed settings:
    
    Directed <- GroupPars$edge!="--"
    
    # Convert labels to expressions:
    if ("edges"%in%as.expression)
    {
      eLabels <- lapply(eLabels,function(x)if (x=="") x else as.expression(parse(text=x)))
    }    # Convert labels to expressions:
    if ("nodes"%in%as.expression)
    {
      Labels <- lapply(Labels,function(x)if (x=="") x else as.expression(parse(text=x)))
    }
    # Restore layout function:
    if (!is.null(layoutFun)) Layout <- layoutFun
    
    # Overwrite node and edge labels:
    if (!missing(nodeLabels))
    {
      nLab <- nodeLabels[object@Vars$name %in% GroupVars$name]
    } else nLab <- Labels
    # Overwrite node and edge labels:
    if (!missing(edgeLabels))
    {
      eLab <- edgeLabels[object@Pars$group==gr]
    } else eLab <- eLabels
    
    
    ### WITHIN - BETWEEN FRAMEWORK ###
    CircleEdgeEnd <- rep(FALSE, nrow(Edgelist))
    
    if (any(c('Between','Within')%in%GroupPars$BetweenWithin))
    {
      if (all(GroupPars$BetweenWithin == 'Within'))
      {
        # WITHIN CLUSTER SETUP
        BetweenPars <- object@Pars[object@Pars$group == gsub("Within$","Between",gr),]
#         BetweenInts <- BetweenPars$rhs[BetweenPars$edge == 'int']
        BetweenVars <- unique(c(BetweenPars$lhs,BetweenPars$rhs))
        CircleEdgeEnd[GroupPars$rhs %in% BetweenVars & GroupPars$edge %in% c('->','~>')] <- TRUE

      } else if (all(GroupPars$BetweenWithin == 'Between'))
      {
        # BETWEEN CLUSTER SETUP
        WithinPars <- object@Pars[object@Pars$group == gsub("Between$","Within",gr),]
        WithinVars <- unique(c(WithinPars$lhs,WithinPars$rhs))
        Shape[Labels %in% WithinVars ] <- shapeLat
        
      } else stop("BetweenWithin not only 'Between' or 'Within'.")
    }
    
    ### Threshold setup ###
  
    bars <- list()
    length(bars) <- nN
    barSide <- rep(1, nN)
    
    if (thresholds)
    {
      if (missing(thresholdColor) & missing(edge.color)) 
      {
        tColor <- rep("border", nN)
      }
      if (nrow(GroupThresh) > 0)
      {
        for (node in unique(match(GroupThresh$lhs,GroupVars$name)))
        {
#           node <- which(Labels==GroupThresh$lhs[i])
          # Compute side:
          IntSide <- 1
          if (layout=="tree")
          {
            if (rotation%in%c(1,3))
            {
              barSide[node] <- ifelse(Layout[node,2]>mean(Layout[,2]),3,1)
            } else {
              barSide[node] <- ifelse(Layout[node,1]>mean(Layout[,1]),4,2)
            }
          } else {
            barSide[node] <- sum((atan2(scale(Layout[,1])[node],scale(Layout[,2])[node])+pi)%%(2*pi) > c(0,pi/2,pi,1.5*pi))
          }
          bars[[node]] <- pnorm(GroupThresh$est[GroupThresh$lhs == GroupVars$name[node]])
        }
      }
    }
    
    # curveScale:
    curveScale <- ! layout %in% c('tree','tree2','tree3')     
#     curveScale <- TRUE
    
    ### RUN QGRAPH ###
    
    qgraphRes[[which(Groups==gr)]] <- qgraph(Edgelist,
                                             labels=nLab,
                                             bidirectional=Bidir,
                                             directed=Directed,
                                             shape=Shape,
                                             layout=Layout,
                                             lty=lty,
                                             loopRotation=loopRotation,
                                             curve=Curve,
                                             edge.labels=eLab,
                                             mar=Mar,
                                             vsize = vSize,
                                             vsize2 = vSize2,
                                             edge.color=eColor,
                                             groups=NodeGroups2,
                                             color=Vcolors,
                                             residuals=LoopAsResid,
                                             residScale = residScale,
                                             residEdge = isResid,
                                             edgelist = TRUE,
                                             curveDefault = curveDefault,
                                             knots = GroupPars$knot,
#                                              curvePivot = curvePivot,
                                             aspect = layoutSplit,
                                             CircleEdgeEnd = CircleEdgeEnd,
                                             curveScale = curveScale,
                                             bars = bars,
                                             barSide = barSide,
                                             barColor = tColor,
                                             barLength = thresholdSize,
                                             barsAtSide = ThreshAtSide,
                                             edge.label.cex = edge.label.cex,
                                             edgeConnectPoints = ECP,
                                             ...)

#     if (thresholds)
#     {
#       # Overwrite color to white if bg is dark (temporary solution)
#       if (missing(thresholdColor) & missing(edge.color)) 
#       {
#         if (mean(col2rgb(qgraphRes[[which(Groups==gr)]]$plotOptions$background)/255) <= 0.5) tColor <- rep("white",length(tColor))
#       }
#       if (nrow(GroupThresh) > 0)
#       {
#         for (i in 1:nrow(GroupThresh))
#         {
#           node <- which(Labels==GroupThresh$lhs[i])
#           # Compute side:
#           IntSide <- 1
#           if (layout=="tree")
#           {
#             if (rotation%in%c(1,3))
#             {
#               IntSide <- ifelse(Layout[node,2]>mean(Layout[,2]),3,1)
#             } else {
#               IntSide <- ifelse(Layout[node,1]>mean(Layout[,1]),4,2)
#             }
#           } else {
#             IntSide <- sum((atan2(qgraphRes[[which(Groups==gr)]]$layout[node,1],qgraphRes[[which(Groups==gr)]]$layout[node,2])+pi)%%(2*pi) > c(0,pi/2,pi,1.5*pi))
#           }
#           IntInNode(qgraphRes[[which(Groups==gr)]]$layout[node,,drop=FALSE],vSize[node],Shape[node],pnorm(GroupThresh$est[i]),width=0.5,triangles=FALSE,col=tColor[i],IntSide,!ThreshAtSide)
#         }
#       }
#     }
    
    if (title)
    {
      #       if (length(Groups)==1) title("Path Diagram",line=3) else title(paste0("Path Diagram for group '",gr,"'"),line=3)
      title(gr, col.main=title.color, adj = title.adj, outer = TRUE, cex.main = title.cex, line = title.line)
    }
  }
  par(ask=askOrig)
  if (length(qgraphRes)==1) qgraphRes <- qgraphRes[[1]]
  invisible(qgraphRes)
}
