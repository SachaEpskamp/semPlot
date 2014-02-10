\name{semPaths}
\alias{semPaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot path diagram for SEM models.
}
\description{
This function creates a path diagram of a SEM model (or general linear model), which is then plotted using \code{\link[qgraph]{qgraph}}. Currently many different SEM programs and packages are supported. Please see my website (www.sachaepskamp.com) for more details on which packages are supported and what is supported for each package.
}
\usage{
semPaths(object, what = "paths", whatLabels, style, layout = "tree", 
  intercepts = TRUE, residuals = TRUE, thresholds = TRUE, intStyle = "multi", 
  rotation = 1, curve, curvature = 1, nCharNodes = 3, nCharEdges = 3, sizeMan = 5,
   sizeLat = 8,  sizeInt = 2, sizeMan2, sizeLat2, sizeInt2, shapeMan, shapeLat, 
  shapeInt = "triangle", ask, mar, title, title.color = "black", title.adj = 0.1, 
  title.line = -1, title.cex = 0.8, include, combineGroups = FALSE, manifests, 
  latents, groups, color, residScale, gui = FALSE, allVars = FALSE, edge.color, 
  reorder = TRUE, structural = FALSE, ThreshAtSide = FALSE, thresholdColor, 
  thresholdSize = 0.5, fixedStyle = 2, freeStyle = 1, 
  as.expression = character(0), optimizeLatRes = FALSE, inheritColor = TRUE, 
  levels, nodeLabels, edgeLabels, pastel = FALSE, rainbowStart = 0, intAtSide, 
  springLevels = FALSE, nDigits = 2, exoVar, exoCov = TRUE, centerLevels = TRUE, 
  panelGroups = FALSE, layoutSplit = FALSE, measurementLayout = "tree", subScale, 
  subScale2, subRes = 4, subLinks, modelOpts = list(),
  curveAdjacent = '<->', edge.label.cex = 0.6,  cardinal = "none", equalizeManifests = FALSE, 
  covAtResiduals = TRUE, bifactor, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{object}{
A \code{"semPlotModel"} object or any of the input types that can be used in \code{\link{semPlotModel}} directly. 
}

  \item{what}{
What should the edges indicate in the path diagram? This function uses \code{\link{grepl}} to allow fuzzy matching and is not case sensitive. E.g., \code{par} will also match \code{Parameters}.
  \describe{
    \item{\code{path}, \code{diagram} or \code{mod}}{This will display the model as an unweighted network (gray edges by default).} 
      \item{\code{est} or \code{par}}{This will display the parameter estimates as weighted edges.}
    \item{\code{stand} or \code{std}}{This will display the standardized parameter estimates, if available, as weighted edges.}
      \item{\code{eq} or \code{cons}}{This is the same graph as \code{path}. except that parameters with equality constraints are now colored. Parameters with the same color are constrained to be equal.}
          \item{\code{col}}{This will create an unweighted graph of the path diagram, where edges are colored with a mix of the colors of connected nodes.}
  }
  }
  
    \item{whatLabels}{
What should the edge labels indicate in the path diagram? This function uses \code{\link{grepl}} to allow fuzzy matching and is not case sensitive. E.g., \code{par} will also match \code{Parameters}. Default depends on the \code{what} argument, defaulting to the respective elements in the list below for values of \code{what} in the list above.
  \describe{
    \item{\code{name}, \code{label}, \code{path} or \code{diagram}}{This will display the edge names as labels.} 
    \item{\code{est} or \code{par}}{This will display the parameter estimate in edge labels.}
    \item{\code{stand} or \code{std}}{This will display the standardized parameter estimate in edge labels.}
        \item{\code{eq} or \code{cons}}{This will display the parameter number in edge labels. 0 indicates the parameter is fixed, parameters with the same parameter number are constrained to be equal.}
        \item{\code{no}, \code{omit}, \code{hide} or \code{invisible}}{Hides edge labels.}
     }
}
  

  \item{style}{
The style to use. Currently only indicates what the (residual) variances look like. Use \code{"ram"}, \code{"mx"} or \code{"OpenMx"} for double headed selfloops and \code{"lisrel"} for single headed edges with no node as origin. Defaults to \code{"ram"} unless the input is a lisrel model.
}

  \item{layout}{
A string indicating how the nodes should be placed. Similar to the 'layout' argument in \code{\link[qgraph]{qgraph}}. Can be one of the following strings.
  \describe{
    \item{tree}{The integrated tree-like layout. Places exogenous variables at the top and endogenous variables at the bottom. See 'details' for more details.}
    \item{circle}{The same layout as "tree", except that afterwards the horizontal levels of the layout are placed in circles. Especially useful for models with a large number of manifest variables and a relatively small number of latent variables.}
    \item{spring}{Calls the "spring" layout in \code{\link[qgraph]{qgraph}}, which uses the Fruchterman-reingold algorithm (Fruchterman & Reingold, 1991).}
    \item{tree2}{Calls the \code{\link[igraph]{layout.reingold.tilford}} function from the igraph package (Csardi & Nepusz, 2006), which uses the Reingold-Tilford algorithm (Reingold & Tilford, 1981). Before calling the algorithm roots are chosen and a slightly modified version of the graph is used to produce consistent results. See 'details'.}
    \item{circle2}{The same layout as "tree2", except that afterwards the horizontal levels of the layout are placed in circles.}
    \item{Other options}{If the assigned value is not in this list it is sent to \code{\link[qgraph]{qgraph}}. This allows for manual specification of the layout as well as using functions found in the 'igraph;  library.}
     }
}
  \item{intercepts}{
  Logical, should intercepts be included in the path diagram?
}

  \item{residuals}{
  Logical, should residuals (and variances) be included in the path diagram?
}

  \item{thresholds}{
  Logical, should thresholds be included in the path diagram?
}

  \item{intStyle}{
  Style of the intercepts. \code{"multi"} plots a separate unit vector node for each intercept and \code{"single"} plots a single unit vector node. Currently, \code{"single"} is not well supported and might lead to unexpected results.
}
  \item{rotation}{
  An integer indicating the rotation of the layout when "tree" or "tree2" layout is used. 1, 2, 3 and 4 indicate that exogenous variables are placed at the top, left side, bottom and right side respectively.
}
  \item{curve}{
  The curvature of the edges. In tree layouts this argument only curves the edges that are between nodes on the same level. e.g., correlations between exogenous manifest variables.
}
  \item{curvature}{Sets the strength of scaling in curvature for curved edges at the same horizontal level in tree layouts. The curve will be set to \code{curve + curvature * n / max(n)}, where \code{n} is the number of nodes in between the two connected nodes.}
  \item{nCharNodes}{
 Number of characters to abbreviate node labels to (using \code{\link[base]{abbreviate}}). Set to 0 to omit abbreviation.
}
  \item{nCharEdges}{
 Number of characters to abbreviate edge labels to (using \code{\link[base]{abbreviate}}). Set to 0 to omit abbreviation.
}
  \item{sizeMan}{
  Width of the manifest nodes, sent to the 'vsize' argument in \code{\link[qgraph]{qgraph}}.
}
  \item{sizeLat}{
  Width of the latent nodes, sent to the 'vsize' argument in \code{\link[qgraph]{qgraph}}.
}
  \item{sizeInt}{
  Width of the unit vector nodes, sent to the 'vsize' argument in \code{\link[qgraph]{qgraph}}.
}
  \item{sizeMan2}{
  Height of the manifest nodes, sent to the 'vsize2' argument in \code{\link[qgraph]{qgraph}}.
}
  \item{sizeLat2}{
  Height of the latent nodes, sent to the 'vsize2' argument in \code{\link[qgraph]{qgraph}}.
}
  \item{sizeInt2}{
  Height of the unit vector nodes, sent to the 'vsize2' argument in \code{\link[qgraph]{qgraph}}.
}
  \item{shapeMan}{
  Shape of the manifest nodes, sent to the 'shape' argument in \code{\link[qgraph]{qgraph}}. Defaults to \code{"square"} or \code{"rectangle"} if width and height differ.
}
  \item{shapeLat}{
  Shape of the latent nodes, sent to the 'shape' argument in \code{\link[qgraph]{qgraph}}. Defaults to \code{"circle"} or \code{"ellipse"} if width and height differ.
}
  \item{shapeInt}{
  Shape of the constant nodes, sent to the 'shape' argument in \code{\link[qgraph]{qgraph}}. Defaults to \code{"triangle"}.
}
  \item{ask}{
Specifies the 'ask' parameter in \code{\link[graphics]{par}}. Defaults to TRUE if multiple groups are in the model.
}
  \item{mar}{
Same as the 'mar' argument in \code{\link[qgraph]{qgraph}}. By default this argument is based on the values of 'rotation', 'style' and 'title'.
}
  \item{title}{
Logical, should titles be plotted of the group names above each plot?
}
  \item{title.color}{
Color of the titles.
}
  \item{title.adj}{
Adjustment of title as used by \code{'adj'} in \code{par}.
}
  \item{title.line}{
Line of title as used by \code{'line'} in \code{title}.
}
  \item{title.cex}{
Size of title as used by \code{'cex.main'} in \code{par}.
}
  \item{include}{
Integer vector indicating which groups should be included in the output. e.g., to only plot a diagram for the first group use \code{include = 1}.
}
  \item{combineGroups}{
Logical. If TRUE all groups are combined in the same path diagram.
}
  \item{manifests}{
A character vector in which every element is the name of a manifest variable in the model. This argument can be used to overwrite the order in which nodes are plotted in the graph if \code{reorder = FALSE}
}
  \item{latents}{
A character vector in which every element is the name of a latent variable in the model. This argument can be used to overwrite the order in which nodes are plotted in the graph if \code{reorder = FALSE}
}
  \item{groups}{
  Groups nodes that should be colored the same, similar to the 'groups' argument in \code{\link[qgraph]{qgraph}} with a few exceptions. Should be a list containing in each element the names (instead of numbers as in qgraph) of nodes that belong together. Nodes that are indicated to belong to a group will be assigned the same color, as given by the 'color' argument. Nodes not belonging to a group will be assigned the color "", which indicates that they will inherit a mix of the colors of connected nodes (or white, if no connected nodes are colored.)
  
  In addition, this argument can be assigned a single character: "manifests", "latents" or "both" to make a single group for each manifest, latent or both manifest and latent variables. e.g., \code{groups = "latents"} will color each latent variable uniquely, and color all manifest variables a mixture of the colors of latents they load on.
}
  \item{color}{
Controls the color of nodes. Similar to 'color' in \code{\link[qgraph]{qgraph}}. A color vector indicating the color for each group, a single color character indicating the color for all nodes or a color vector indicating the color for each node separately. Can also be a list contaning one or more of the following elements (using fuzzy matching):
\describe{
\item{man}{The colors for manifest nodes}
\item{lat}{The colors for latent nodes}
\item{int}{The color for intercepts}
}
}
  \item{residScale}{
The size of residual edges if \code{style = "lisrel"}. Defaults to two times the value of 'sizeMan'.
}
  \item{gui}{
Not yet implemented.
}
  \item{allVars}{
Logical. If TRUE all variables are plotted in the path diagrams for each group. If FALSE only variables are plotted that are used in the group.
}
  \item{edge.color}{
A value indicating the color of all edges or a vector indicating the color of each edge. Useful for manually overwriting edge colors. 
}
  \item{reorder}{
Logical. Should manifest variables be reordered to be near latent factors they are most connected to in the "tree" layout? If FALSE manifest variables are placed in the order they appear in the Pars.
}
  \item{structural}{
Logical. Set this to FALSE to only show the structural model (omit all manifest variables.)
}
  \item{ThreshAtSide}{
 Logical. If TRUE, thresholds are plotted as small lines at the side of manifest nodes, otherwise they are plotted as lines inside the nodes.
}
  \item{thresholdColor}{
Color of the threshold lines. Defaults to "black"
}

  \item{thresholdSize}{
Size of threshold bars relative to the size of the node. 
}

  \item{fixedStyle}{
A vector of length one or two specifying the color and line type (same as 'lty' in \code{\link[graphics]{par}}) of fixed parameters. Can be both character and numeric. If one of the elements encodes a color it is used to overwrite the color of fixed edges, and if an element can be coerced to a numeric it is used to encode the line type.

For example, \code{ fixedStyle = c("red",3)} specifies that all fixed parameters should be visualized with a red edge with \code{lty=3}
}
  \item{freeStyle}{
Same as 'fixedStyle' but for free parameters instead.
}
  \item{as.expression}{
A character vector indicating which labels should be treated as an \code{\link[base]{expression}}, so that mathematical notation and Greek letters can be used in the path diagram. If this vector contains \code{"nodes"} all node labels are converted to expressions, and if this vector contains \code{"edges"} all node labels are converted to expressions. Defaults to \code{"edges"} only if the input is a Lisrel model.
}
  \item{optimizeLatRes}{
Logical. If this is TRUE, the angle of the incoming residuals on latent variables is attempted to be optimally chosen so its position conflicts with the least amount of connected edges.
}
  \item{inheritColor}{
Logical, should uncolored nodes obtain a mix of connected colored nodes? Defaults to \code{TRUE}.
}
\item{levels}{A numeric vector usually of length 4. Controls the relative vertical position of variable levels (exogenous and endogenous latents and manifests) under default rotation in tree and circle layouts. This can be used to control the spacing between these levels. e.g., \code{c(1,5,6,7)} will create more space between endogenous manifests and latents.}

\item{nodeLabels}{A vector or list to manually overwrite the node labels. Can include expressions.}

\item{edgeLabels}{A vector or list to manually overwrite the edge labels. Can include expressions.}

\item{pastel}{Logical, should default colors (for groups or edge equality constraints) be chosen from pastel colors? If TRUE then \code{\link[colorspace]{rainbow_hcl}} is used.}

\item{rainbowStart}{A number between 0 and 1 indicating the offset used in rainbow functions for default node coloring.}
  \item{exoVar}{
Should variances of truely exogenous variables (no incomming directed edge) be plotted? Defaults to \code{TRUE} unless \code{style = "lisrel"}.
}
  \item{intAtSide}{
Logical to control if intercepts should be plotted to the side of manifest nodes or at the bottom/top. Defaults only to FALSE if 'residuals=FALSE'.
}
  \item{springLevels}{
Logical indicating if the placement on horizontal levels with \code{tree3} layout should be determined by a force embedded algorithm.
}
  \item{nDigits}{
Number of digits to round numeric values to.
}

  \item{exoCov}{
Should covariances between truely exogenous variables (no incomming directed edge) be plotted? Defaults to \code{TRUE}.
}
 \item{centerLevels}{
Only used if \code{layout} is set to \code{"tree2"}, should each level be centered? Defaults to \code{TRUE}
}
  \item{panelGroups}{
Logical to automatically create a panel plot of multiple group models. Defaults to FALSE.
}
  \item{layoutSplit}{
Logical that can be used to split computing of layout between structural and measurment models. This is very useful in more complicated models where the structural part is best shown by using a spring layout.
}
  \item{measurementLayout}{
Logical indicating the layout algorithm to use for measurement models if \code{layoutSplit = TRUE} (the structural model will obtain a layout given by the \code{layout} argument).
}
  \item{subScale}{
Width of submodels (measurment models) if \code{layoutSplit = TRUE}.
}
  \item{subScale2}{
Height of submodels (measurment models) if \code{layoutSplit = TRUE}.
}
  \item{subRes}{
Integer indicating the resolution of which measurment models can be rotated around their corresponding latent variable. The default, 4, indicates that they can be placed only to polar coordinates. Set to 360 to allow every angle of rotation.
}
  \item{subLinks}{
Vector of variables to link to. Currently not well supported so avoid using this argument.
}

  \item{modelOpts}{
A lists containing arguments sent to \code{\link{semPlotModel}} in case the input is not of class \code{semPlotModel}.
}

  \item{curveAdjacent}{What edges between adjacent horizontal nodes be curved? Can be \code{'<->'} or \code{'cov'} to indicate bidirectional covariances, \code{'->'} or \code{'reg'} for directed regressions or a vector containing both.}
  
  \item{edge.label.cex}{Controls the font size of the edge labels. Same as in \code{\link[qgraph]{qgraph}} except that the default is now 0.8.}

  \item{cardinal}{Should edges in a tree layout connect to the four cardinal points of one of the borders of the node rather than point to the center of the node? Can be set to \code{TRUE} or \code{"all"} to enamble this behavior for all edges and \code{FALSE} or \code{"none"} to disable this behavior for all edges. Alternatively a vector with strings can be specified in which each string specifies a certain group of edges. Fuzzy matching is used on the strings \code{"exo"} for edges with the first node being exogenous (or indicator of exogenous latent), \code{endo} for edges with first node being endogeonous, \code{manifest} for edges connected to any manifest node, \code{latent} for edges connected to any latent node, \code{cov} for covariances, \code{reg} for regressions, \code{load} for factor-loadings, \code{source} for only the start of an edge and \code{end} for only the end of a node. These strings can be combined at will. For example, \code{cardinal = c("exo cov","load end")} (the default) or equivelantly \code{cardinal = c("exogenous covariances","source of loadings")} will only cardenalize the edges that repressent exogenous covariances or the end of factor loadings. }

\item{equalizeManifests}{Logical. Should the distances between manifest nodes in the \code{tree1} layout be equalized? Defaults to \code{TRUE}}

\item{covAtResiduals}{Logical, should covariances be drawn at the start of residuals when \code{style="lisrel"} is used? Defaults to TRUE.}


\item{bifactor}{
A string vector containing the name(s) of the general bifactor(s). This will automatically create a bifactor plot.
}


  \item{\dots}{
  Arguments sent to the \code{\link[qgraph]{qgraph}} function. These arguments can further control the output of the graph. Some usefull arguments in drawing path diagrams are:
  \describe{
  \item{edge.width}{
    Scales the edge width and arrow size of the plot. These can also be manually set using 'esize' and 'asize'.
  }
    \item{node.width}{
    Scales the width of nodes and also the height if shapes circle and square are used. Can also be a vector with scalar for each node.
  }
    \item{node.height}{
    Scales the height of nodes. Can also be a vector with scalar for each node. Not used with circle and square shapes.
  }
  
    \item{esize}{Size of the largest edge (or what it would be if there was an edge with weight maximum). Defaults to: max((-1/72)*(nNodes)+5.35,1) for weighted graphs and 2 for unweighted graphs. In directed graphs these values are halved.}
    
    \item{asize}{Size of the arrowhead. Defaults to 2 for graphs with more than 10 nodes and 2 to smaller graphs.}
  
  \item{minimum}{Edges with absolute weights under this value are omitted. Defaults to 0 for graphs with less than 50 nodes or 0.1 for larger graphs.}
\item{maximum}{qgraph regards the highest of the maximum or highest absolute edge weight as the highest weight to scale the edge widths too. To compare several graphs, set this argument to a higher value than any edge weight in the graphs (typically 1 for correlations).}
\item{cut}{In weighted graphs, this argument can be used to cut the scaling of edges in width and color saturation. Edges with absolute weights over this value will have the strongest color intensity and become wider the stronger they are, and edges with absolute weights under this value will have the smallest width and become vaguer the weaker the weight. If this is set to NULL, no cutoff is used and all edges vary in width and color. Defaults to NULL for graphs with less then 50 nodes and 0.3 to larger graphs.}

\item{details}{Logical indicating if minimum, maximum and cutoff score should be printed under the graph. Defaults to FALSE.}

\item{mar}{A vector of the form c(bottom, left, top, right) which gives the margins. Works similar to the argument in par(). Defaults to c(3,3,3,3)}
  \item{filetype}{A character containing the file type to save the output in. "R" outputs in a new R window, "pdf" creates a pdf file. "svg" creates a svg file (requires RSVGTipsDevice). "tex" creates LaTeX code for the graph (requires tikzDevice). 'jpg', 'tiff' and 'png' can also be used. If this is given any other string (e.g. filetype="") no device is opened. Defaults to 'R' if the current device is the NULL-device or no new device if there already is an open device. A function such as \code{x11} can also be used}
	\item{filename}{Name of the file without extension}
	\item{width}{Width of the plot, in inches}
	\item{height}{Height of the plot, in inches}
  \item{normalize}{Logical, should the plot be normalized to the plot size. If TRUE (default) border width, vertex size, edge width and arrow sizes are adjusted to look the same for all sizes of the plot, corresponding to what they would look in a 7 by 7 inches plot if normalize is FALSE.}
  
      \item{DoNotPlot}{Runs qgraph but does not plot. Useful for saving the output (i.e. layout) without plotting}
  \item{plot}{Logical. Should a new plot be made? Defaults to TRUE. Set to FALSE to add the graph to the existing plot.}
	\item{rescale}{Logical. Defines if the layout should be rescaled to fit the -1 to 1 x and y area. Defaults to TRUE. Can best be used in combination with plot=FALSE.}
  
    \item{label.cex}{Scalar on the label size.}
    \item{label.color}{Character containing the color of the labels, defaults to "black"}
    \item{borders}{Logical indicating if borders should be plotted, defaults to TRUE.}
  \item{border.color}{Color vector indicating colors of the borders. Is repeated if length is equal to 1. Defaults to "black"}
    \item{border.width}{Controls the width of the border. Defaults to 2 and is comparable to 'lwd' argument in 'points'.}
    
        \item{polygonList}{ A list contaning named lists for each element to include polygons to lookup in the \code{shape} arguments. Each element must be named as they are used in \code{shape} and contain a list with elements \code{x} and \code{y} contaning the coordinates of the polygon. By default \code{ellipse} and \code{heart} are added to this list. These polygons are scaled according to \code{vsize} and \code{vsize2}}
  
        \item{vTrans}{Transparency of the nodes, must be an integer between 0 and 255, 255 indicating no transparency. Defaults to 255}
        
          \item{label.prop}{Controls the proportion of the width of the node that the label rescales to. Defaults to 0.
  9.}
    \item{label.norm}{A single string that is used to normalize label size. If the width of the label is lower than the width of the hypothetical label given by this argument the width of label given by this argument is used instead. Defaults to "OOO" so that every label up to three characters has the same fontsize.}
    \item{label.scale}{Logical indicating if labels should be scaled to fit the node. Defaults to TRUE.}
\item{label.font}{Integer specifying the label font of nodes. Can be a vector with value for each node}


   \item{posCol}{Color of positive edges. Can be a vector of two to indicate color of edges under 'cut' value and color of edges over 'cut' value. If 'fade is set to TRUE the first color will be faded the weaker the edge weight is. If this is only one element this color will also be used for edges stronger than the 'cut' value. Defaults to c("#009900","darkgreen")}
\item{negCol}{Color of negative edges. Can be a vector of two to indicate color of edges under 'cut' value and color of edges over 'cut' value. If 'fade is set to TRUE the first color will be faded the weaker the edge weight is. If this is only one element this color will also be used for edges stronger than the 'cut' value. Defaults to c("#BF0000","red")}
\item{unCol}{Color to indicate the default edge color of unweighted graphs. Defaults to "#808080".}
\item{colFactor}{Exponent of transformation in color intensity of relative strength. Defaults to 1 for linear behavior.}
  \item{trans}{In weighted graphs: logical indicating if the edges should fade to white (FALSE) or become more transparent (TRUE; use this only if you use a background). In directed graphs this is a value between 0 and 1 indicating the level of transparency. (also used as 'transparency')}
    \item{fade}{if TRUE (default) and if 'edge.color' is assigned, transparency will be added to edges that are not transparent (or for which no transparency has been assigned) relative to the edge strength, similar if 'trans' is set to TRUE.}
  \item{loop}{This can be used to scale the size of the loop. defaults to 1.}
        
         \item{curvePivot}{Quantile to pivot curves on. This can be used to, rather than round edges, make straight edges as curves with "knicks" in them. Can be logical or numeric. \code{FALSE} (default) indicates no pivoting in the curved edges, a number indicates the quantile (and one minus this value as quantile) on which to pivot curved edges and \code{TRUE} indicates a value of 0.1.}
        \item{curvePivotShape}{The shape of the curve around the pivots, as used in \code{xspline}. Defaults to \code{0.25}.}
        
            \item{edge.label.bg}{Either a logical or character vector/matrix. Indicates the background behind edge labels. If TRUE (default) a white background is plotted behind each edge label. If FALSE no background is plotted behind edge labels. Can also be a single color character, a vector or matrix of color vectors for each edge.}
    \item{edge.label.position}{Vetor of numbers between 0 and 1 controlling the relative position of each edge label. Defaults to 0.5 for placing edge labels at the middle of the edge.}
\item{edge.label.font}{Integer specifying the label font of edges. Can be a vector or matrix with value for each node}

  \item{layout.par}{A list of arguments passed to \code{\link{qgraph.layout.fruchtermanreingold}} when layout="spring" or to an igraph function when such a function is assigned to 'layout'}
  
    \item{bg}{If this is TRUE, a background is plotted in which node colors cast a light of that color on a black background. Can also be a character containing the color of the background Defaults to FALSE}
  \item{bgcontrol}{The higher this is, the less light each node gives if bg=TRUE. Defaults to 6.}
	\item{bgres}{square root of the number of pixels used in bg=TRUE, defaults to 100.}
  
    \item{pty}{See 'par'}
    
    \item{font}{Integer specifying the default font for node and edge labels}
    
    
      \item{arrows}{A logical indicating if arrows should be drawn, or a number indicating how much arrows should be drawn on each edge. If this is TRUE, a simple arrow is plotted, if this is a number, arrows are put in the middle of the edges.}
  \item{arrowAngle}{Angle of the arrowhead, in radians. Defaults to pi/8 for unweighted graphs and pi/4 for weighted graphs.}
	\item{asize}{Size of the arrowhead. Defaults to 2 for graphs with more than 10 nodes and 2 to smaller graphs.}
	\item{open}{Logical indicating if open (TRUE) or closed (FALSE) arrowheads should be drawn.}
  
    \item{weighted}{Logical that can be used to force either a weighted graph (TRUE) or an unweighted graph(FALSE).}
  
    \item{XKCD}{If set to TRUE the graph is plotted in XKCD style based on http://stackoverflow.com/a/12680841/567015.}
    
  }
}
}


\details{
The default \code{"tree"} layout under default rotation places the nodes in one of four horizontal levels. At the top the exogenous manifest variables, under that the exogenous latent variables, under that the endogenous latent variables and at the bottom the endogenous manifest variables. If one of these kinds of variables does not exist its level is omitted. Afterwards, the \code{rotation} argument will rotate the graph and the \code{"circle"} layout will make the layout circular using these levels as nested circles.

If not manually set (see \code{\link{semPlotModel-edit}}), \code{semPath} will automatically try to set the endogenous and exogenous variables, such that the resulting layout looks good. A latent variable is identified as \emph{exogenous} if it is not on the right hand side of a directed edge (\code{->} or \code{~>}) with another latent variable as node of origin. A manifest variable is set as \emph{exogenous} if it is only connected, in any way, to exogenous latent variables and if it is not the right hand side (dependent variable) of a regression edge (\code{~>}). If all variables are set to exogenous this way, they are all set to endogenous for consistency in the layouts. Afterwards, manifest variables only used in formative measurement models (only outgoing directed edges to latents) are set to exogenous again so that MIMIC models are displayed properly.

Intercepts are placed on the same level as the variable, either on the left or right side of the node (pointing outward from the center). Residuals for manifest variables are placed at the top or bottom (for exogenous and endogenous manifests respectively). Residuals of latents are placed at the bottom or top respectively for exogenous and endogenous variables, but is switched if the latent is not connected to a manifest. Residuals for the leftmost and rightmost latent are placed at the left and right side respectively, or diagonal if the latent is connected to an intercept.

The \code{"tree2"} and \code{"circle2"} layouts call the \code{\link[igraph]{layout.reingold.tilford}} function from the \code{igraph} package. As roots are used the first available variables of the following list:
\itemize{
\item Intercepts of exogenous manifests
\item Exogenous manifest
\item Intercepts of exogenous latents
\item Exogenous latents
\item Interceots of endogenous latents
\item Endogenous latents
\item Intercepts of endogenous manifests
\item The endogenous manifest with the most outgoing edges (this should not be possible by default, but can be manually set)
\item The most connected endogenous manigest.
}
To compute an optimal layout \code{\link[igraph]{layout.reingold.tilford}} is run on a slightly altered version of the path diagram. In this version, the direction of edges from all intercepts that are not roots is reversed, the direction of all edges leading to exogenous manifests is reversed and all bidirectional edges are removed.
}

\value{
A \code{"qgraph"} object as returned by \code{\link[qgraph]{qgraph}}. This object can be used to alter the graph (such as manually redefining the layout) and to plot the graph again with different arguments.

If there are multiple groups a list is returned with a "qgraph" object for each path diagram that has been produced.
}
\references{

Fruchterman, T. & Reingold, E. (1991). Graph drawing by force-directed placement. Software - Pract. Exp. 21, 1129-1164.

Reingold, E and Tilford, J (1981). Tidier drawing of trees. IEEE Trans. on Softw. Eng., SE-7(2):223-228.

 Csardi G, Nepusz T (2006). The igraph software package for complex network research, InterJournal, Complex Systems 1695. http://igraph.sf.net
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
 \code{\link[qgraph]{qgraph}}
 \code{\link{semPlotModel}}
 \code{\link{semPlotModel-class}}
 \code{\link{semCors}}
 \code{\link{lisrelModel}}
 \code{\link{semSyntax}}
}

\examples{ # Regression analysis with interaction effects ----------------------------

# A silly dataset:
X <- rnorm(100)
Y <- rnorm(100)
Z <- rnorm(1)*X + rnorm(1)*Y + rnorm(1)*X*Y
DF <- data.frame(X,Y,Z)

# Regression including interaction:
res <- lm(Z ~ X*Y, data = DF)

# Path diagram:
semPaths(res, intAtSide=TRUE)

# Standardized estimates:
semPaths(res,"std","hide", intAtSide=TRUE)

# Simple CFA -------------------------------------------
library("lavaan")
example(cfa)

semPaths(fit, 'std', 'est', curvePivot = FALSE, curveAdjacent = TRUE, style = "lisrel")


# MIMIC model ----------------------------------------------------
## Lavaan
\dontrun{
library("lavaan")

# Example 5.8 from mplus user guide:
Data <- read.table("http://www.statmodel.com/usersguide/chap5/ex5.8.dat")
names(Data) <- c(paste("y", 1:6, sep=""),
                 paste("x", 1:3, sep=""))

# Model:
model.Lavaan <- 'f1 =~ y1 + y2 + y3
f2 =~ y4 + y5 + y6
f1 + f2 ~ x1 + x2 + x3 '

# Run Lavaan:
library("lavaan")
fit <- lavaan:::cfa(model.Lavaan, data=Data, std.lv=TRUE)

# Plot path diagram:
semPaths(fit,title=FALSE, curvePivot = TRUE)

# Omit exogenous covariances:
semPaths(fit,title=FALSE, curvePivot = TRUE, exoVar = FALSE, exoCov = FALSE)

# Standardized parameters:
semPaths(fit,"std",edge.label.cex=0.5, curvePivot = TRUE, exoVar = FALSE, 
  exoCov = FALSE)

}
## Mplus

# Same model, now using mplus output:
outfile <- tempfile(fileext=".out")
download.file("http://www.statmodel.com/usersguide/chap5/ex5.8.out",outfile)

# Plot model:
semPaths(outfile,intercepts=FALSE)
# Note that mplus did not report the fixed variances of the exogenous variables.




# Thresholds  -----------------------------------------------------
\dontrun{
## Lavaan

# Example 5.8 from mplus user guide:
Data <- read.table("http://www.statmodel.com/usersguide/chap5/ex5.2.dat")
names(Data) <- c("u1","u2","u3","u4","u5","u6")
Data <- as.data.frame(lapply(Data, ordered))

# Lavaan model:
model <- ' f1 =~ u1 + u2 + u3; f2 =~ u4 + u5 + u6 '

# Run Lavaan:
fit <- lavaan::cfa(model, data=Data)

# Plot path diagram:
semPaths(fit,intercepts=FALSE)
}

## Mplus

# Same model, now using mplus output:
outfile <- tempfile(fileext=".out")
download.file("http://www.statmodel.com/usersguide/chap5/ex5.2.out",outfile)

# Plot model:
semPaths(outfile)




# OpenMx ----------------------------------------------------------
\dontrun{
# To isntall OpenMx see:
# http://openmx.psyc.virginia.edu/

library("OpenMx")

# Example from mxRun help page:
# Create and run the 1-factor CFA on the openmx.psyc.virginia.edu front page
data(demoOneFactor)  # load the demoOneFactor dataframe
manifests <- names(demoOneFactor) # set the manifest to the 5 demo variables
latents   <- c("G")  # define 1 latent variable
model <- mxModel("One Factor", type="RAM",
  manifestVars = manifests,
	latentVars   = latents,
	mxPath(from=latents  , to=manifests),
	mxPath(from=manifests, arrows=2),
	mxPath(from=latents  , arrows=2, free=FALSE, values=1.0),
	mxData(cov(demoOneFactor), type="cov", numObs=500)
)
model <- mxRun(model) #run model, returning the result

# Plot with colors from OpenMx front page:
semPaths(model, color = list(
              lat = rgb(245, 253, 118, maxColorValue = 255), 
              man = rgb(155, 253, 175, maxColorValue = 255)),
          mar = c(10, 5, 10, 5))

## Factor Analysis:
source("http://openmx.psyc.virginia.edu/svn/trunk/demo/TwoFactorModel_PathCov.R")
semPaths(twoFactorFit, layout = "tree2")





# Multi-group analysis -------------------------------------------
## LISREL:
# Download measurment invariance example:
modFile <- tempfile(fileext=".OUT")
download.file("http://sachaepskamp.com/files/mi1.OUT",modFile)
layout(t(1:2))
semPaths(modFile,"eq",ask=FALSE, intAtSide = TRUE, mar = c(8, 1, 5, 1))
# Color indicates equality constraints.
}
}
