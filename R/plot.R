

#' Compute font size of col/rownames in complexHeatmap
#'
#' @param n Single numeric value. Number of column/rows
#' @param ... Other parameters passed to gpar
#'
#' @return a gpar object
#' @export
#'
#' @examples
#' autoGparFontSizeMatrix(5)
autoGparFontSizeMatrix<-function(n,...){
	require("grid")
	n=max(n,50)
	n=min(n,1000)
	return(gpar(fontsize=1/n*600,...))
}

#' Multiple ggplot one one page
#'
#'
#' @param ... Several plot from ggplot. Each given as an argument.
#' @param plotlist Plots given as a list. Each element contains a plot from grid or inherited package.
#' @param cols Number of columns. If layout is NULL, then use 'cols' to determine layout
#' @param layout A 2D matrix of numeric value (from one the the number of plot) indicating the layout to be plotted. Matrix()
#'
#' @return NULL
#' @export
#'
#' @examples
#' require(ggplot2)
#' require(ComplexHeatmap)
#'
#' p1<-ggplot(data = data.frame(x=1:5,y=1:5),aes(x=x,y=y))+geom_point()+ggtitle("p1")
#' p2<-qplot(c("a","a","b","b","b","c"))+ggtitle("p2")
#' p3<-qplot(c("a","a","b","b","b"))+ggtitle("p3")
#' p4<-qplot(c("a","b","c","d"))+ggtitle("p4")
#'
#' multiplot(p1,p2,p3,p4)
#' plotList<-list(p1,p2,p3,p4)
#' multiplot(plotlist = plotList)
#' multiplot(plotlist = plotList,cols = 2)
#'
#' layout<-matrix(data = c(
#' 	4,3,
#' 	2,1
#' ),ncol = 2,byrow = TRUE)
#' multiplot(plotlist = plotList,layout = layout)
#' layout<-matrix(data = c(
#' 	1,2,3,
#' 	4,0,0
#' ),ncol = 3,byrow = TRUE)
#' multiplot(plotlist = plotList,layout = layout)

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
	require(grid)

	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)

	numPlots = length(plots)

	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
										 ncol = cols, nrow = ceiling(numPlots/cols))
	}
	if (numPlots==1) {
		print(plots[[1]])

	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
																			layout.pos.col = matchidx$col))
		}
	}
}


#' Default Color scale of ggplot2
#'
#' @param n Single integer. Number of colors to be returned.
#' @param h numeric vector of 2 value. Hue range (see `hcl`)
#'
#' @return A character vector of colors.
#' @export
#'
#' @examples
#' ggplotColours()
#' ggplotColours(5)
#' ggplotColours(5,)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
	if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
	hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


#' Convert color from additive to subtracting mixing
#'
#' @param color A vector of 3 numeric containing rgb values (from 0 to 1), or a single character of the color name or hex code
#'
#' @return a vector containing rgb values if returnHex=FALSE, otherwise a color as a character format.
#' @export
#'
#' @examples
#' convertColorAdd2Sub(c(0,0,0))
#' convertColorAdd2Sub(c(0,0,0),returnHex = FALSE)
#'
#' color<-rgb(red=.1,green=.1,blue=.1)
#' invertedColor<-convertColorAdd2Sub(color)
#'
#' plotPalette(c(color,invertedColor))
convertColorAdd2Sub<-function(color,returnHex=TRUE){
	if(is.character(color)){
		color<-col2rgb(color)[,1]/255
	}
	if(length(color)<3) stop("if blue or green is null, the color argument must contain 3 numeric values.")
	red<-color[1]
	green<-color[2]
	blue<-color[3]

	newRed=mean(c(1-green,1-blue))
	newGreen=mean(c(1-red,1-blue))
	newBlue=mean(c(1-red,1-green))
	if(returnHex){
		return(rgb(newRed,newGreen,newBlue))
	}else{
		return(c("red"=newRed,"green"=newGreen,"blue"=newBlue))
	}
}


#' Best theoretical color palette (wrapper for qualpal)
#'
#' @param n The number of colors to generate.
#' @param colorspace A color space to generate colors from. See qualpalr::qualpal
#' @param cvd Color vision deficiency adaptation. Use cvd_severity to set the severity of color vision deficiency to adapt to. Permissible values are "protan", "deutan", and "tritan".
#' @param cvd_severity Severity of color vision deficiency to adapt to. Can take any value from 0, for normal vision (the default), and 1, for dichromatic vision.
#' @param n_threads The number of threads to use, provided to setThreadOptions if non-null.
#'
#' @return Colors in hex format.
#' @export
#'
#' @examples
#' mostDistantColor(3)
#'
mostDistantColor<-function(n, colorspace = "rainbow", cvd = c("protan", "deutan","tritan"), cvd_severity = 0, n_threads = NULL){
	require("qualpalr")
	if(n==1) return("#000000")
	qualpal(n=n, colorspace = colorspace, cvd = cvd, cvd_severity = cvd_severity, n_threads = n_threads)$hex
}


#' Compute density value for the point of a a 2D-space
#'
#' @param mat numeric matrix of point coordinates. Each row is a point, 1st column = x, 2nd column=y.
#' @param eps radius of the eps-neighborhood, i.e., bandwidth of the uniform kernel).
#'
#' @return A vector of numeric value of the length equal to the number of points representing the density
#' @export
#'
#' @examples
#' coord<-matrix(rnorm(2000),ncol = 2)
#' pointdensity.nrd(coord)
pointdensity.nrd<-function(mat,eps=1){
	if(!is.matrix(mat)) mat<-as.matrix(mat)
	require(dbscan)
	require(MASS)
	pointdensity(apply(mat,2,function(d) d/bandwidth.nrd(d)),eps = eps,type = "density")
}


#' log10(x+1) continuous scale for ggplot2
#'
#' @return A transformation object.
#' @export
#'
#' @examples
#' ggplot(data.frame(x=c(0,10,100,1000),y=1:4),mapping = aes(x=x,y=y))+geom_point()+scale_x_continuous(trans = log10plus1())
log10plus1<-function(){
	require(scales)
	scales::trans_new(name = "log10plus1",transform = function(x) log10(x+1),inverse = function(x) 10^x-1 ,domain=c(0,Inf))
}

#' Plot colors
#'
#' @param colorScale Character vector, Colors to plot.
#' @param continuousStep Number of plotted color with intermediate between the given colors.
#'
#' @return NULL, draw in current graphic device.
#' @export
#'
#' @examples
#' plotPalette(c("#001D6E","white","#E6001F"))
#' plotPalette(c("#001D6E","white","#E6001F"),continuousStep = 50)
plotPalette<-function(colorScale,continuousStep=NULL){
	if(is.null(continuousStep)){
		image(1:length(colorScale), 1, as.matrix(1:length(colorScale)),
					col=colorScale,
					xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
		if(!is.null(names(colorScale))) axis(1,1:length(colorScale),names(colorScale))
	}else{
		require(circlize)
		br = round(seq(1,continuousStep,length.out = length(colorScale)))
		cols = colorRamp2(breaks = br,colors = colorScale)(1:continuousStep)
		image(1:continuousStep, 1, as.matrix(1:continuousStep),
					col=cols,
					xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
	}
}


#' Compute a color scale function from numeric values by interpolation
#'
#' @param colors A character vector containing the colors.
#' @param values A numeric vector of the value to has to be mapped to colors.
#' @param useProb Logical. Use quantile probability to map the colors. Else the min and max of values will be mapped to first and last color and interpolated continuously.
#' @param probs A numeric vector (between 0 and 1) same length as color or NULL. Quantile probability of the values that will be mapped to colors.
#' @param minProb A numeric value (between 0 and 1). If `useProb=TRUE` and `probs=NULL` this will be the quantile of the value for the first color, quantile will be mapped continuously as to the maxProb.
#' @param maxProb A numeric value (between 0 and 1).
#' @param midColorIs0 Logical. Force that 0 return the midColor.
#' @param returnColorFun Logical.Return converted values to colors or the scale function.
#' @param returnGGscale Logical. Return a ggplot2 gradiantn scale.
#' @param geomAes "fill" or "color". Ggplot layer that will receive the scale.
#' @param geomArgument list of additional argument to pass to the ggplot2 gradiantn scale.
#'
#' @return A vector of colors, or a function if `returnColorFun=TRUE` or a ggplot scale if `returnGGscale=TRUE`.
#' @export
#'
#' @examples
#'
#' values=sort(rnorm(100))
#'
#' plotPalette(computeColorScaleFun(colors = c("black","red"),values = values,returnColorFun = FALSE))
#' plotPalette(computeColorScaleFun(colors = c("blue","white","red"),values = values,returnColorFun = FALSE))
#' plotPalette(computeColorScaleFun(colors = c("blue","white","red"),values = values,returnColorFun = FALSE,midColorIs0 = TRUE))
#' plotPalette(computeColorScaleFun(colors = c("blue","white","red"),values = values,returnColorFun = FALSE,useProb = TRUE))
#' plotPalette(computeColorScaleFun(colors = c("blue","white","red"),values = values,returnColorFun = FALSE,useProb = TRUE,probs = c(.25,.5,.75)))
#'
#' colorFun<-computeColorScaleFun(colors = c("blue","white","red"),values = values,returnColorFun = TRUE,useProb = TRUE)
#' plotPalette(c(colorFun(-1),colorFun(0),colorFun(1)))
#'
#' dat<-data.frame(x=rnorm(10),y=rnorm(10),expr=rnorm(10))
#' ggplot(dat,aes(x=x,y=y,fill=expr))+
#' 	geom_point(size=5,shape=21)+theme_bw()+
#' 	computeColorScaleFun(colors = c("blue","white","red"),values = dat$expr,returnGGscale = TRUE,useProb = TRUE,geomAes = "fill")
computeColorScaleFun<-function(colors,values,useProb=FALSE,probs = NULL,minProb=0.05,maxProb=0.95,
															 midColorIs0=FALSE,returnColorFun=TRUE,returnGGscale=FALSE,geomAes="fill",geomArgument=list()){
	require(circlize)
	if(!useProb){
		breaks = seq(min(values),max(values),length.out = length(colors))
	}else{
		if(is.null(probs)){
			probs=seq(minProb,maxProb,length.out = length(colors))
		}
		breaks<-quantile(values,probs=probs)
	}
	if(midColorIs0 & (length(colors) %% 2 == 1)){
		breaks[ceiling(length(breaks)/2)]<-0
	}
	colorFun<-colorRamp2(breaks=breaks,colors=colors)
	if(returnGGscale){
		require(ggplot2)
		scaledBreaks<-linearScale(values,c(0,1),returnFunction = TRUE)(breaks)
		if(scaledBreaks[1]>0){
			scaledBreaks<-c(0,scaledBreaks)
			colors<-c(colors[1],colors)
		}
		if(scaledBreaks[length(scaledBreaks)]<1){
			scaledBreaks<-c(scaledBreaks,1)
			colors<-c(colors,colors[length(colors)])
		}

		geomArgument$values <- scaledBreaks
		geomArgument$colors <- colors
		return(do.call(paste0("scale_",geomAes,"_gradientn"),geomArgument))
	}
	if(returnColorFun){
		return(colorFun)
	}else{
		return(colorFun(values))
	}
}



#' Generate a list of value/color mapping
#'
#' @param annots Dataframe. Can contain factors or numeric. Contains the values that has to be mapped.
#' @param colorScales
#' List or NULL. Precomputed color scales. Color scales will be only generated for the features not described.
#' Must be in the format of a list named by columns of `annots`.
#' Each element contains the colors at breaks for continuous values or a mapping function if `returnContinuousFun=TRUE` (a function that return a color for a given numeric value).
#' In the case of factors, the colors named to its corresponding level.
#' @param discreteFuns A list functions that take a single integer n and return n colors.
#' If several functions are provided it will be used for each factor column successively.
#' @param returnContinuousFun Logical. Return a mapping function for continuous values instead of a vector of colors.
#' @param continuousPalettes A list of color vector. If several vector are provided it will be used for each numerical column successively.
#' @param ... Parameters passed to `computeColorScaleFun`.
#'
#' @return A list describing the color scale of each column of `annots`, in the same format than the argument `colorScales`
#' @export
#'
#' @examples
#' data("iris")
#'
#' genColorsForAnnots(iris)
#'
#' precomputedColorScale<-list(Species=c("setosa"="red","versicolor"="blue","virginica"="grey"))
#'
#' genColorsForAnnots(iris,colorScales = precomputedColorScale)
#'
#' colorScales<-genColorsForAnnots(iris,returnContinuousFun = TRUE)
#' colorScales$Sepal.Length(4.5)
#' colorScales$Species
#'
#' library(ComplexHeatmap)
#' Heatmap(rowScale(t(iris[,c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width")]),center = TRUE,scaled = TRUE),
#' 				top_annotation = genTopAnnot(iris["Species"],colorScales=colorScales["Species"]))


genColorsForAnnots<-function(annots, colorScales=NULL,discreteFuns = list(mostDistantColor2),returnContinuousFun=FALSE ,
														 continuousPalettes=list(c("#440154","#6BA75B","#FDE725"),c("#2EB538","#1D1D1B","#DC0900"),c("#FFFFC6","#FF821B","#950961")),...){

	if(is.null(colnames(annots))) stop("annots must have colnames.")
	annotNames <- colorScalesToGen <- colnames(annots)
	newColorScales<-list()

	if(!is.null(colorScales)){
		for(colorScaleName in names(colorScales)){
			if(!colorScaleName%in%colorScalesToGen) stop("Condition '",colorScaleName,"' does not match with existing condition names")
			colorScale<-colorScales[[colorScaleName]]
			annotVect<-annots[,colorScaleName]
			if(!is.null(names(colorScale))){ #factors
				if(is.numeric(annotVect)){
					warning(colorScaleName," is numeric but encoded as factors (color vector has names). It will be converted to factors.")
					annots[,colorScaleName]<-as.factor(as.character(sampleAnnot[,colorScaleName]))
					annotVect<-annots[,colorScaleName]
				}else if(!is.factor(annotVect)){
					stop(colorScaleName," is not factors or numeric, please check the sample annotation table.")
				}
				if(sum(!levels(annotVect) %in% names(colorScale))>0) stop("Levels of ",colorScaleName," are existing in sample annotation table but not in ",colorScaleFile)
				newColorScales[[colorScaleName]]<-colorScale[levels(annotVect)]
			}else{ #numeric
				if(!is.numeric(annotVect)) stop(colorScaleName," is not numeric but encoded as numeric (color vector has no names)")
				if(is.function(colorScale) & !returnContinuousFun) stop("You must not provide function in colorScales if returnContinuousFun=FALSE")
				if(!is.function(colorScale) & returnContinuousFun){
					newColorScales[[colorScaleName]]<-computeColorScaleFun(colorScale,values = annotVect,returnColorFun = TRUE,...)
				}else{
					newColorScales[[colorScaleName]]<-colorScale
				}

			}
		}
		colorScalesToGen<-setdiff(colorScalesToGen,names(newColorScales))
	}
	cN<-1
	cF<-1
	for(colorScaleName in colorScalesToGen){
		annotVect<-annots[,colorScaleName]
		if(is.numeric(annotVect)){
			if(returnContinuousFun){
				newColorScales[[colorScaleName]]<-computeColorScaleFun(continuousPalettes[[cN]],values = annotVect,returnColorFun = TRUE,...)
			}else{
				newColorScales[[colorScaleName]]<-continuousPalettes[[cN]]
			}
			cN<-cN+1
			if(cN>length(continuousPalettes)) cN<-1
		}else{
			annots[,colorScaleName]<-as.factor(as.character(annots[,colorScaleName]))
			annotVect<-annots[,colorScaleName]
			newColorScales[[colorScaleName]]<-discreteFuns[[cF]](nlevels(annotVect))
			names(newColorScales[[colorScaleName]])<-levels(annotVect)
			cF<-cF+1
			if(cF>length(discreteFuns)) cF<-1
		}
	}
	newColorScales
}

#' Add line break between factors and remove line in the middle if x axis is discrete.
#'
#' @param gg ggplot object
#' @param borderSize Single numeric value. Size width of the line.
#' @param borderColor Single character value. Color of the line.
#'
#' @return A ggplot objet.
#' @export
#'
#' @examples
#' g<-ggplot(data.frame(x=c("A","A","B","B","B","C")),aes(x=x))+geom_bar()
#' g
#' ggBorderedFactors(g)
#' ggBorderedFactors(g,borderColor="white",borderSize=1.5)
ggBorderedFactors<-function(gg,borderSize=.75,borderColor="grey75"){
	require(ggplot2)
	nX<-nlevels(as.factor(gg$data[,quo_name(gg$mapping$x)]))
	gg+geom_vline(xintercept = seq(1.5,nX -0.5, 1),size=borderSize,color=borderColor)+
		scale_x_discrete(expand = c(0,0.5, 0, 0.5))+
		theme(
			panel.grid.major.x = element_line(colour = NA),
			panel.grid.minor.x = element_line(colour = NA),
		)
}


#' Plot a PCA in 2D
#'
#' @param pca PCA object created by `PCA` function.
#' @param comp A vector of 2 integer. Principal Component to be plotted.
#' @param colorBy A vector of factor or numeric same size as number of samples in PCA. Color each point on the projection by its corresponding value.
#' @param plotVars Logical. Plot variable loadings instead of coordinates of samples on PCs.
#' @param pointSize Single numeric. Size of points.
#' @param plotText Logical. Plot sample/variable name instead of point.
#' @param ratioPerPercentVar Logical. Ratio of the plot is corresponding to the percentage of explained variation from each PC.
#' @param main Single character. Main title of the graph.
#' @param ellipse Logical. Add stat_ellipse to the projection (eventually several ellipses if colorBy is a factor).
#' @param colorScale A vector of colors. Equal to the desired number of breaks for continuous values or number of levels for factors.
#' @param returnGraph Logical. Return the graph as a ggplot object instead of printing it.
#' @param outlierLabel Logical. Add a `ggrepel` label for point in low-density area of the plot.
#' @param outlierLabelThres Numeric. Density value threshold below that the labelled is displayed if `outlierLabel=TRUE`.
#' @param outlierLabelSize Numeric. Size of outlier labels.
#' @param customRatio Single numeric. Custom ratio for the plot.
#' @param ... Arguments passed to `proj2d`
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' data(iris)
#' pca<-PCA(iris[,1:4],transpose = FALSE,scale = TRUE,center = TRUE)
#' pca2d(pca)
#' pca2d(pca,comp = c(1,3))
#'
#' pca2d(pca,colorBy = iris$Species,ratioPerPercentVar = TRUE)
#' pca2d(pca,colorBy = iris$Species,colorScale=c("blue","white","red"))
#'
#' pca2d(pca,plotVars = TRUE,plotText = TRUE,pointSize = 5)
#'
#' pca2d(pca,outlierLabel = TRUE,colorBy = iris$Species)
#' pca2d(pca,outlierLabel = TRUE,colorBy = iris$Species,ellipse = TRUE,outlierLabelThres = .2)
#'
#' pca2d(pca,outlierLabel = TRUE,colorBy = iris$Species,ellipse = TRUE,outlierLabelThres = .2,returnGraph = TRUE)+
#' 	theme_dark()
pca2d<-function(pca, comp=1:2,colorBy=NULL, plotVars = FALSE, pointSize=2, plotText=FALSE,ratioPerPercentVar=FALSE,main=NULL,
								ellipse=FALSE,colorScale=NULL,returnGraph=FALSE, outlierLabel=FALSE,outlierLabelThres=NULL,
								outlierLabelSize=3,customRatio=NULL,...){
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(comp)!=2) stop("You must give a vector of 2 integer for comp parameter");
	percentVar <- pca$percentVar
	if(ratioPerPercentVar){
		customRatio<-percentVar[comp[2]]/percentVar[comp[1]]
	}
	dt<-pca[[ifelse(plotVars, "rotation","x")]]
	graph<-proj2d(coord = dt,axis = comp,colorBy = colorBy,returnGraph = T,colorScale = colorScale,main = main,plotText = plotText,
								ellipse = ellipse,pointSize=pointSize,customRatio=customRatio,...)+
		xlab(paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance")) +
		ylab(paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance"))
	if(outlierLabel){
		require(ggrepel)
		if(is.null(rownames(dt))) rownames(dt)<-as.character(1:nrow(dt))
		dtOutlier<-data.frame(x=dt[,comp[1]],y=dt[,comp[2]],label=rownames(dt))
		dens<-pointdensity.nrd(dtOutlier[,1:2]);
		if(is.null(outlierLabelThres)) outlierLabelThres=Mode(dens)/3
		dtOutlier<-dtOutlier[dens/max(dens)<outlierLabelThres,]
		graph<-graph+geom_text_repel(data = dtOutlier,mapping = aes(x=x,y=y,label=label),
																 size=outlierLabelSize,fontface="italic",inherit.aes = FALSE)
	}

	if(returnGraph){
		return(graph)
	}else{
		print(graph)
	}
}

#' Plot a 2D projection.
#'
#' @param coord Matrix or dataframe containing at least two vectors of numeric corresponding to the x/y coordinates.
#' @param colorBy A vector of factor or numeric same size as number of samples in PCA. Color each point on the projection by its corresponding value.
#' @param axis  A vector of 2 integer. The column index of `coord` containing respectively the x and y coordinates.
#' @param pointSize  Single numeric. Size of points.
#' @param plotText Logical. Plot sample instead of point.
#' @param main Single character. Main title of the graph.
#' @param alpha Single numeric (from 0 to 1). Opacity of the points.
#' @param ellipse Logical. Add stat_ellipse to the projection (eventually several ellipses if colorBy is a factor).
#' @param emph Single character or `NULL`. Must be a level of `colorBy` if it is a factor. Emphasize a particular value on the plot.
#' @param colorScale A vector of colors. Equal to the desired number of breaks for continuous values or number of levels for factors.
#' @param returnGraph Logical. Return the graph as a ggplot object instead of printing it.
#' @param legendTitle Character. The legend title on the plot.
#' @param axis.names Vector of two characters. Axis names (x, y) showed on the plot.
#' @param na.color Character. Color of NA values if `ColorBy` is provided.
#' @param na.bg Logical. Put the NA points behind the others in term of layers.
#' @param plotFactorsCentroids Logical. If `ColorBy` is a factor Display the name of the factor at the centroid of its points.
#' @param pointStroke Numeric. Width of points border.
#' @param strokeColor Character. Color of points border.
#' @param funAutoColorScale A function that return n colors with n as a first argument. Responsible for color mapping of the levels.
#' @param fixedCoord Logical. The ratio of the plot is representative of the ratio between x and y values.
#' @param plotLabelRepel Logical. Plot the rownames of the coordinates as a ggrepel label of each point.
#' @param labelSize Numeric. Size of the labels if `plotLabelRepel=TRUE`.
#' @param sizeProp2Dens Logical. Size of point is inversely proportional to the 2D density of its area.
#' @param densEps Numeric. Radius of the eps-neighborhood, i.e., bandwidth of the uniform kernel). Used if `sizeProp2Dens=TRUE`.
#' @param nnMatrix A matrix of integer, number of row equal to the number of points, as returned by `make.umap(..., ,ret_nn = TRUE)`. Plot segments between the neighbors.
#' @param nnSegmentParam List of arguments for customizing the nn segments. Passed to `geom_segment`.
#' @param useScatterMore Logical. Use `geom_scatter` for making the plot. Render the points as a pixel image, difficult to modify afterward but quicker plotting for very large datasets.
#' @param customRatio Single numeric. Custom ratio for the plot.
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' data(iris)
#'
#' proj2d(iris)
#' proj2d(iris,axis = c(1,3))
#'
#' proj2d(iris,pointSize = .1)
#' proj2d(iris,axis.names = c("Width of sepal","Length of sepal") ,main = "Two features from iris")
#'
#'
#' proj2d(iris,colorBy = iris$Species,ellipse = TRUE,plotFactorsCentroids = TRUE,legendTitle = "Species")
#' proj2d(iris,colorBy = iris$Species,pointStroke = 3,strokeColor = "grey")
#'
#' proj2d(iris,colorBy = iris$Species,emph = "setosa",na.color = "white")
#' proj2d(iris,colorBy = iris$Species,returnGraph = TRUE)+
#' 	geom_vline(xintercept = 6)
#'
#' proj2d(iris,colorBy = iris$Species,colorScale = c("blue","white","red"))
#' proj2d(iris,colorBy = iris$Species,funAutoColorScale = mostDistantColor)
#' proj2d(iris,colorBy = iris$Petal.Length,legendTitle = "Petal.Length",fixedCoord = TRUE)
#' proj2d(iris,colorBy = iris$Petal.Length,legendTitle = "Petal.Length",colorScale = c("blue","white","red"))
#'
#' proj2d(iris,plotLabelRepel = TRUE)
#' proj2d(iris,sizeProp2Dens = TRUE)
#' proj2d(iris,sizeProp2Dens = TRUE,densEps = 2)
#'
#' umap<-make.umap(scale(iris[,c(1:4)]),ret_nn = TRUE,transpose = FALSE,n_neighbors = nrow(iris))
#'
#' proj2d(umap)
#' proj2d(umap,colorBy = iris$Species ,nnMatrix = umap$nn$euclidean$idx[,1:3],fixedCoord = TRUE)
#' proj2d(umap,colorBy = iris$Species,useScatterMore = T,fixedCoord = T,pointSize = 5)
proj2d<-function(coord,colorBy=NULL,axis=1:2, pointSize=3, plotText=FALSE,main=NULL,alpha=9/10,
								 ellipse=FALSE,emph=NULL,colorScale=NULL,returnGraph=FALSE,legendTitle="Values",axis.names=NULL,
								 na.color="grey50",na.bg=TRUE,plotFactorsCentroids=FALSE,
								 pointStroke=1/8,strokeColor="black",funAutoColorScale=ggplotColours,fixedCoord=FALSE,plotLabelRepel=FALSE,labelSize=3,
								 sizeProp2Dens=FALSE,densEps=1,nnMatrix=NULL,nnSegmentParam=list(alpha=.75,size=.1),useScatterMore=FALSE,customRatio=NULL){
	coord<-data.frame(coord)
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(axis)!=2) stop("You must give a vector of 2 integer for axis parameter");
	if((!is.null(axis.names)) & length(axis.names)!=2) stop("Error, if given axis.names parameter must contain 2 values.");
	if(!is.null(customRatio)) fixedCoord<-TRUE

	d <- data.frame(lab=rownames(coord),Axis1=coord[,axis[1]], Axis2=coord[,axis[2]],sizeMultiplier=rep(pointSize,nrow(coord)));
	if(sizeProp2Dens){
		dens<-pointdensity.nrd(d[,c("Axis1","Axis2")],eps = densEps)
		d$sizeMultiplier<-d$sizeMultiplier*(1-dens/max(dens))*2
	}
	if(is.null(colorBy)){
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab))
	}else{
		if(is.data.frame(colorBy) | is.matrix(colorBy)){
			if(!is.null(colnames(colorBy))) legendTitle<- colnames(colorBy)[1]
			colorBy<-colorBy[,1]
		}

		if(is.character(colorBy)) colorBy<-as.factor(colorBy)
		d$colorBy<-colorBy;
		colorByIsFactor<-is.factor(colorBy)
		if(!is.null(emph)){
			if(!emph%in%levels(colorBy)) stop("emph value not in levels of colorBy")
			d$colorBy<-as.character(d$colorBy)
			d$colorBy[which(d$colorBy!=emph)]<-NA
			d$colorBy<-as.factor(d$colorBy)
			d$sizeMultiplier[which(d$colorBy==emph)]<-d$sizeMultiplier[which(d$colorBy==emph)]
		}
		if(na.bg){
			indexNA<-which(is.na(d$colorBy))
			indexNotNA<-which(!is.na(d$colorBy))
			if(length(indexNA)>0){
				tempd<-d
				tempd[1:length(indexNA),]<-d[indexNA,]
				tempd[(length(indexNA)+1):nrow(d),]<-d[indexNotNA,]
				d<-tempd
			}
		}
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab,color=colorBy,fill=colorBy))
		if(is.null(colorScale)){
			colorScale<-c("grey","red")
			if(colorByIsFactor) colorScale<-funAutoColorScale(nlevels(colorBy))
		}
		funColorScaleFill<-paste0("scale_fill_",ifelse(colorByIsFactor,"manual","gradientn"))
		funColorScaleColor<-paste0("scale_color_",ifelse(colorByIsFactor,"manual","gradientn"))
		paramColorScale<-list("name"=legendTitle,na.value=na.color)
		paramColorScaleType<-ifelse(colorByIsFactor,"values","colors");paramColorScale[[paramColorScaleType]]<-colorScale
		graph<-graph+do.call(funColorScaleFill,paramColorScale)
		graph<-graph+do.call(funColorScaleColor,paramColorScale)
	}
	if(!is.null(nnMatrix)){
		if(!is.matrix(nnMatrix)) stop("If nnMatrix is not null, it should be a matrix !")
		retainCoord<-as.matrix(coord[,axis])
		segmentsCoord<-data.frame(do.call("rbind",lapply(1:nrow(retainCoord),function(i){
			t(vapply(nnMatrix[i,],FUN.VALUE = numeric(4),function(x){
				c(retainCoord[i,],retainCoord[x,])
			}))
		})));colnames(segmentsCoord)<-c("x","y","xend","yend")
		graph<-graph+do.call(
			"geom_segment",
			c(list(data = segmentsCoord,mapping = aes(x=x,y=y,xend=xend,yend=yend),inherit.aes = FALSE),nnSegmentParam)
		)
	}
	if(plotText){
		graph<-graph+geom_text(alpha=alpha,size=d$sizeMultiplier)
	}else{
		if(useScatterMore){
			require(scattermore)
			ratioX<-ratioY<-1
			if(fixedCoord){
				if(is.null(customRatio)){
					rangeA1<-range(d$Axis1)
					rangeA2<-range(d$Axis2)
					ratio<-(rangeA1[2]-rangeA1[1])/(rangeA2[2]-rangeA2[1])
				}else{
					ratio<-1/(customRatio*.7)
				}
				if(ratio>1){
					ratioX<-ratio
				}else{
					ratioY<-1/ratio
				}
			}
			if(is.null(colorBy)){
				graph<-graph+geom_scattermore(pointsize=d$sizeMultiplier,pixels=c(1000*ratioX,1000*ratioY))
			}else{
				graph<-graph+geom_scattermore(pointsize=d$sizeMultiplier,mapping = aes(color=colorBy),pixels=c(1000*ratioX,1000*ratioY))
			}
		}else{
			if(is.null(colorBy)){
				graph<-graph+geom_point(stroke=pointStroke,colour = strokeColor,shape=21,alpha=alpha,fill="black",size=d$sizeMultiplier)
			}else{
				graph<-graph+geom_point(stroke=pointStroke,colour = strokeColor,shape=21,alpha=alpha,size=d$sizeMultiplier)
			}
		}
	}
	if(plotLabelRepel){
		require(ggrepel)
		graph<-graph+geom_text_repel(color="black",size=labelSize)
	}
	if(is.null(axis.names)){
		if(is.null(colnames(coord))){
			axis.names<-paste0("Axis",axis)
		}else{
			axis.names<-colnames(coord[,c(axis[1],axis[2])])
		}
	}
	graph<-graph+xlab(axis.names[1]) + ylab(axis.names[2])
	if(!is.null(main)) graph <- graph + ggtitle(main)
	if(ellipse) graph<-graph+stat_ellipse(size=.5)
	if(fixedCoord)graph <- graph + coord_fixed(ratio = ifelse(is.null(customRatio),1,customRatio))
	if(plotFactorsCentroids){
		samplePerFactor<-lapply(levels(colorBy),function(x) which(colorBy==x))
		names(samplePerFactor)<-levels(colorBy)
		centroids<-data.frame(t(sapply(samplePerFactor,function(x) colMeans(d[x,c("Axis1","Axis2")]))),groupName=names(samplePerFactor))
		graph<-graph+geom_label(data = centroids,mapping = aes(x=Axis1, y=Axis2, label = groupName),inherit.aes = FALSE)
	}
	graph<-graph+theme(
		panel.background = element_rect(fill = NA,colour="black"),
		panel.grid.major = element_line(colour = NA),
		panel.grid.minor = element_line(colour = NA)
	) + guides(size="none")
	if(returnGraph){
		return(graph)
	}else{
		print(graph)
	}
}

#' Plot a 1D projection.
#'
#' @param variable Numeric vector. X coordinates.
#' @param colorBy A vector of factor or numeric same size as number of samples in PCA. Color each point on the projection by its corresponding value.
#' @param pointSize Single numeric. Size of points.
#' @param plotText Logical. Plot sample names instead of points.
#' @param main Single character. Main title of the graph.
#' @param alpha Single numeric (from 0 to 1). Opacity of the points.
#' @param ellipse Logical. Add stat_ellipse to the projection (eventually several ellipses if colorBy is a factor).
#' @param emph Single character or `NULL`. Must be a level of `colorBy` if it is a factor. Emphasize a particular value on the plot.
#' @param colorScale A vector of colors. Equal to the desired number of breaks for continuous values or number of levels for factors.
#' @param returnGraph Logical. Return the graph as a ggplot object instead of printing it.
#' @param legendTitle Character. The legend title on the plot.
#' @param variable.name Single characters. Axis name showed on the plot.
#' @param na.color Character. Color of NA values if `ColorBy` is provided.
#' @param na.bg Logical. Put the NA points behind the others in term of layers.
#' @param plotFactorsCentroids Logical. If `ColorBy` is a factor Display the name of the factor at the centroid of its points.
#' @param pointStroke Numeric. Width of points border.
#' @param strokeColor Character. Color of points border.
#' @param funAutoColorScale A function that return n colors with n as a first argument. Responsible for color mapping of the levels.
#' @param plotLabelRepel Logical. Plot the rownames of the coordinates as a ggrepel label of each point.
#' @param labelSize Numeric. Size of the labels if `plotLabelRepel=TRUE`.
#' @param sizeProp2Dens Logical. Size of point is inversely proportional to the density of its area.
#' @param densEps Numeric. Radius of the eps-neighborhood, i.e., bandwidth of the uniform kernel). Used if `sizeProp2Dens=TRUE`.
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' data(iris)
#'
#' proj1d(iris$Sepal.Length)
#'
#' proj1d(iris$Sepal.Length,colorBy=iris$Species)
#'
#' proj1d(iris$Sepal.Length,pointSize = .1)
#' proj1d(iris$Sepal.Length,variable.name = "Length of sepal" ,main = "Two features from iris")
#'
#'
#' proj1d(iris$Sepal.Length,colorBy = iris$Species,ellipse = TRUE,plotFactorsCentroids = TRUE,legendTitle = "Species")
#' proj1d(iris$Sepal.Length,colorBy = iris$Species,pointStroke = 3,strokeColor = "grey")
#'
#' proj1d(iris$Sepal.Length,colorBy = iris$Species,emph = "setosa",na.color = "white")
#' proj1d(iris$Sepal.Length,colorBy = iris$Species,returnGraph = TRUE)+
#' 	geom_vline(xintercept = 6)
#'
#' proj1d(iris$Sepal.Length,colorBy = iris$Species,colorScale = c("blue","white","red"))
#' proj1d(iris$Sepal.Length,colorBy = iris$Species,funAutoColorScale = mostDistantColor)
#' proj1d(iris$Sepal.Length,colorBy = iris$Petal.Length,legendTitle = "Petal.Length",colorScale = c("blue","white","red"))
#'
#' proj1d(iris$Sepal.Length,plotLabelRepel = TRUE)
#' proj1d(iris$Sepal.Length,sizeProp2Dens = TRUE)
#' proj1d(iris$Sepal.Length,sizeProp2Dens = TRUE,densEps = 2)
proj1d<-function(variable,colorBy=NULL, pointSize=3, plotText=FALSE,main=NULL,alpha=9/10,
								 ellipse=FALSE,emph=NULL,colorScale=NULL,returnGraph=FALSE,legendTitle="Values",variable.name="x",
								 na.color="grey50",na.bg=TRUE,plotFactorsCentroids=FALSE,
								 pointStroke=1/8,strokeColor="black",funAutoColorScale=ggplotColours,plotLabelRepel=FALSE,labelSize=3,
								 sizeProp2Dens=FALSE,densEps=1){
	if(!is.numeric(variable)) stop("'variable' must be a numeric vector");
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(is.null(names(variable))) names(variable)<-as.character(1:length(variable))
	d <- data.frame(lab=names(variable),Axis1=variable, Axis2=rep(0,length(variable)),sizeMultiplier=rep(pointSize,length(variable)));
	if(sizeProp2Dens){
		dens<-pointdensity.nrd(d[,c("Axis1","Axis2")],eps = densEps)
		d$sizeMultiplier<-d$sizeMultiplier*(1-dens/max(dens))*2
	}
	if(is.null(colorBy)){
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab))
	}else{
		if(is.data.frame(colorBy) | is.matrix(colorBy)){
			if(!is.null(colnames(colorBy))) legendTitle<- colnames(colorBy)[1]
			colorBy<-colorBy[,1]
		}

		if(is.character(colorBy)) colorBy<-as.factor(colorBy)
		d$colorBy<-colorBy;
		colorByIsFactor<-is.factor(colorBy)
		if(!is.null(emph)){
			if(!emph%in%levels(colorBy)) stop("emph value not in levels of colorBy")
			d$colorBy<-as.character(d$colorBy)
			d$colorBy[which(d$colorBy!=emph)]<-NA
			d$colorBy<-as.factor(d$colorBy)
			d$sizeMultiplier[which(d$colorBy==emph)]<-d$sizeMultiplier[which(d$colorBy==emph)]
		}
		if(na.bg){
			indexNA<-which(is.na(d$colorBy))
			indexNotNA<-which(!is.na(d$colorBy))
			if(length(indexNA)>0){
				tempd<-d
				tempd[1:length(indexNA),]<-d[indexNA,]
				tempd[(length(indexNA)+1):nrow(d),]<-d[indexNotNA,]
				d<-tempd
			}
		}
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab,color=colorBy,fill=colorBy))
		if(is.null(colorScale)){
			colorScale<-c("grey","red")
			if(colorByIsFactor) colorScale<-funAutoColorScale(nlevels(colorBy))
		}
		funColorScaleFill<-paste0("scale_fill_",ifelse(colorByIsFactor,"manual","gradientn"))
		funColorScaleColor<-paste0("scale_color_",ifelse(colorByIsFactor,"manual","gradientn"))
		paramColorScale<-list("name"=legendTitle,na.value=na.color)
		paramColorScaleType<-ifelse(colorByIsFactor,"values","colors");paramColorScale[[paramColorScaleType]]<-colorScale
		graph<-graph+do.call(funColorScaleFill,paramColorScale)
		graph<-graph+do.call(funColorScaleColor,paramColorScale)
	}
	if(plotText){
		graph<-graph+geom_text(alpha=alpha,size=d$sizeMultiplier)
	}else{
		if(is.null(colorBy)){
			graph<-graph+geom_point(stroke=pointStroke,colour = strokeColor,shape=21,alpha=alpha,fill="black",size=d$sizeMultiplier)
		}else{
			graph<-graph+geom_point(stroke=pointStroke,colour = strokeColor,shape=21,alpha=alpha,size=d$sizeMultiplier)
		}
	}
	if(plotLabelRepel){
		require(ggrepel)
		graph<-graph+geom_text_repel(color="black",size=labelSize)
	}
	graph<-graph+xlab(variable.name)
	if(!is.null(main)) graph <- graph + ggtitle(main)
	if(ellipse) graph<-graph+stat_ellipse(size=.5)
	if(plotFactorsCentroids){
		samplePerFactor<-lapply(levels(colorBy),function(x) which(colorBy==x))
		names(samplePerFactor)<-levels(colorBy)
		centroids<-data.frame(t(sapply(samplePerFactor,function(x) colMeans(d[x,c("Axis1","Axis2")]))),groupName=names(samplePerFactor))
		graph<-graph+geom_label(data = centroids,mapping = aes(x=Axis1, y=Axis2, label = groupName),inherit.aes = FALSE)
	}
	graph<-graph+theme(panel.background = element_blank(),
										 axis.text.y = element_blank(),
										 axis.ticks.y = element_blank(),
										 axis.title.y = element_blank()) + guides(size="none")
	if(returnGraph){
		return(graph)
	}else{
		print(graph)
	}
}


#' Wrapper for proj2d with highly contrasted color scale.
#'
#' @export
#'
#' @examples
#' data(iris)
#' proj2dWideColor(iris,colorBy = iris$Petal.Length)
proj2dWideColor<-function(coord, axis=1:2,colorBy=NULL, pointSize=3, plotText=FALSE,main=NULL,alpha=9/10,
													ellipse=FALSE,emph=NULL,returnGraph=FALSE,legendTitle="Values",axis.names=NULL,
													na.color="grey50",na.bg=TRUE,funAutoColorScale=ggplotColours,...){
	proj2d(coord, axis=axis,colorBy=colorBy, pointSize=pointSize, plotText=plotText,main=main,alpha=alpha,
				 ellipse=ellipse,emph=emph,colorScale=c("white","yellow","red","purple"),returnGraph=returnGraph,
				 legendTitle=legendTitle,axis.names=axis.names,na.color=na.color,na.bg=na.bg,funAutoColorScale=ggplotColours,...)
}

#' Plot a 3D projection.
#'
#' @param coord Matrix or dataframe containing at least three vectors of numeric corresponding to the x/y/z coordinates.
#' @param axis A vector of 3 integer. The column index of `coord` containing respectively the x and y coordinates.
#' @param colorBy A vector of factor or numeric same size as number of samples in PCA. Color each point on the projection by its corresponding value.
#' @param pointSize Single numeric. Size of points.
#' @param plotText Logical. Plot sample instead of point.
#' @param colorScale A vector of colors. Equal to the desired number of breaks for continuous values or number of levels for factors.
#' @param na.color Character. Color of NA values if `ColorBy` is provided.
#' @param na.bg Logical. Put the NA points behind the others in term of layers.
#' @param alpha Single numeric (from 0 to 1). Opacity of the points.
#' @param ... Other arguments passed to `plot3d`.
#'
#' @return Plot the projection in a new rgl window.
#' @export
#'
#' @examples
#' data(iris)
#'
#' proj3d(iris,pointSize = .05)
#' proj3d(iris,axis = c(1,3,4),pointSize = .05)
#'
#'
#' proj3d(iris,colorBy = iris$Species,pointSize = .05)
#'
#' proj3d(iris,colorBy = iris$Species,colorScale = c("blue","white","red"),pointSize = .05)
proj3d<-function(coord, axis=1:3, colorBy=NULL, pointSize=5, plotText=FALSE,colorScale=NULL,na.color="grey50",na.bg=TRUE,alpha=1,...){
	if(!(is.factor(colorBy)|is.numeric(colorBy)|is.null(colorBy))) stop("Error, colorBy must be numeric, factor or null.")
	if(!require("rgl")) stop("You must install rgl");
	if(!require("circlize")) stop("You must install circlize");
	if(length(axis)!=3) stop("You must give a vector of 3 integer for axis parameter")

	if(is.null(colorBy)){ colors="black"}
	else{
		if(na.bg){
			indexNA<-which(is.na(colorBy))
			indexNotNA<-which(!is.na(colorBy))
			if(length(indexNA)>0){
				tempc<-coord;tempc[1:length(indexNA),]<-coord[indexNA,];tempc[(length(indexNA)+1):nrow(coord),]<-coord[indexNotNA,];coord<-tempc;
				tempg<-colorBy;tempg[1:length(indexNA)]<-colorBy[indexNA];tempg[(length(indexNA)+1):length(colorBy)]<-colorBy[indexNotNA];colorBy<-tempg;
			}
		}
		if(is.factor(colorBy)){
			if(is.null(colorScale)) colorScale<-rainbow(nlevels(colorBy))
			hashCol<-colorScale
			names(hashCol)<-levels(colorBy)
			colors<-hashCol[colorBy]
		}else{
			if(is.null(colorScale)) colorScale<-c("grey","red")
			funCol<-colorRamp2(seq(min(colorBy),max(colorBy),length.out=length(colorScale)),colorScale)
			colors<-funCol(colorBy)
		}
	}
	colors[is.na(colors)]<-na.color
	plot3d(coord[,axis[1]],coord[,axis[2]],coord[,axis[3]],
				 xlab=paste0("Axis",axis[1]),
				 ylab=paste0("Axis",axis[2]),
				 zlab=paste0("Axis",axis[3]),
				 type="n",...)
	if(plotText){
		text3d(coord[,axis[1]],coord[,axis[2]],coord[,axis[3]],texts=rownames(coord),cex=pointSize,col=colors,alpha=alpha)
	}else{
		spheres3d(coord[,axis[1]],coord[,axis[2]],coord[,axis[3]],col=colors,radius=pointSize,alpha=alpha)
	}
	if(!is.null(colorBy)){
		if(is.factor(colorBy)) legend3d("topright", legend = names(hashCol), pch = 16, col = hashCol, cex=1, inset=c(0.02))
		if(is.numeric(colorBy)) legend3d("topright", legend = c(min(colorBy),max(colorBy)), pch = 16, col = colorScale, cex=1, inset=c(0.02))
	}
}


#' Plot a density-based distribution of a feature.
#'
#' @param x A numeric vector.
#' @param returnGraph Logical. Logical. Return the graph as a ggplot object instead of printing it.
#' @param ... Parameters passed to `qplot`
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' qplotDensity(rnorm(10000))
qplotDensity<-function(x,returnGraph=FALSE,...){
	require(ggplot2)
	if(class(x)=="data.frame" | class(x)=="list") x<-unlist(x)
	if(class(x)=="matrix") x<-as.vector(x)
	g<-qplot(x,geom="density",...)
	if(!returnGraph){
		print(g)
	}else{
		g
	}
}

#' Quick simple barplot with heights provided in ggplot format.
#'
#' @param y Vector of numeric.
#' @param returnGraph Logical. Print the ggplot object or return it.
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' qplotBarplot(1:5)
#' y<-1:5
#' names(y)<-intToUtf8(65:69,multiple = TRUE)
#' qplotBarplot(y)
qplotBarplot<-function(y,returnGraph=FALSE){
	require(ggplot2)
	if(is.null(names(y))){
		aesX = 1:length(y)
	}else{
		aesX=names(y)
	}
	g<-ggplot(data.frame(x=aesX,y=y),aes(x=x,y=y))+
		geom_bar(stat="identity")
	if(!returnGraph){
		print(g)
	}else{
		g
	}
}

#' Quick plot of values with x values corresponding to index.
#'
#' @param x Vector of numeric.
#' @param returnGraph Logical. Print the ggplot object or return it.
#' @param geom Geom to draw. Defaults to "point"
#' @param ...  Parameters passed to `qplot`
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' qplotAutoX(1:5)
qplotAutoX<-function(x,returnGraph=FALSE,geom="point",...){
	require(ggplot2)
	g<-qplot(x=1:length(x),y=x,geom=geom,...)
	if(!returnGraph){
		print(g)
	}else{
		g
	}
}

#' Scree plot of explained variance per PC from a PCA
#'
#' @param pca PCA object created by `PCA` function.
#' @param nPC Integer. Number of principal component to be plotted.
#' @param returnGraph Logical. Print the ggplot object or return it.
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' data(iris)
#' pca<-PCA(iris[,1:4],transpose = FALSE,scale = TRUE,center = TRUE)
#' barplotPercentVar(pca)
#' barplotPercentVar(pca, nPC = 2)
barplotPercentVar<-function(pca, nPC=length(pca$percentVar),returnGraph=FALSE,...){
	g<-qplotBarplot(pca$percentVar[1:nPC]*100,returnGraph=TRUE)+ylab("% variance explained")+xlab("Principal component")+
		scale_x_continuous(breaks = seq(1,nPC,by = 2))+
		theme(
			panel.background = element_rect(fill = NA,colour="black"),
			panel.grid.major = element_line(colour = "grey"),
			panel.grid.minor = element_line(colour = NA)
		)
	if(!returnGraph){
		print(g)
	}else{
		g
	}
}

#' Plot a double arrow in a grid plot.
#'
#' @param x Numeric. X coordinate of the arrow.
#' @param y Numeric. Y coordinate of the arrow.
#' @param width Numeric. X coordinate of the arrow.
#' @param height Numeric. X coordinate of the arrow.
#' @param just A string or numeric vector specifying the justification of the viewport relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification. Possible string values are: "left", "right", "centre", "center", "bottom", and "top". For numeric values, 0 means left alignment and 1 means right alignment.
#' @param gp An object of class "gpar", typically the output from a call to the function gpar. This is basically a list of graphical parameter settings
#' @param ... Other parameters passed to pushViewport.
#'
#' @return Plot in the current graphical device.
#' @export
#'
#' @examples
#' filledDoubleArrow()
#' qplot(1:5)
#' filledDoubleArrow()
filledDoubleArrow<-function(x=0,y=0,width=1,height=1,just = c("left", "bottom"),gp=gpar(col="black"),...){
	require(grid)
	pushViewport(viewport(x=x,y=y,width = width, height = height,just=just,...))
	grid.polygon(x = c(0,.2,.2,.8,.8,1,.8,.8,.2,.2,0),y=c(.5,.7,.6,.6,.7,.5,.3,.4,.4,.3,.5),gp = gp)
	popViewport()
}


#' Pot the expression of one or several genes
#'
#' @param expr Numeric 2D matrix. Each row is gene and each column a sample. Row must be named by genes.
#' @param group Factor vector, same length as number of column in expr. Experimental group attributed to each sample.
#' @param log10Plus1yScale Logical or NULL. Use a `log10(x+1)` scale. By default `FALSE` if one gene and `TRUE` if several.
#' @param violin Logical. Plot violin plots.
#' @param boxplot Logical. Plot violin plots.
#' @param dotplot Logical. Plot dot plots.
#' @param violinArgs List. Additional arguments given to `geom_violin`.
#' @param boxplotArgs List. Additional arguments given to `geom_boxplot`.
#' @param dotplotArgs List. Additional arguments given to `geom_beeswarm`.
#' @param colorScale A list of color. Must be the same length as number of levels in group.
#' @param defaultGroupName Character. Displayed in legend title.
#' @param dodge.width Numeric. Width of individual distribution element (violin/boxplot/dotplot).
#' @param returnGraph Logical. Print the ggplot object or return it.
#'
#' @return Plot in the current graphical device or a ggplot object if `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' exprMat<-rbind(
#' 	rnegbin(30,theta = 10,mu = 10),
#' 	rnegbin(30,theta = 1,mu = 10),
#' 	rnegbin(30,theta = 1,mu = 1000)
#' )
#' rownames(exprMat)<-c("gene1","gene2","gene3")
#'
#' group=c(rep("A",10),rep("B",10),rep("C",10))
#'
#' plotExpression(exprMat)
#' plotExpression(exprMat["gene1",])
#' plotExpression(exprMat["gene1",],group = group)
#' plotExpression(exprMat,group = group)
#'
#' plotExpression(exprMat,group = group,log10Plus1yScale = FALSE)
#' plotExpression(exprMat,group = group,dodge.width = .5,defaultGroupName = "Letter")
#' plotExpression(exprMat,boxplot = FALSE,violin = FALSE,dotplot = TRUE)
#' plotExpression(exprMat,group = group,colorScale = c("blue","white","red"))
#' plotExpression(exprMat,group = group,returnGraph = TRUE)+
#' 	geom_point()
#'
#' plotExpression(exprMat,group = group,violinArgs = list(scale="area"))

plotExpression<-function(expr,group=NULL,log10Plus1yScale=NULL,violin=TRUE,boxplot=TRUE,dotplot=FALSE,
												 violinArgs=list(),boxplotArgs=list(),dotplotArgs=list(),colorScale=mostDistantColor2,
												 defaultGroupName="group",dodge.width=.9,returnGraph=FALSE){
	require(ggplot2)
	require(reshape2)
	require(ggbeeswarm)
	barplotGraph=greyGraph=coloredGraph=FALSE

	if(is.vector(expr))expr<-t(as.matrix(expr))
	if(is.vector(group) | is.factor(group)){
		group = data.frame(group=group,stringsAsFactors = TRUE)
		colnames(group)<-defaultGroupName
		rownames(group)<-colnames(expr)
	}

	if(!is.matrix(expr)) expr<-as.matrix(expr)
	if(is.null(log10Plus1yScale)) log10Plus1yScale<-nrow(expr)>1 #if more than one gene log10Plus1yScale is turned on

	if(is.null(group)){
		if(nrow(expr)<2){
			if(is.null(colnames(expr))) colnames(expr)<-as.character(1:ncol(expr))
			ggData<-data.frame(expression=expr[1,],sample=factor(colnames(expr),levels = colnames(expr)[order(expr[1,],decreasing = TRUE)]))
			ggData$sample<-factor(ggData$sample)
			g<-ggplot(ggData,mapping = aes(x=sample,y=expression))+geom_bar(stat = "identity")+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3))
			barplotGraph<-TRUE
		}else{
			ggData<-reshape2::melt(expr,value.name="expression",varnames=c("gene","sample"))
			ggData$gene<-factor(ggData$gene,levels = rownames(expr))
			g<-ggplot(ggData,mapping = aes(x=gene,y=expression))+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold.italic"))
			greyGraph<-TRUE
		}
	}else{ #group on
		if(ncol(expr)!=nrow(group)) stop("expr should have the same number of sample than group")
		if(ncol(group)>1) stop("Multiple group are not allowed in the same time")
		groupName <- colnames(group)
		if(nrow(expr)==1){
			ggData<-data.frame(expression=expr[1,],group)
			g<-ggplot(ggData,mapping = aes_string(x=groupName,y="expression"))+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold"))
			greyGraph<-TRUE
		}else{
			coloredGraph=TRUE
			ggData<-reshape2::melt(data.frame(t(expr),group),value.name="expression",variable.name="gene",id.vars = groupName)
			if(violin){
				factorSampling<-table(group[,1])
				factor2drop<-names(factorSampling)[factorSampling<3]
				if(length(factor2drop)>1) warning(paste0(factor2drop,collapse = " "),
																					" were dropped (n<3 is not compatible with violin plot). You can deactivate violin layer by setting violin argument to FALSE")
				ggData<-ggData[!ggData[,groupName]%in%factor2drop,] #drop levels where n < 3

			}
			colors <- if(is.function(colorScale)) colorScale(nlevels(group[,1])) else colorScale
			g<-ggplot(ggData,mapping = aes_string(x="gene",y="expression",fill=groupName))+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold.italic"))+
				scale_fill_manual(values = colors[!levels(group[,1])%in%factor2drop])
		}
	}
	if(greyGraph){
		if(is.null(violinArgs$fill)) violinArgs$fill<-"grey50"
		if(is.null(violinArgs$scale)) violinArgs$scale<-"width"
		if(is.null(boxplotArgs$width)) boxplotArgs$width<-.2
	}
	if(coloredGraph){
		if(is.null(violinArgs$scale)) violinArgs$scale<-"width"
		if(is.null(boxplotArgs$width)) boxplotArgs$width<-.2
		if(is.null(violinArgs$position)) violinArgs$position<-position_dodge(preserve = "total",width = dodge.width)
		if(is.null(boxplotArgs$position)) boxplotArgs$position<-position_dodge(preserve = "total",width = dodge.width)
		if(is.null(dotplotArgs$dodge.width)) dotplotArgs$dodge.width<-dodge.width
	}
	if(!barplotGraph){
		g<-ggBorderedFactors(g,borderColor="black",borderSize = .5)
		if(violin) g<-g+do.call("geom_violin",violinArgs)
		if(boxplot) g<-g+do.call("geom_boxplot",boxplotArgs)
		if(dotplot)	g<-g+do.call("geom_beeswarm",dotplotArgs)
	}
	if(log10Plus1yScale){
		maxExpr<-max(expr)
		breaks<-c(0,2,round(10^(seq(1,nchar(maxExpr),0.5))))
		#breaks<-c(0,rbind(breaks/2,breaks)) #intelacing 1,10,100... and 5,50,500...
		if(maxExpr<breaks[length(breaks)-1]) breaks<-breaks[1:length(breaks)-1]

		g<-g+scale_y_continuous(trans=log10plus1(),limits = c(breaks[1],breaks[length(breaks)]),breaks = breaks,minor_breaks = NULL)
	}
	if(!returnGraph){
		print(g)
	}else{
		return(g)
	}
}


GSDA.HeatmapAnnot<-function(contributions,maxGeneContribAtOneSide=3,width=unit(3,"cm"),fontsizeFactor=400){
	AnnotationFunction(fun = function(index, k, n) {
		pushViewport(viewport(xscale = c(0,10), yscale = c(0.5, length(index) + 0.5)))
		grid.lines(.5,c(0,1))
		i<-length(index)
		for(sel in index){
			contrib<-sort(contribPerPathway[[sel]])

			left<-names(contrib)[contrib<0]
			right<-names(contrib)[contrib>0]
			left<-left[1:min(maxGeneContribAtOneSide,length(left))]
			right<-right[max(1,(length(right)-maxGeneContribAtOneSide)+1):length(right)]
			if(!is.null(left[1])){
				grid.text(paste0(left,collapse = "  "),just = "right",
									x=4.5,y=i,default.units = "native",gp=gpar(fontsize=1/length(index)*fontsizeFactor,fontface="italic"))
			}
			if(!is.null(right[1])){
				grid.text(paste0(right,collapse = "  "),just = "left",
									x=5.5,y=i,default.units = "native",gp=gpar(fontsize=1/length(index)*fontsizeFactor,fontface="italic"))
			}
			i<-i-1
		}
		popViewport()
	},
	var_import = list(contribPerPathway = contributions, maxGeneContribAtOneSide=maxGeneContribAtOneSide,width=width,
										fontsizeFactor=fontsizeFactor),
	subsetable = FALSE,
	width = width,which = "row"
	)
}

viewKEGG<-function(x,pathway,corrIdGenes=NULL,species="Human",speciesData=NULL,directory=getwd(),...){
	require(pathview)

	blacklist<-c("hsa04215 Apoptosis - multiple species")
	if(pathway%in%blacklist){
		warning(pathway," is blacklisted as it contains issues in vizualisation, it will not be rendered.")
		return(NULL)
	}

	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}
	if(is.null(speciesData)) speciesData<-getSpeciesData2(species)
	if(is.null(corrIdGenes)) corrIdGenes<-speciesData$GeneIdTable
	entrezId<-ConvertKey(keyList = names(x),tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID" );
	notNA<-which(!is.na(entrezId))
	if(length(notNA)>0){
		dat<-x[notNA];
		names(dat)<-entrezId[notNA]
		dat<-dat[takefirst(names(dat),returnIndex = T)]


		pathview(gene.data = dat, pathway.id = pathway, species = speciesData$kegg,kegg.native=TRUE,
						 low="#4B9AD5",mid="white",high="#FAB517",na.col="grey75",kegg.dir=directory,...)
	}else{
		warning("no entrez id were found")
	}
}



genTopAnnot<-function(annot,colorScales=NULL,border = TRUE,...){
	require(ComplexHeatmap)
	require(circlize)

	if(is.factor(annot)|is.data.frame(annot)) annot<-droplevels(annot)
	if(is.vector(annot) | is.factor(annot)){
		annot=data.frame(Annotation=annot)
		if(is.list(colorScales)) colnames(annot)<-names(colorScales)[1]
		if( (!is.null(colorScales)) & !is.list(colorScales)) colorScales<-list("Annotation"=colorScales)
	}
	colorScales<-genColorsForAnnots(annots = annot,colorScales = colorScales,returnContinuousFun = TRUE,...)
	HeatmapAnnotation(df=annot,col = colorScales,border = border)
}


heatmap.DM3<-function(matrix,preSet="expression",clustering_distance_rows=NULL,clustering_distance_columns="euclidean",clustering_method_columns="ward.D2",
											clustering_method_rows="ward.D2",autoFontSizeRow=TRUE,autoFontSizeColumn=TRUE,scale=FALSE,center=TRUE,returnHeatmap=FALSE,name=NULL,
											additionnalRowNamesGpar=NULL,additionnalColNamesGpar=list(),border=TRUE,
											colorScale=NULL,colorScaleFun=NULL,midColorIs0=NULL,probs=NULL,useProb=TRUE,minProb=0.05,maxProb=0.95,
											cluster_rows=NULL,cluster_columns=NULL,sampleAnnot=NULL,colorAnnot=NULL,showGrid=NULL,gparGrid=gpar(col="black"),showValues=FALSE,Nsignif=3,
											column_dend_reorder = FALSE, row_dend_reorder=FALSE, ...){
	require("ComplexHeatmap")
	require("circlize")

	args<-list()
	allowedPreSet<-c("default","expression","correlation","distance")
	if(!preSet%in%allowedPreSet){
		stop(paste0("preSet must equal to one of this value: ",paste0(allowedPreSet,collapse = ", ")))
	}

	if(preSet=="expression"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows = corrDistBicor
		if(is.null(name)) name="centered log expression"
		if(is.null(midColorIs0)) midColorIs0=TRUE
		if(is.null(colorScale)) colorScale=	c("darkblue","white","red2")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list(fontface="italic")
	}else if(preSet=="correlation"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows ="euclidean"
		if(is.null(clustering_distance_columns)) clustering_distance_columns ="euclidean"
		if(is.null(name)) name="Pearson\ncorrelation"
		if(is.null(midColorIs0)) midColorIs0=TRUE
		if(is.null(colorScale)) colorScale=c("darkblue","white","#FFAA00")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list()
	}else if(preSet=="distance"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows ="euclidean"
		if(is.null(name)) name="Euclidean\ndistance"
		if(is.null(midColorIs0)) midColorIs0=FALSE
		if(is.null(colorScale)) colorScale=c("white","yellow","red","purple")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list()
	}else	if(preSet=="default"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows = corrDistBicor
		if(is.null(name)) name="matrix"
		if(is.null(midColorIs0)) midColorIs0=TRUE
		if(is.null(colorScale)) colorScale=c("#2E3672","#4B9AD5","white","#FAB517","#E5261D")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list()
	}

	if(is.null(cluster_rows)){
		args$clustering_method_rows<-clustering_method_rows
		args$clustering_distance_rows<-clustering_distance_rows
	}else{
		args$cluster_rows<-cluster_rows
	}
	if(is.null(cluster_columns)){
		args$clustering_method_columns<-clustering_method_columns
		args$clustering_distance_columns<-clustering_distance_columns
	}else{
		args$cluster_columns<-cluster_columns
	}
	matrix<-as.matrix(matrix)
	if(min(apply(matrix,1,sd))==0 & (scale | identical(corrDist,clustering_distance_rows)) ){
		warning("some row have a 0 sd. sd-based method (correlation distance, scaling) will be desactivated or switched.")
		scale=FALSE
		if(identical(corrDist,clustering_distance_rows)){
			args$clustering_distance_rows<-"euclidean"
		}
	}

	if(scale | center) matrix<-rowScale(matrix,scaled=scale,center=center)
	if(is.null(colorScaleFun)){
		colorScaleFun<-computeColorScaleFun(colors = colorScale,values = unlist(matrix),useProb = useProb,probs = probs,minProb = minProb,
																				maxProb = maxProb, midColorIs0 = midColorIs0,returnColorFun = TRUE)
	}
	args$col<-colorScaleFun
	if(is.null(showGrid)){
		if(nrow(matrix)*ncol(matrix)<500){
			showGrid = TRUE
		}else{
			showGrid = FALSE
		}
	}
	if(showGrid){
		args$rect_gp = gparGrid
	}
	if(showValues){
		args$cell_fun = function(j, i, x, y, w, h, col) {
			#dark or light background .
			if(colSums(col2rgb(col))<382.5) col="white" else col="black"
			grid.text(as.character(signif(matrix[i,j],Nsignif)),x,y,gp=gpar(col=col))
		}
	}
	if(autoFontSizeRow) args$row_names_gp=do.call("autoGparFontSizeMatrix",c(list(nrow(matrix)),additionnalRowNamesGpar))
	if(autoFontSizeColumn) args$column_names_gp=do.call("autoGparFontSizeMatrix",c(list(ncol(matrix)),additionnalColNamesGpar))

	if(!is.null(sampleAnnot)){
		args$top_annotation<-genTopAnnot(sampleAnnot,colorAnnot)
	}

	args$column_dend_reorder<-column_dend_reorder;args$row_dend_reorder<-row_dend_reorder
	args$matrix<-matrix
	args$name=name
	args$border<-border
	args<-c(args,list(...))

	ht<-do.call("Heatmap",args)
	if(returnHeatmap){
		return(ht)
	}else{
		print(ht)
	}
}

heatmapSquareMatrix<-function(matrix,clustering_distance="euclidean",clustering_method="ward.D2",
															autoFontSizeRow=TRUE,autoFontSizeColumn=TRUE,scale=FALSE,center=FALSE,returnHeatmap=FALSE,name="Pearson\ncorrelation",
															additionnalRowNamesGpar=list(),additionnalColNamesGpar=list(),border=TRUE,
															colorScale=c("darkblue","white","#FFAA00"),colorScaleFun=NULL,midColorIs0=TRUE,probs=NULL,useProb=TRUE,minProb=0.05,maxProb=0.95,
															cluster=NULL,cluster_columns=NULL,sampleAnnot=NULL,colorAnnot=NULL,showGrid=NULL,gparGrid=gpar(col="black"),showValues=FALSE,Nsignif=3, ...){

	require("ComplexHeatmap")
	require("circlize")

	matrix<-as.matrix(matrix)
	args<-list()
	if(is.null(cluster)){
		cluster<-unsupervisedClustering(matrix,transpose = FALSE,method.dist = clustering_distance,method.hclust = clustering_method)
	}
	args$cluster_rows<-cluster
	args$cluster_columns<-cluster

	if(min(apply(matrix,1,sd))==0 & scale ){
		warning("some row have a 0 sd. sd-based method (correlation distance, scaling) will be desactivated or switched.")
		scale=FALSE
	}

	if(scale | center) matrix<-rowScale(matrix,scaled=scale,center=center)
	if(is.null(colorScaleFun)){
		colorScaleFun<-computeColorScaleFun(colors = colorScale,values = unlist(matrix),useProb = useProb,probs = probs,minProb = minProb,
																				maxProb = maxProb, midColorIs0 = midColorIs0,returnColorFun = TRUE)
	}
	args$col<-colorScaleFun
	if(is.null(showGrid)){
		if(nrow(matrix)*ncol(matrix)<500){
			showGrid = TRUE
		}else{
			showGrid = FALSE
		}
	}
	if(showGrid){
		args$rect_gp = gparGrid
	}
	if(showValues){
		args$cell_fun = function(j, i, x, y, w, h, col) {
			#dark or light background .
			if(colSums(col2rgb(col))<382.5) col="white" else col="black"
			grid.text(as.character(signif(matrix[i,j],Nsignif)),x,y,gp=gpar(col=col))
		}
	}
	if(autoFontSizeRow) args$row_names_gp=do.call("autoGparFontSizeMatrix",c(list(nrow(matrix)),additionnalRowNamesGpar))
	if(autoFontSizeColumn) args$column_names_gp=do.call("autoGparFontSizeMatrix",c(list(ncol(matrix)),additionnalColNamesGpar))

	if(!is.null(sampleAnnot)){
		args$top_annotation<-genTopAnnot(sampleAnnot,colorAnnot)
	}

	args$matrix<-matrix
	args$name=name
	args$border<-border
	args<-c(args,list(...))

	ht<-do.call("Heatmap",args)
	if(returnHeatmap){
		return(ht)
	}else{
		print(ht)
	}
}


volcanoPlot.DESeq2<-function(DEresult,formula,downLevel,upLevel,condColumn,padjThreshold,LFCthreshold,topGene=30){
	require(ggrepel)
	require(grid)
	DEresult<-DEresult[!is.na(DEresult$padj),]
	gene2Plot<-order(DEresult$padj)
	gene2Plot<-gene2Plot[DEresult[gene2Plot,"isDE"]!="NONE"]
	gene2Plot<-gene2Plot[1:min(topGene,length(gene2Plot))]
	g<-ggplot(DEresult,aes(x=log2FoldChange,y=-log10(padj)+0.01,color=isDE))+
		geom_point(size=1)+theme_bw()+scale_color_manual(values=c("#3AAA35","grey75","#E40429"))+
		geom_text_repel(data = DEresult[gene2Plot,],aes(x=log2FoldChange,y=-log10(padj),label=gene),
										inherit.aes = FALSE,color="black",fontface = "bold.italic",size=3)+
		ylab("-log10(adjusted pvalue)")+xlab(NULL)+
		geom_vline(xintercept = c(-LFCthreshold,LFCthreshold))+
		geom_hline(yintercept = -log10(padjThreshold)) + guides(color = FALSE)+
		ggtitle("Volcano plot")

	grid.newpage();
	pushViewport(viewport(x = 0, y = 0,
												width = .8, height = .1,
												just = c("left", "bottom")))
	filledDoubleArrow(x=.3,y=1,width = .3,just=c("left","center"),gp = gpar(fill="black"))
	grid.text(label = downLevel,x = .28,y = 1.03,just=c("right","center"),gp = gpar(fontface="bold"))
	grid.text(label = upLevel,x = .62,y = 1.03,just=c("left","center"),gp = gpar(fontface="bold"))
	grid.text(label = "log2(Fold-Change)",x = .45,y = 1.5,just = c("center","center"))
	popViewport()
	pushViewport(viewport(x = .65, y = 1,
												width = .35, height = 1,
												just = c("left", "top"),default.units = "npc"))
	grid.text(label = paste0("Experimental design:\n",formula),x = 0,y = .95,just=c("left","center"))
	grid.text(label = paste0("Results for\n",condColumn,":\n",downLevel," vs ",upLevel),x = .0,y = .8,just=c("left","center"))
	grid.text(label = paste0(sum(DEresult$isDE=="DOWNREG")," downreg. genes"),x = 0,y = .65,just=c("left","center"),
						gp = gpar(col="#3AAA35",fontface="bold"))

	grid.text(label = paste0(sum(DEresult$isDE=="UPREG")," upreg. genes"),x =.0,y = .55,just=c("left","center"),
						gp = gpar(col="#E40429",fontface="bold"))
	grid.text(label = paste0("From ",nrow(DEresult),"\ntested genes"),x = 0,y = .45,just=c("left","center"))
	popViewport()
	main_vp <- viewport(x = 0, y = 1,
											width = .8, height = .9,
											just = c("left", "top"))
	pushViewport(main_vp);print(g,vp=main_vp);popViewport()
}



customUpsetPlot<-function(genePerGroupList,universe=NULL){
	require(ComplexHeatmap)
	if(is.null(universe)) universe<- unique(unlist(genePerGroupList))
	isInGroupMatrix<-list_to_matrix(genePerGroupList,universal_set = universe)
	upsetMatrix<-make_comb_mat(isInGroupMatrix,mode = "intersect")
	upsetMatrix<-upsetMatrix[comb_degree(upsetMatrix) > 1] # retain only intersections of sets

	combsize = comb_size(upsetMatrix)
	setsize = set_size(upsetMatrix)

	#Are the intersections sets (or venn diagramm region) enriched or not ?
	regionEnrich<-sapply(comb_name(upsetMatrix),function(region){
		colOfcomp=which(strsplit(region,split = "")[[1]]=="1")
		intersectionEnrichment(isInGroupMatrix[,colOfcomp])
	})

	enrich_ha = HeatmapAnnotation(
		"enrichment" = anno_barplot(
			regionEnrich, gp = gpar(fill = "black"), height = unit(3, "cm"),axis_param = list(side = "left"),
			ylim = c(0, max(regionEnrich)*1.1)
		),
		annotation_name_side = "left", annotation_name_rot = 0,annotation_name_gp = gpar(fontface="bold"),
		annotation_label = "Enrichment\n(real/expected size)"
	)
	intersect_ha = HeatmapAnnotation(
		"intersection_size" = anno_barplot(
			combsize, gp = gpar(fill = "black"), height = unit(3, "cm"),axis_param = list(side = "left"),
			ylim = c(0, max(combsize)*1.1)
		),
		annotation_name_side = "left", annotation_name_rot = 0,annotation_name_gp = gpar(fontface="bold"),
		annotation_label = "Intersection\nsize"
	)
	set_size_ha = rowAnnotation(
		"set_size" = anno_barplot(
			setsize,gp = gpar(fill = "black"),width = unit(2, "cm"),
			ylim = c(0, max(setsize)*1.3)
		),
		annotation_name_side = "bottom", annotation_name_rot = 0,annotation_name_gp = gpar(fontface="bold"),
		annotation_label = "Set\nsize"
	)


	ht = draw(UpSet(upsetMatrix,top_annotation = intersect_ha,bottom_annotation = enrich_ha,right_annotation = set_size_ha,
									border=TRUE,column_split=comb_degree(upsetMatrix),
									row_names_gp=gpar(fontsize=min(1/max(nchar(rownames(upsetMatrix)))*260,20))#automatic fontsize to avoid out of bound text
	))


	#Offset to counterbalance column split space
	colPerSplit=sapply(column_order(ht),length)
	offsetPerSplit=seq(0,length(colPerSplit)-1)
	offsets<-unlist(lapply(seq_along(colPerSplit),function(i) rep(offsetPerSplit[i],colPerSplit[i])),use.names = FALSE)

	rowOrder = rev(row_order(ht))
	columnOrder = unlist(column_order(ht))

	decorate_annotation("intersection_size", {
		grid.text(combsize[columnOrder], x = unit(seq_along(combsize),"native")+unit(offsets,"mm"),
							y = unit(combsize[columnOrder], "native") + unit(6, "pt"),
							default.units = "native", just = "center", gp = gpar(fontsize = 8))
	})
	decorate_annotation("enrichment", {
		grid.text(round(regionEnrich[columnOrder],2), x = unit(seq_along(regionEnrich),"native")+unit(offsets,"mm"),
							y = unit(regionEnrich[columnOrder], "native") + unit(6, "pt"),
							default.units = "native", just = "center", gp = gpar(fontsize = 8))
	})
	decorate_annotation("set_size", {
		grid.text(round(setsize[rowOrder],2), y = seq_along(setsize), x = unit(setsize[rowOrder], "native") + unit(7, "pt"),
							default.units = "native", just = "center", gp = gpar(fontsize = 10),rot=-90)
	})
}


