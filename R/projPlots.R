
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
#' pca <-
#'     PCA(iris[, seq_len(4)],
#'         transpose = FALSE,
#'         scale = TRUE,
#'         center = TRUE)
#' pca2d(pca)
#' pca2d(pca, comp = c(1, 3))
#'
#' pca2d(pca,
#'       colorBy = iris$Species,
#'       ratioPerPercentVar = TRUE)
#' pca2d(pca,
#'       colorBy = iris$Species,
#'       colorScale = c("blue", "white", "red"))
#'
#' pca2d(pca,
#'       plotVars = TRUE,
#'       plotText = TRUE,
#'       pointSize = 5)
#'
#' pca2d(pca, outlierLabel = TRUE, colorBy = iris$Species)
#' pca2d(
#'     pca,
#'     outlierLabel = TRUE,
#'     colorBy = iris$Species,
#'     ellipse = TRUE,
#'     outlierLabelThres = .2
#' )
#'
#' pca2d(
#'     pca,
#'     outlierLabel = TRUE,
#'     colorBy = iris$Species,
#'     ellipse = TRUE,
#'     outlierLabelThres = .2,
#'     returnGraph = TRUE
#' ) + theme_dark()
pca2d <-
	function(pca,
					 comp = c(1, 2),
					 colorBy = NULL,
					 plotVars = FALSE,
					 pointSize = 2,
					 plotText = FALSE,
					 ratioPerPercentVar = FALSE,
					 main = NULL,
					 ellipse = FALSE,
					 colorScale = NULL,
					 returnGraph = FALSE,
					 outlierLabel = FALSE,
					 outlierLabelThres = NULL,
					 outlierLabelSize = 3,
					 customRatio = NULL,
					 ...) {
		if (length(comp) != 2)
			stop("You must give a vector of 2 integer for comp parameter")

		propExplVar <- pca$propExplVar
		if (ratioPerPercentVar) {
			customRatio <- propExplVar[comp[2]] / propExplVar[comp[1]] * 2
		}
		if (pca$isFastPCA) {
			varianceMesage <- "% variance of computed PC"
		} else{
			varianceMesage <- "% total variance"
		}
		dt <- pca[[ifelse(plotVars, "rotation", "x")]]
		graph <-
			proj2d(
				coord = dt,
				axis = comp,
				colorBy = colorBy,
				returnGraph = TRUE,
				colorScale = colorScale,
				main = main,
				plotText = plotText,
				ellipse = ellipse,
				pointSize = pointSize,
				customRatio = customRatio,
				fixedCoord = FALSE,
				...
			) +
			xlab(paste0("PC", comp[1], ": ", round(propExplVar[comp[1]] * 100), varianceMesage)) +
			ylab(paste0("PC", comp[2], ": ", round(propExplVar[comp[2]] * 100), varianceMesage))
		if (outlierLabel) {
			if (is.null(rownames(dt)))
				rownames(dt) <- as.character(seq_len(nrow(dt)))
			dtOutlier <-
				data.frame(x = dt[, comp[1]],
									 y = dt[, comp[2]],
									 rname = rownames(dt))
			dens <- pointdensity.nrd(dtOutlier[, c(1, 2)])
			if (is.null(outlierLabelThres))
				outlierLabelThres = Mode(dens) / 3
			dtOutlier <-
				dtOutlier[dens / max(dens) < outlierLabelThres,]
			graph <-
				graph + ggrepel::geom_text_repel(
					data = dtOutlier,
					mapping = aes(
						x = .data$x,
						y = .data$y,
						label = .data$rname
					),
					size = outlierLabelSize,
					fontface = "italic",
					inherit.aes = FALSE
				)
		}

		if (returnGraph) {
			return(graph)
		} else{
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
#' proj2d(iris, axis = c(1, 3))
#'
#' proj2d(iris, pointSize = .1)
#' proj2d(iris,
#'        axis.names = c("Width of sepal", "Length of sepal") ,
#'        main = "Two features from iris")
#'
#'
#' proj2d(
#'     iris,
#'     colorBy = iris$Species,
#'     ellipse = TRUE,
#'     plotFactorsCentroids = TRUE,
#'     legendTitle = "Species"
#' )
#' proj2d(
#'     iris,
#'     colorBy = iris$Species,
#'     pointStroke = 3,
#'     strokeColor = "grey"
#' )
#'
#' proj2d(
#'     iris,
#'     colorBy = iris$Species,
#'     emph = "setosa",
#'     na.color = "white"
#' )
#' proj2d(iris, colorBy = iris$Species, returnGraph = TRUE) +
#'     geom_vline(xintercept = 6)
#'
#' proj2d(iris,
#'        colorBy = iris$Species,
#'        colorScale = c("blue", "white", "red"))
#' proj2d(iris,
#'        colorBy = iris$Species,
#'        funAutoColorScale = mostDistantColor)
#' proj2d(
#'     iris,
#'     colorBy = iris$Petal.Length,
#'     legendTitle = "Petal.Length",
#'     fixedCoord = FALSE
#' )
#' proj2d(
#'     iris,
#'     colorBy = iris$Petal.Length,
#'     legendTitle = "Petal.Length",
#'     colorScale = c("blue", "white", "red")
#' )
#'
#' proj2d(iris, plotLabelRepel = TRUE)
#' proj2d(iris, sizeProp2Dens = TRUE)
#' proj2d(iris, sizeProp2Dens = TRUE, densEps = 2)
#'
#' umap <-
#'     UMAP(
#'         scale(iris[, c(seq_len(4))]),
#'         ret_nn = TRUE,
#'         transpose = FALSE,
#'         n_neighbors = nrow(iris)
#'     )
#'
#' proj2d(umap)
#' proj2d(umap,
#'        colorBy = iris$Species ,
#'        nnMatrix = umap$nn$euclidean$idx[, seq_len(3)])
#' proj2d(
#'     umap,
#'     colorBy = iris$Species,
#'     useScatterMore = TRUE,
#'     pointSize = 5
#' )
proj2d <-
	function(coord,
					 colorBy = NULL,
					 axis = c(1, 2),
					 pointSize = 3,
					 plotText = FALSE,
					 main = NULL,
					 alpha = 9 / 10,
					 ellipse = FALSE,
					 emph = NULL,
					 colorScale = NULL,
					 returnGraph = FALSE,
					 legendTitle = "Values",
					 axis.names = NULL,
					 na.color = "grey50",
					 na.bg = TRUE,
					 plotFactorsCentroids = FALSE,
					 pointStroke = 1 / 8,
					 strokeColor = "black",
					 funAutoColorScale = oobColors,
					 fixedCoord = TRUE,
					 plotLabelRepel = FALSE,
					 labelSize = 3,
					 sizeProp2Dens = FALSE,
					 densEps = 1,
					 nnMatrix = NULL,
					 nnSegmentParam = list(alpha = .75, size = .1),
					 useScatterMore = FALSE,
					 customRatio = NULL) {
		coord <- data.frame(coord)
		if (length(axis) != 2)
			stop("You must give a vector of 2 integer for axis parameter")

		if ((!is.null(axis.names)) &
				length(axis.names) != 2)
			stop("Error, if given axis.names parameter must contain 2 values.")

		if (!is.null(customRatio))
			fixedCoord <- TRUE

		if (is.null(rownames(coord)))
			rownames(coord) <- as.character(seq_len(nrow(coord)))
		d <-
			data.frame(
				lab = rownames(coord),
				Axis1 = coord[, axis[1]],
				Axis2 = coord[, axis[2]],
				sizeMultiplier = rep(pointSize, nrow(coord))
			)

		if (sizeProp2Dens) {
			dens <- pointdensity.nrd(d[, c("Axis1", "Axis2")], eps = densEps)
			d$sizeMultiplier <-
				d$sizeMultiplier * (1 - dens / max(dens)) * 2
		}
		if (is.null(colorBy)) {
			graph <- ggplot(data = d,
											mapping = aes(
												x = .data$Axis1,
												y = .data$Axis2,
												label = .data$lab
											))
		} else{
			if (is.data.frame(colorBy) | is.matrix(colorBy)) {
				if (!is.null(colnames(colorBy)))
					legendTitle <- colnames(colorBy)[1]
				colorBy <- colorBy[, 1]
			}
			if(is.logical(colorBy)) {
				colorBy <- as.factor(colorBy)
			}
			if (is.character(colorBy))
				colorBy <- as.factor(colorBy)
			d$colorBy <- colorBy

			colorByIsFactor <- is.factor(colorBy)
			if (!is.null(emph)) {
				if (!emph %in% levels(colorBy))
					stop("emph value not in levels of colorBy")
				d$colorBy <- as.character(d$colorBy)
				d$colorBy[which(d$colorBy != emph)] <- NA
				d$colorBy <- as.factor(d$colorBy)
				d$sizeMultiplier[which(d$colorBy == emph)] <-
					d$sizeMultiplier[which(d$colorBy == emph)]
			}
			if (na.bg) {
				indexNA <- which(is.na(d$colorBy))
				indexNotNA <- which(!is.na(d$colorBy))
				if (length(indexNA) > 0) {
					tempd <- d
					tempd[seq_along(indexNA),] <- d[indexNA,]
					tempd[seq(length(indexNA) + 1, nrow(d), 1),] <-
						d[indexNotNA,]
					d <- tempd
				}
			}
			graph <-
				ggplot(
					data = d,
					mapping = aes(
						x = .data$Axis1,
						y = .data$Axis2,
						label = .data$lab,
						color = .data$colorBy,
						fill = .data$colorBy
					)
				)
			if (is.null(colorScale)) {
				colorScale <- c("grey", "red")
				if (colorByIsFactor)
					colorScale <-
						funAutoColorScale(nlevels(colorBy))
			}
			funColorScaleFill <-
				paste0("scale_fill_",
							 ifelse(colorByIsFactor, "manual", "gradientn"))
			funColorScaleColor <-
				paste0("scale_color_",
							 ifelse(colorByIsFactor, "manual", "gradientn"))
			paramColorScale <-
				list("name" = legendTitle, na.value = na.color)
			paramColorScaleType <-
				ifelse(colorByIsFactor, "values", "colors")
			paramColorScale[[paramColorScaleType]] <- colorScale
			graph <-
				graph + do.call(funColorScaleFill, paramColorScale)
			graph <-
				graph + do.call(funColorScaleColor, paramColorScale)
		}
		if (!is.null(nnMatrix)) {
			if (!is.matrix(nnMatrix))
				stop("If nnMatrix is not null, it should be a matrix !")
			retainCoord <- as.matrix(coord[, axis])
			segmentsCoord <-
				data.frame(do.call("rbind", lapply(seq_len(
					nrow(retainCoord)
				), function(i) {
					t(vapply(nnMatrix[i,], FUN.VALUE = numeric(4), function(x) {
						c(retainCoord[i,], retainCoord[x,])
					}))
				})))
			colnames(segmentsCoord) <- c("x", "y", "xend", "yend")
			graph <- graph + do.call("geom_segment",
															 c(
															 	list(
															 		data = segmentsCoord,
															 		mapping = aes(
															 			x = .data$x,
															 			y = .data$y,
															 			xend = .data$xend,
															 			yend = .data$yend
															 		),
															 		inherit.aes = FALSE
															 	),
															 	nnSegmentParam
															 ))
		}
		if (plotText) {
			graph <- graph + geom_text(alpha = alpha, size = d$sizeMultiplier)
		} else{
			if (useScatterMore) {
				ratioX <- ratioY <- 1
				if (fixedCoord) {
					if (is.null(customRatio)) {
						rangeA1 <- range(d$Axis1)
						rangeA2 <- range(d$Axis2)
						ratio <-
							(rangeA1[2] - rangeA1[1]) / (rangeA2[2] - rangeA2[1])
					} else{
						ratio <- 1 / (customRatio)
					}
					if (ratio > 1) {
						ratioX <- ratio
					} else{
						ratioY <- 1 / ratio
					}
				}
				if (is.null(colorBy)) {
					graph <-
						graph + geom_scattermore(
							pointsize = pointSize * 2,
							pixels = c(1000 * ratioX, 1000 * ratioY)
						)
				} else{
					graph <-
						graph + geom_scattermore(
							pointsize = pointSize * 2,
							mapping = aes(color = .data$colorBy),
							pixels = c(1000 * ratioX, 1000 * ratioY)
						)
				}
			} else{
				if (is.null(colorBy)) {
					graph <-
						graph + geom_point(
							stroke = pointStroke,
							colour = strokeColor,
							shape = 21,
							alpha = alpha,
							fill = "black",
							size = d$sizeMultiplier
						)
				} else{
					graph <-
						graph + geom_point(
							stroke = pointStroke,
							colour = strokeColor,
							shape = 21,
							alpha = alpha,
							size = d$sizeMultiplier
						)
				}
			}
		}
		if (plotLabelRepel) {
			graph <- graph + geom_text_repel(color = "black", size = labelSize)
		}
		if (is.null(axis.names)) {
			if (is.null(colnames(coord))) {
				axis.names <- paste0("Axis", axis)
			} else{
				axis.names <- colnames(coord[, c(axis[1], axis[2])])
			}
		}
		graph <- graph + xlab(axis.names[1]) + ylab(axis.names[2])
		if (!is.null(main))
			graph <- graph + ggtitle(main)
		if (ellipse)
			graph <- graph + stat_ellipse(size = .5)
		if (fixedCoord)
			graph <-
			graph + coord_fixed(ratio = ifelse(is.null(customRatio), 1, customRatio))
		if (plotFactorsCentroids) {
			samplePerFactor <-
				lapply(levels(colorBy), function(x)
					which(colorBy == x))
			names(samplePerFactor) <- levels(colorBy)
			centroids <-
				data.frame(t(sapply(samplePerFactor, function(x)
					colMeans(d[x, c("Axis1", "Axis2")]))), groupName = names(samplePerFactor))
			graph <-
				graph + geom_label(
					data = centroids,
					mapping = aes(
						x = .data$Axis1,
						y = .data$Axis2,
						label = .data$groupName
					),
					inherit.aes = FALSE
				)
		}
		graph <- graph + theme(
			panel.background = element_rect(fill = NA, colour = "black"),
			panel.grid.major = element_line(colour = NA),
			panel.grid.minor = element_line(colour = NA)
		) + guides(size = "none")
		if (returnGraph) {
			return(graph)
		} else{
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
#' proj1d(iris$Sepal.Length, colorBy = iris$Species)
#'
#' proj1d(iris$Sepal.Length, pointSize = .1)
#' proj1d(iris$Sepal.Length, variable.name = "Length of sepal" , main = "Two features from iris")
#'
#'
#' proj1d(
#'     iris$Sepal.Length,
#'     colorBy = iris$Species,
#'     plotFactorsCentroids = TRUE,
#'     legendTitle = "Species"
#' )
#' proj1d(
#'     iris$Sepal.Length,
#'     colorBy = iris$Species,
#'     pointStroke = 3,
#'     strokeColor = "grey"
#' )
#'
#' proj1d(
#'     iris$Sepal.Length,
#'     colorBy = iris$Species,
#'     emph = "setosa",
#'     na.color = "white"
#' )
#' proj1d(iris$Sepal.Length,
#'        colorBy = iris$Species,
#'        returnGraph = TRUE) +
#'     geom_vline(xintercept = 6)
#'
#' proj1d(
#'     iris$Sepal.Length,
#'     colorBy = iris$Species,
#'     colorScale = c("blue", "white", "red")
#' )
#' proj1d(iris$Sepal.Length,
#'        colorBy = iris$Species,
#'        funAutoColorScale = mostDistantColor)
#' proj1d(
#'     iris$Sepal.Length,
#'     colorBy = iris$Petal.Length,
#'     legendTitle = "Petal.Length",
#'     colorScale = c("blue", "white", "red")
#' )
#'
#' proj1d(iris$Sepal.Length, plotLabelRepel = TRUE)
#' proj1d(iris$Sepal.Length, sizeProp2Dens = TRUE)
#' proj1d(iris$Sepal.Length,
#'        sizeProp2Dens = TRUE,
#'        densEps = 2)
proj1d <-
	function(variable,
					 colorBy = NULL,
					 pointSize = 3,
					 plotText = FALSE,
					 main = NULL,
					 alpha = 9 / 10,
					 ellipse = FALSE,
					 emph = NULL,
					 colorScale = NULL,
					 returnGraph = FALSE,
					 legendTitle = "Values",
					 variable.name = "x",
					 na.color = "grey50",
					 na.bg = TRUE,
					 plotFactorsCentroids = FALSE,
					 pointStroke = 1 / 8,
					 strokeColor = "black",
					 funAutoColorScale = oobColors,
					 plotLabelRepel = FALSE,
					 labelSize = 3,
					 sizeProp2Dens = FALSE,
					 densEps = 1) {
		if (!is.numeric(variable))
			stop("'variable' must be a numeric vector")

		if (is.null(names(variable)))
			names(variable) <- as.character(seq_along(variable))
		d <-
			data.frame(
				lab = names(variable),
				Axis1 = variable,
				Axis2 = rep(0, length(variable)),
				sizeMultiplier = rep(pointSize, length(variable))
			)

		if (sizeProp2Dens) {
			dens <- pointdensity.nrd(d[, c("Axis1", "Axis2")], eps = densEps)
			d$sizeMultiplier <-
				d$sizeMultiplier * (1 - dens / max(dens)) * 2
		}
		if (is.null(colorBy)) {
			graph <- ggplot(data = d,
											mapping = aes(
												x = .data$Axis1,
												y = .data$Axis2,
												label = .data$lab
											))
		} else{
			if (is.data.frame(colorBy) | is.matrix(colorBy)) {
				if (!is.null(colnames(colorBy)))
					legendTitle <- colnames(colorBy)[1]
				colorBy <- colorBy[, 1]
			}

			if(is.logical(colorBy)) {
				colorBy <- as.factor(colorBy)
			}
			if (is.character(colorBy))
				colorBy <- as.factor(colorBy)
			d$colorBy <- colorBy

			colorByIsFactor <- is.factor(colorBy)
			if (!is.null(emph)) {
				if (!emph %in% levels(colorBy))
					stop("emph value not in levels of colorBy")
				d$colorBy <- as.character(d$colorBy)
				d$colorBy[which(d$colorBy != emph)] <- NA
				d$colorBy <- as.factor(d$colorBy)
				d$sizeMultiplier[which(d$colorBy == emph)] <-
					d$sizeMultiplier[which(d$colorBy == emph)]
			}
			if (na.bg) {
				indexNA <- which(is.na(d$colorBy))
				indexNotNA <- which(!is.na(d$colorBy))
				if (length(indexNA) > 0) {
					tempd <- d
					tempd[seq_along(indexNA),] <- d[indexNA,]
					tempd[(length(indexNA) + 1):nrow(d),] <-
						d[indexNotNA,]
					d <- tempd
				}
			}
			graph <-
				ggplot(
					data = d,
					mapping = aes(
						x = .data$Axis1,
						y = .data$Axis2,
						label = .data$lab,
						color = .data$colorBy,
						fill = .data$colorBy
					)
				)
			if (is.null(colorScale)) {
				colorScale <- c("grey", "red")
				if (colorByIsFactor)
					colorScale <-
						funAutoColorScale(nlevels(colorBy))
			}
			funColorScaleFill <-
				paste0("scale_fill_",
							 ifelse(colorByIsFactor, "manual", "gradientn"))
			funColorScaleColor <-
				paste0("scale_color_",
							 ifelse(colorByIsFactor, "manual", "gradientn"))
			paramColorScale <-
				list("name" = legendTitle, na.value = na.color)
			paramColorScaleType <-
				ifelse(colorByIsFactor, "values", "colors")
			paramColorScale[[paramColorScaleType]] <- colorScale
			graph <-
				graph + do.call(funColorScaleFill, paramColorScale)
			graph <-
				graph + do.call(funColorScaleColor, paramColorScale)
		}
		if (plotText) {
			graph <- graph + geom_text(alpha = alpha, size = d$sizeMultiplier)
		} else{
			if (is.null(colorBy)) {
				graph <-
					graph + geom_point(
						stroke = pointStroke,
						colour = strokeColor,
						shape = 21,
						alpha = alpha,
						fill = "black",
						size = d$sizeMultiplier
					)
			} else{
				graph <-
					graph + geom_point(
						stroke = pointStroke,
						colour = strokeColor,
						shape = 21,
						alpha = alpha,
						size = d$sizeMultiplier
					)
			}
		}
		if (plotLabelRepel) {
			graph <- graph + geom_text_repel(color = "black", size = labelSize)
		}
		graph <- graph + xlab(variable.name)
		if (!is.null(main))
			graph <- graph + ggtitle(main)
		if (ellipse)
			graph <- graph + stat_ellipse(size = .5)
		if (plotFactorsCentroids) {
			samplePerFactor <-
				lapply(levels(colorBy), function(x)
					which(colorBy == x))
			names(samplePerFactor) <- levels(colorBy)
			centroids <-
				data.frame(t(sapply(samplePerFactor, function(x)
					colMeans(d[x, c("Axis1", "Axis2")]))), groupName = names(samplePerFactor))
			graph <-
				graph + geom_label(
					data = centroids,
					mapping = aes(
						x = .data$Axis1,
						y = .data$Axis2,
						label = .data$groupName
					),
					inherit.aes = FALSE
				)
		}
		graph <- graph + theme(
			panel.background = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.title.y = element_blank()
		) + guides(size = "none")
		if (returnGraph) {
			return(graph)
		} else{
			print(graph)
		}
	}


#' Wrapper for proj2d with highly contrasted color scale.
#'
#' @export
#' @inheritParams proj2d
#' @param ... additional arguments passed to `oob::proj2d`
#'
#' @examples
#' data(iris)
#' proj2dWideColor(iris,colorBy = iris$Petal.Length)
proj2dWideColor <-
	function(coord,
					 axis = c(1, 2),
					 colorBy = NULL,
					 pointSize = 3,
					 plotText = FALSE,
					 main = NULL,
					 alpha = 9 / 10,
					 ellipse = FALSE,
					 emph = NULL,
					 returnGraph = FALSE,
					 legendTitle = "Values",
					 axis.names = NULL,
					 na.color = "grey50",
					 na.bg = TRUE,
					 funAutoColorScale = oobColors,
					 ...) {
		proj2d(
			coord,
			axis = axis,
			colorBy = colorBy,
			pointSize = pointSize,
			plotText = plotText,
			main = main,
			alpha = alpha,
			ellipse = ellipse,
			emph = emph,
			colorScale = c("white", "yellow", "red", "purple"),
			returnGraph = returnGraph,
			legendTitle = legendTitle,
			axis.names = axis.names,
			na.color = na.color,
			na.bg = na.bg,
			funAutoColorScale = ggplotColours,
			...
		)
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
proj3d <-
	function(coord,
					 axis = c(1, 2, 3),
					 colorBy = NULL,
					 pointSize = 5,
					 plotText = FALSE,
					 colorScale = NULL,
					 na.color = "grey50",
					 na.bg = TRUE,
					 alpha = 1,
					 ...) {
		if (!(is.factor(colorBy) |
					is.numeric(colorBy) |
					is.null(colorBy)))
			stop("Error, colorBy must be numeric, factor or null.")
		if (length(axis) != 3)
			stop("You must give a vector of 3 integer for axis parameter")

		if (is.null(colorBy)) {
			colors = "black"
		}
		else{
			if (na.bg) {
				indexNA <- which(is.na(colorBy))
				indexNotNA <- which(!is.na(colorBy))
				if (length(indexNA) > 0) {
					tempc <-
						coord
					tempc[seq_along(indexNA),] <-
						coord[indexNA,]
					tempc[(length(indexNA) + 1):nrow(coord),] <-
						coord[indexNotNA,]
					coord <- tempc

					tempg <-
						colorBy
					tempg[seq_along(indexNA)] <-
						colorBy[indexNA]
					tempg[(length(indexNA) + 1):length(colorBy)] <-
						colorBy[indexNotNA]
					colorBy <- tempg

				}
			}
			if (is.factor(colorBy)) {
				if (is.null(colorScale))
					colorScale <- rainbow(nlevels(colorBy))
				hashCol <- colorScale
				names(hashCol) <- levels(colorBy)
				colors <- hashCol[colorBy]
			} else{
				if (is.null(colorScale))
					colorScale <- c("grey", "red")
				funCol <-
					circlize::colorRamp2(seq(
						min(colorBy),
						max(colorBy),
						length.out = length(colorScale)
					), colorScale)
				colors <- funCol(colorBy)
			}
		}
		colors[is.na(colors)] <- na.color
		plot3d(
			coord[, axis[1]],
			coord[, axis[2]],
			coord[, axis[3]],
			xlab = paste0("Axis", axis[1]),
			ylab = paste0("Axis", axis[2]),
			zlab = paste0("Axis", axis[3]),
			type = "n",
			...
		)
		if (plotText) {
			text3d(
				coord[, axis[1]],
				coord[, axis[2]],
				coord[, axis[3]],
				texts = rownames(coord),
				cex = pointSize,
				col = colors,
				alpha = alpha
			)
		} else{
			spheres3d(
				coord[, axis[1]],
				coord[, axis[2]],
				coord[, axis[3]],
				col = colors,
				radius = pointSize,
				alpha = alpha
			)
		}
		if (!is.null(colorBy)) {
			if (is.factor(colorBy))
				legend3d(
					"topright",
					legend = names(hashCol),
					pch = 16,
					col = hashCol,
					cex = 1,
					inset = c(0.02)
				)
			if (is.numeric(colorBy))
				legend3d(
					"topright",
					legend = c(min(colorBy), max(colorBy)),
					pch = 16,
					col = colorScale,
					cex = 1,
					inset = c(0.02)
				)
		}
	}
