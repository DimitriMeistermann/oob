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
autoGparFontSizeMatrix <- function(n, ...) {
    n <- max(n, 50)
    n <- min(n, 1000)
    return(gpar(fontsize = 1 / n * 600, ...))
}

#' Multiple ggplot one one page
#'
#' @param ... Several plot from ggplot. Each given as an argument.
#' @param plotlist Plots given as a list. Each element contains a plot from grid
#'   or inherited package.
#' @param cols Number of columns. If layout is NULL, then use 'cols' to
#'   determine layout
#' @param layout A 2D matrix of numeric value (from one the the number of plot)
#'   indicating the layout to be plotted. Matrix()
#'
#' @return NULL, print the plots in the current device
#' @export
#'
#' @examples
#' p1<-ggplot(data = data.frame(x=seq_len(5),y=seq_len(5)),
#'     aes(x=x,y=y))+geom_point()+ggtitle("p1")
#' p2<-oobqplot(c("a","a","b","b","b","c"))+ggtitle("p2")
#' p3<-oobqplot(c("a","a","b","b","b"))+ggtitle("p3")
#' p4<-oobqplot(c("a","b","c","d"))+ggtitle("p4")
#'
#' multiplot(p1,p2,p3,p4)
#' plotList<-list(p1,p2,p3,p4)
#' multiplot(plotlist = plotList)
#' multiplot(plotlist = plotList,cols = 2)
#'
#' layout<-matrix(data = c(
#'     4,3,
#'     2,1
#' ),ncol = 2,byrow = TRUE)
#' multiplot(plotlist = plotList,layout = layout)
#' layout<-matrix(data = c(
#'     1,2,3,
#'     4,0,0
#' ),ncol = 3,byrow = TRUE)
#' multiplot(plotlist = plotList,layout = layout)

multiplot <- function(...,
                    plotlist = NULL,
                    cols = 1,
                    layout = NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots <- length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                        ncol = cols,
                        nrow = ceiling(numPlots / cols))
    }
    if (numPlots == 1) {
        print(plots[[1]])
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (i in seq_len(numPlots)) {
            # Get the i,j matrix positions of
            # the regions that contain this subplot
            matchidx <-
                as.data.frame(which(layout == i, arr.ind = TRUE))

            print(
                plots[[i]],
                vp = viewport(
                    layout.pos.row = matchidx$row,
                    layout.pos.col = matchidx$col
                )
            )
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
ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
    if ((diff(h) %% 360) < 1)
        h[2] <- h[2] - 360 / n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


#' Convert color from additive to subtracting mixing
#'
#' @param color A vector of 3 numeric containing rgb values (from 0 to 1), or a
#'   single character of the color name or hex code
#' @param returnHex Logical value indicating if the color should be returned as
#'   a hex code or as a vector of rgb values.
#'
#' @return a vector containing rgb values if returnHex=FALSE, otherwise a color
#'   as a character format.
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
convertColorAdd2Sub <- function(color, returnHex = TRUE) {
    if (is.character(color)) {
        color <- col2rgb(color)[, 1] / 255
    }
    if (length(color) < 3)
        stop("if blue or green is null, the color argument must contain ",
            "3 numeric values.")
    red <- color[1]
    green <- color[2]
    blue <- color[3]

    newRed <- mean(c(1 - green, 1 - blue))
    newGreen <- mean(c(1 - red, 1 - blue))
    newBlue <- mean(c(1 - red, 1 - green))
    if (returnHex) {
        return(rgb(newRed, newGreen, newBlue))
    } else{
        return(c(
            "red" = newRed,
            "green" = newGreen,
            "blue" = newBlue
        ))
    }
}


#' Compute density value for the point of a a 2D-space
#'
#' @param mat numeric matrix of point coordinates. Each row is a point, 1st
#'   column = x, 2nd column=y.
#' @param eps radius of the eps-neighborhood, i.e., bandwidth of the uniform
#'   kernel).
#'
#' @return A vector of numeric value of the length equal to the number of points
#'   representing the density
#' @export
#'
#' @examples
#' coord<-matrix(rnorm(2000),ncol = 2)
#' pointdensity.nrd(coord)
pointdensity.nrd <- function(mat, eps = 1) {
    if (!is.matrix(mat))
        mat <- as.matrix(mat)
    dbscan::pointdensity(apply(mat, 2, function(d)
        d / max(MASS::bandwidth.nrd(d),1e-6)),
        eps = eps,
        type = "density")
}


#' log10(x+1) continuous scale for ggplot2
#'
#' @return A transformation object.
#' @export
#'
#' @examples
#' ggplot(data.frame(x=c(0,10,100,1000),y=seq_len(4)),
#'     mapping = aes(x=x,y=y))+geom_point()+
#'     scale_x_continuous(trans = log10plus1())
log10plus1 <- function() {
    scales::trans_new(
        name = "log10plus1",
        transform = function(x)
            log10(x + 1),
        inverse = function(x)
            10 ^ x - 1 ,
        domain = c(0, Inf)
    )
}

#' Plot colors
#'
#' @param colorScale Character vector, Colors to plot.
#' @param continuousStep Number of plotted color with intermediate between the
#'   given colors.
#'
#' @return NULL, draw in current graphic device.
#' @export
#'
#' @examples
#' plotPalette(c("#001D6E","white","#E6001F"))
#' plotPalette(c("#001D6E","white","#E6001F"),continuousStep = 50)
plotPalette <- function(colorScale, continuousStep = NULL) {
    if (is.null(continuousStep)) {
        image(
            seq_along(colorScale),
            1,
            as.matrix(seq_along(colorScale)),
            col = colorScale,
            xlab = "",
            ylab = "",
            xaxt = "n",
            yaxt = "n",
            bty = "n"
        )
        if (!is.null(names(colorScale)))
            axis(1, seq_along(colorScale), names(colorScale))
    } else{
        br <- round(seq(1, continuousStep, length.out = length(colorScale)))
        cols <- colorRamp2(breaks = br,
                        colors = colorScale)(seq_len(continuousStep))
        image(
            seq_len(continuousStep),
            1,
            as.matrix(seq_len(continuousStep)),
            col = cols,
            xlab = "",
            ylab = "",
            xaxt = "n",
            yaxt = "n",
            bty = "n"
        )
    }
}


#' Plot a double arrow in a grid plot.
#'
#' @param x Numeric. X coordinate of the arrow.
#' @param y Numeric. Y coordinate of the arrow.
#' @param width Numeric. X coordinate of the arrow.
#' @param height Numeric. X coordinate of the arrow.
#' @param just A string or numeric vector specifying the justification of the
#'   viewport relative to its (x, y) location. If there are two values, the
#'   first value specifies horizontal justification and the second value
#'   specifies vertical justification. Possible string values are: "left",
#'   "right", "centre", "center", "bottom", and "top". For numeric values, 0
#'   means left alignment and 1 means right alignment.
#' @param gp An object of class "gpar", typically the output from a call to the
#'   function gpar. This is basically a list of graphical parameter settings
#' @param ... Other parameters passed to pushViewport.
#'
#' @return Plot in the current graphical device.
#' @export
#'
#' @examples
#' filledDoubleArrow()
#' oobqplot(seq_len(5))
#' filledDoubleArrow()
filledDoubleArrow <-
    function(x = 0,
            y = 0,
            width = 1,
            height = 1,
            just = c("left", "bottom"),
            gp = gpar(col = "black"),
            ...) {
        pushViewport(viewport(
            x = x,
            y = y,
            width = width,
            height = height,
            just = just,
            ...
        ))
        grid.polygon(
            x = c(0, .2, .2, .8, .8, 1, .8, .8, .2, .2, 0),
            y = c(.5, .7, .6, .6, .7, .5, .3, .4, .4, .3, .5),
            gp = gp
        )
        popViewport()
    }


#' Plot the expression of one or several genes
#'
#' @param expr Numeric 2D matrix. Each row is gene and each column a sample. Row
#'   must be named by genes.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param group Factor vector, same length as number of column in expr.
#'   Experimental group attributed to each sample.
#' @param log10Plus1yScale Logical or NULL. Use a `log10(x+1)` scale. By default
#'   `FALSE` if one gene and `TRUE` if several.
#' @param violin Logical. Plot violin plots.
#' @param boxplot Logical. Plot violin plots.
#' @param dotplot Logical. Plot dot plots.
#' @param violinArgs List. Additional arguments given to `geom_violin`.
#' @param boxplotArgs List. Additional arguments given to `geom_boxplot`.
#' @param dotplotArgs List. Additional arguments given to `geom_beeswarm`.
#' @param colorScale A list of color. Must be the same length as number of
#'   levels in group.
#' @param legendTitle Character. Displayed in legend title.
#' @param dodge.width Numeric. Width of individual distribution element
#'   (violin/boxplot/dotplot).
#' @param returnGraph Logical. Print the ggplot object or return it.
#' @param sce_assay Integer or character,
#'   if `data` is a `SummarizedExperiment` object, the assay name to use.
#'
#' @return Plot in the current graphical device or a ggplot object if
#'   `returnGraph=TRUE`.
#' @export
#'
#' @examples
#' exprMat<-rbind(
#'     MASS::rnegbin(30,theta = 10,mu = 10),
#'     MASS::rnegbin(30,theta = 1,mu = 10),
#'     MASS::rnegbin(30,theta = 1,mu = 1000)
#' )
#' rownames(exprMat)<-c("gene1","gene2","gene3")
#'
#' group=c(rep("A",10),rep("B",10),rep("C",10))
#'
#' plotExpr(exprMat)
#' plotExpr(exprMat["gene1",])
#' plotExpr(exprMat["gene1",],group = group)
#' plotExpr(exprMat,group = group)
#'
#' plotExpr(exprMat,group = group,log10Plus1yScale = FALSE)
#' plotExpr(exprMat,group = group,dodge.width = .5,legendTitle = "Letter")
#' plotExpr(exprMat,boxplot = FALSE,violin = FALSE,dotplot = TRUE)
#' plotExpr(exprMat,group = group,colorScale = c("blue","white","red"))
#' plotExpr(exprMat,group = group,returnGraph = TRUE)+
#'     geom_point()
#'
#' plotExpr(exprMat,group = group,violinArgs = list(scale="area"))
#' sce <- SingleCellExperiment(assays = list(counts = exprMat))
#' plotExpr(sce["gene1",])

plotExpr <-
    function(expr,
            group = NULL,
            log10Plus1yScale = NULL,
            violin = TRUE,
            boxplot = TRUE,
            dotplot = FALSE,
            violinArgs = list(),
            boxplotArgs = list(),
            dotplotArgs = list(),
            colorScale = oobColors,
            legendTitle = "group",
            dodge.width = .9,
            returnGraph = FALSE,
            sce_assay = 1) {

        if (inherits(expr, "SummarizedExperiment")) {
            expr <- assay(expr, sce_assay)
        }
        barplotGraph <- greyGraph <- coloredGraph <- FALSE

        if (is.vector(expr))
            expr <- t(as.matrix(expr))
        if (is.vector(group) | is.factor(group)) {
            group <- data.frame(group = group, stringsAsFactors = TRUE)
            rownames(group) <- colnames(expr)
        }

        if (!is.matrix(expr))
            expr <- as.matrix(expr)
        if (is.null(log10Plus1yScale))
            log10Plus1yScale <-
                nrow(expr) > 1
        #if more than one gene log10Plus1yScale is turned on

        if (is.null(group)) {
            if (nrow(expr) < 2) {
                if (is.null(colnames(expr)))
                    colnames(expr) <-
                        as.character(seq_len(ncol(expr)))
                ggData <-
                    data.frame(
                        expression = expr[1,],
                        sample = factor(
                            colnames(expr), levels =
                                colnames(expr)[order(expr[1,],
                                    decreasing = TRUE)]
                        )
                    )
                ggData$sample <- factor(ggData$sample)
                g <-
                    ggplot(ggData, mapping =
                            aes(x = .data$sample, y = .data$expression)) +
                    geom_bar(stat = "identity") +
                    theme_bw() + theme(axis.text.x = element_text(
                        angle = 90,
                        hjust = 1,
                        vjust = .3
                    ))
                barplotGraph <- TRUE
            } else{
                ggData <-
                    reshape2::melt(expr,
                                value.name = "expression",
                                varnames = c("gene", "sample"))
                ggData$gene <-
                    factor(ggData$gene, levels = rownames(expr))
                g <-
                    ggplot(ggData, mapping = aes(
                        x = .data$gene,
                        y = .data$expression)
                        ) +
                    theme_bw() + theme(axis.text.x = element_text(
                        angle = 90,
                        hjust = 1,
                        vjust = .3,
                        face = "bold.italic"
                    ))
                greyGraph <- TRUE
            }
        } else{
            #group on
            if (ncol(expr) != nrow(group))
                stop("expr should have the same number of sample than group")
            if (ncol(group) > 1)
                stop("Multiple group are not allowed in the same time")
            groupName <- colnames(group)
            if (nrow(expr) == 1) {
                ggData <- data.frame(expression = expr[1,], group)
                g <-
                    ggplot(ggData, mapping = aes(x = .data[[groupName]],
                                            y = .data$expression)) +
                    theme_bw() + theme(axis.text.x = element_text(
                        angle = 90,
                        hjust = 1,
                        vjust = .3,
                        face = "bold"
                    ))
                greyGraph <- TRUE
            } else{
                coloredGraph <- TRUE
                ggData <-
                    reshape2::melt(
                        data.frame(t(expr), group),
                        value.name = "expression",
                        variable.name = "gene",
                        id.vars = groupName
                    )
                # Used below (both in and out of the violin branch) to filter `colors` --
                # default to "drop nothing" so non-violin plots (boxplot/dotplot, which have
                # no n<3 incompatibility) aren't affected by this violin-specific rule.
                factor2drop <- character(0)
                if (violin) {
                    factorSampling <- table(group[, 1])
                    factor2drop <-
                        names(factorSampling)[factorSampling < 3]
                    if (length(factor2drop) > 1)
                        warning(
                            paste0(factor2drop, collapse = " "),
                            " were dropped (n<3 is not compatible with violin",
                            " plot). You can deactivate violin layer by ",
                            "setting violin argument to FALSE"
                        )
                    ggData <-
                        ggData[!ggData[, groupName] %in% factor2drop,]
                    #drop levels where n < 3
                }
                colors <-
                    if (is.function(colorScale))
                        colorScale(nlevels(group[, 1]))
                else
                    colorScale
                g <-
                    ggplot(ggData,
                        mapping = aes(
                            x = .data$gene,
                            y = .data$expression,
                            fill = .data[[groupName]]
                        )) +
                    theme_bw() + theme(axis.text.x = element_text(
                        angle = 90,
                        hjust = 1,
                        vjust = .3,
                        face = "bold.italic"
                    )) +
                    scale_fill_manual(
                        values = colors[!levels(group[, 1]) %in% factor2drop],
                        name =legendTitle
                    )
            }
        }
        if (greyGraph) {
            if (is.null(violinArgs$fill))
                violinArgs$fill <- "grey50"
            if (is.null(violinArgs$scale))
                violinArgs$scale <- "width"
            if (is.null(boxplotArgs$width))
                boxplotArgs$width <- .2
        }
        if (coloredGraph) {
            if (is.null(violinArgs$scale))
                violinArgs$scale <- "width"
            if (is.null(boxplotArgs$width))
                boxplotArgs$width <- .2
            if (is.null(violinArgs$position))
                violinArgs$position <-
                    position_dodge(preserve = "total", width = dodge.width)
            if (is.null(boxplotArgs$position))
                boxplotArgs$position <-
                    position_dodge(preserve = "total", width = dodge.width)
            if (is.null(dotplotArgs$dodge.width))
                dotplotArgs$dodge.width <- dodge.width
        }
        if (!barplotGraph) {
            g <- ggBorderedFactors(g,
                        borderColor = "black",
                        borderSize = .5)
            if (violin)
                g <- g + do.call("geom_violin", violinArgs)
            if (boxplot)
                g <- g + do.call("geom_boxplot", boxplotArgs)
            if (dotplot)
                g <- g + do.call("geom_beeswarm", dotplotArgs)
        }
        if (log10Plus1yScale) {
            maxExpr <- max(expr)
            ncharMaxExpr <- nchar(round(maxExpr))
            breaks <- c(0, 2, round(10 ^ (seq(
                1, ncharMaxExpr, 0.5
            ))))
            #breaks<-c(0,rbind(breaks/2,breaks))
            #intelacing 1,10,100... and 5,50,500...
            if (maxExpr < breaks[length(breaks) - 1])
                breaks <- breaks[seq_along(breaks) - 1]
            g <-
                g + scale_y_continuous(
                    trans = log10plus1(),
                    limits = c(breaks[1], breaks[length(breaks)]),
                    breaks = breaks,
                    minor_breaks = NULL
                )
        }
        if (!returnGraph) {
            print(g)
        } else{
            return(g)
        }
    }


#' Volcano plot with additional annotation for interpreting DE genes.
#'
#' @param DEresult Dataframe that contains at least those columns:
#' - padj (adjusted p-value)
#' - isDE (a character vector equal to "NONE" if the gene is not DE,
#'  "DOWNREG" or "UPREG" if DE).
#' - log2FoldChange
#'   Row must be named by genes.
#' @param formula Character. Design formula given to DESeq2.
#' @param downLevel Character. Condition considered as the reference. If a gene
#'   is more expressed in this condition, LFC < 0.
#' @param upLevel Character. Condition considered as the target group. If a gene
#'   is more expressed in this condition, LFC > 0.
#' @param condColumn Character. Name of the experimental variable that have been
#'   used for differential expression.
#' @param padjThreshold Numeric. Significance threshold of the adjusted p-value.
#' @param LFCthreshold Numeric. Significance threshold of the Log2 Fold-Change.
#' @param topGene Integer. Number of gene name to be shown on the plot. Genes
#'   names are plotted from the most significant.
#'
#' @return Plot in the current graphical device.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' volcanoPlot.DESeq2(DEgenesPrime_Naive,formula = "~culture_media+Run",
#'     condColumn = "culture_media",downLevel = "KSR+FGF2",upLevel = "T2iLGO")
volcanoPlot.DESeq2 <-
    function(DEresult,
            formula,
            downLevel,
            upLevel,
            condColumn,
            padjThreshold = 0.05,
            LFCthreshold = 1,
            topGene = 30) {
        DEresult <- DEresult[!is.na(DEresult$padj),]
        gene2Plot <- order(DEresult$padj)
        gene2Plot <-
            gene2Plot[DEresult[gene2Plot, "isDE"] != "NONE"]
        gene2Plot <-
            gene2Plot[seq_len(min(topGene, length(gene2Plot)))]
        g <-
            ggplot(DEresult,
                aes(
                    x = .data$log2FoldChange,
                    y = -log10(.data$padj),
                    color = .data$isDE
                )) +
            geom_point(size = 1) + theme_bw() +
            scale_color_manual(values = c("#3AAA35", "grey75", "#E40429")) +
            geom_text_repel(
                data = DEresult[gene2Plot,],
                aes(
                    x = .data$log2FoldChange,
                    y = -log10(.data$padj),
                    label = rownames(DEresult)[gene2Plot]
                ),
                inherit.aes = FALSE,
                color = "black",
                fontface = "bold.italic",
                size = 3
            ) +
            ylab("-log10(adjusted pvalue)") + xlab(NULL) +
            geom_vline(xintercept = c(-LFCthreshold, LFCthreshold)) +
            geom_hline(yintercept = -log10(padjThreshold)) +
            guides(color = "none") +
            ggtitle("Volcano plot")

        grid.newpage()

        pushViewport(viewport(
            x = 0,
            y = 0,
            width = .8,
            height = .1,
            just = c("left", "bottom")
        ))
        filledDoubleArrow(
            x = .3,
            y = 1,
            width = .3,
            just = c("left", "center"),
            gp = gpar(fill = "black")
        )
        grid.text(
            label = downLevel,
            x = .28,
            y = 1.03,
            just = c("right", "center"),
            gp = gpar(fontface = "bold")
        )
        grid.text(
            label = upLevel,
            x = .62,
            y = 1.03,
            just = c("left", "center"),
            gp = gpar(fontface = "bold")
        )
        grid.text(
            label = "log2(Fold-Change)",
            x = .45,
            y = 1.5,
            just = c("center", "center")
        )
        popViewport()
        pushViewport(viewport(
            x = .65,
            y = 1,
            width = .35,
            height = 1,
            just = c("left", "top"),
            default.units = "npc"
        ))
        grid.text(
            label = paste0("Experimental design:\n", formula),
            x = 0,
            y = .95,
            just = c("left", "center")
        )
        grid.text(
            label = paste0("Results for\n", condColumn, ":\n",
                        downLevel, " vs ", upLevel),
            x = .0,
            y = .8,
            just = c("left", "center")
        )
        grid.text(
            label = paste0(sum(DEresult$isDE == "DOWNREG"), " downreg. genes"),
            x = 0,
            y = .65,
            just = c("left", "center"),
            gp = gpar(col = "#3AAA35", fontface = "bold")
        )
        grid.text(
            label = paste0(sum(DEresult$isDE == "UPREG"), " upreg. genes"),
            x = .0,
            y = .55,
            just = c("left", "center"),
            gp = gpar(col = "#E40429", fontface = "bold")
        )
        grid.text(
            label = paste0("From ", nrow(DEresult), "\ntested genes"),
            x = 0,
            y = .45,
            just = c("left", "center")
        )
        popViewport()
        main_vp <- viewport(
            x = 0,
            y = 1,
            width = .8,
            height = .9,
            just = c("left", "top")
        )
        pushViewport(main_vp)
        print(g, vp = main_vp)
        popViewport()
    }


#' Colors for a qualitative scale
#'
#' @param n Number of wanted colors
#'
#' @return A vector of colors
#' @export
#'
#' @examples
#' oobColors() |> plotPalette()
#' oobColors(n=5) |> plotPalette()
#' oobColors(n=40) |> plotPalette()
oobColors <- function(n = 20) {
	myCOlors <- c(
		"#E52421",
		"#66B32E",
		"#2A4B9B",
		"#6EC6D9",
		"#F3E600",
		"#A6529A",
		"#7C1623",
		"#006633",
		"#29235C",
		"#0084BC",
		"#E6007E",
		"#F49600",
		"#E3E3E3",
		"#626F72",
		"#040505",
		"#E74B65",
		"#95B37F",
		"#683C11",
		"#F8BAA0",
		"#DD8144"
	)
	if (n <= 20) {
		return(myCOlors[seq_len(n)])
	} else{
		return(
			extendColorPalette(
				n = n,
				colors = myCOlors,
				sortColorIn = TRUE,
				sortColorOut = TRUE
			)
		)
	}
}



#' Complex Heatmap wrapper optimized for RNA-Seq analyses...
#'
#' @inheritParams ComplexHeatmap::Heatmap
#' @param matrix A matrix. Either numeric or character.
#'   If it is a simple vector, it will be converted to a one-column matrix.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param preSet A value from `"expr"`, `"cor"`, `"dist"` or `NULL`. Change
#'   other arguments given a specific preset (default preSet if NULL).
#' @param autoFontSizeRow Logical, should row names font size automatically
#'   adjusted to the number of row?
#' @param autoFontSizeColumn Logical, should column names font size
#'   automatically adjusted to the number of columns?
#' @param scale Logical. Divide rows of `matrix` by their standard deviation. If
#'   NULL determined by preSet.
#' @param center Logical. Subtract rows of `matrix` by their average. If NULL
#'   determined by preSet.
#' @param returnHeatmap Logical, return the plot as a Heatmap object or print it
#'   in the current graphical device.
#' @param additionnalRowNamesGpar List. Additional parameter passed to `gpar`
#'   for row names.
#' @param additionnalColNamesGpar List. Additional parameter passed to `gpar`
#'   for column names.
#' @param border Logical. Whether draw border. The value can be logical or a
#'   string of color.
#' @param colorScale A vector of colors that will be used for mapping colors to
#'   the main heatmap.
#' @param colorScaleFun A function that map values to colors. Used for the main
#'   heatmap. If not NULL this will supersede the use of the `colorScale`
#'   argument.
#' @param midColorIs0 Logical. Force that 0 is the midColor.  If NULL turned on
#'   if the matr.
#' @param probs A numeric vector (between 0 and 1) same length as color or NULL.
#'   Quantile probability of the values that will be mapped to colors.
#' @param useProb Logical. Use quantile probability to map the colors. Else the
#'   min and max of values will be mapped to first and last color and
#'   interpolated continuously.
#' @param minProb A numeric value (between 0 and 1). If `useProb=TRUE` and
#'   `probs=NULL` this will be the quantile of the value for the first color,
#'   quantile will be mapped continuously as to the maxProb.
#' @param maxProb A numeric value (between 0 and 1).
#' @param colData A vector of factor, character, numeric or logical. Or, a
#'   dataframe of any of these type of value. The annotation that will be
#'   displayed on the heatmap.
#' @param colorAnnot List or NULL. Precomputed color scales for the `colData`.
#'   Color scales will be only generated for the features not described. Must be
#'   in the format of a list named by columns of `annots`. Each element contains
#'   the colors at breaks for continuous values. In the case of factors, the
#'   colors are named to their corresponding level or in the order of the
#'   levels.
#' @param showGrid Logical. Draw a border of each individual square on the
#'   heatmap. If NULL automatically true if number of values < 500.
#' @param gparGrid Gpar object of the heatmap grid if `showGrid`.
#' @param showValues Logical. Show values from the matrix in the middle of each
#'   square of the heatmap.
#' @param Nsignif Integer. Number of significant digits showed if `showValues`.
#' @param squareHt Logical or NULL. Apply clustering columns on rows. If NULL
#'   automatically turned TRUE if `ncol==nrow` and col/rownames are the same.
#' @param sce_assay Integer or character, if `data`
#'   is a `SummarizedExperiment` related object, the assay name to use.
#' @param ... Other parameters passed to `Heatmap`.
#'
#' @return A Heatmap object if `returnHeatmap` or print the Heatmap in the
#'   current graphical device.
#' @export
#'
#' @seealso [genTopAnnot()], [genRowAnnot()]
#'
#' @details
#'
#' A preSet attributes a list of default values for each argument. However, even
#' if a preSet is selected, arguments precised by the user precede the preSet. #
#' Default arguments ## preSet is `NULL`
#' ```
#' clustering_distance_rows = covDist #see covDist for more details
#' clustering_distance_columns = covDist
#' name="matrix"
#' colorScale=c("#2E3672","#4B9AD5","white","#FAB517","#E5261D")
#' center=TRUE
#' scale=FALSE
#' ```
#' ## preSet is `"expr"` (expression)
#' ```
#' clustering_distance_rows = covDist
#' clustering_distance_columns = covDist
#' name="centered log expression"
#' colorScale=    c("darkblue","white","red2")
#' additionnalRowNamesGpar=list(fontface="italic")
#' center=TRUE
#' scale=FALSE
#' ```
#' ## preSet is `"cor"` (correlation)
#' ```
#' clustering_distance_rows ="euclidean"
#' clustering_distance_columns ="euclidean"
#' name="Pearson correlation"
#' colorScale=c("darkblue","white","#FFAA00")
#' center=FALSE
#' scale=FALSE
#' ```
#' ## preSet is `"dist"` (distance)
#' ```
#' clustering_distance_rows ="euclidean"
#' name="Euclidean distance"
#' colorScale=c("white","yellow","red","purple")
#' center=FALSE
#' scale=FALSE
#' ```
#'
#' ## preSet is `"vanilla"` (don't transform value, same as default
#' ComplexHeatmap)
#' ```
#' clustering_distance_rows ="euclidean"
#' name="matrix"
#' colorScale=c("#2E3672","#4B9AD5","white","#FAB517","#E5261D")
#' center=FALSE
#' scale=FALSE
#' ```
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("DEgenesPrime_Naive")
#'
#' bestDE <- rownames(DEgenesPrime_Naive)[whichTop(DEgenesPrime_Naive$pvalue,
#'                                           decreasing = FALSE,
#'                                           top = 50)]
#' heatmap.DM(
#'     matrix(rnorm(50), ncol = 5),
#'     preSet = NULL,
#'     showValues = TRUE,
#'     Nsignif = 2
#' )
#'
#' heatmap.DM(bulkLogCounts[bestDE, ],
#'   colData = sampleAnnot[, c("culture_media", "line")])
#' heatmap.DM(
#'   bulkLogCounts[
#'     bestDE[seq_len(5)],
#'     rownames(sampleAnnot)[sampleAnnot$culture_media %in%
#'       c("T2iLGO","KSR+FGF2")],
#'   ]
#' )
#'
#' corDat <- cor(bulkLogCounts)
#' heatmap.DM(corDat, preSet = "cor")
#' heatmap.DM(
#'     corDat,
#'     preSet = "cor",
#'     center = TRUE,
#'     colorScaleFun = circlize::colorRamp2(c(-0.2, 0, 0.2),
#'       c("blue", "white", "red"))
#' )
#' sce <- SingleCellExperiment(assays = list(counts = bulkLogCounts),
#'    colData = sampleAnnot)
#' heatmap.DM(sce[bestDE[seq_len(5)],], colData = c("line", "culture_media"))
heatmap.DM <-
    function(matrix,
            preSet = "expr",
            clustering_distance_rows = NULL,
            clustering_distance_columns = NULL,
            clustering_method_columns = "ward.D2",
            clustering_method_rows = "ward.D2",
            autoFontSizeRow = TRUE,
            autoFontSizeColumn = TRUE,
            scale = NULL,
            center = NULL,
            returnHeatmap = FALSE,
            name = NULL,
            additionnalRowNamesGpar = NULL,
            additionnalColNamesGpar = list(),
            border = TRUE,
            colorScale = NULL,
            colorScaleFun = NULL,
            midColorIs0 = NULL,
            probs = NULL,
            useProb = TRUE,
            minProb = 0.05,
            maxProb = 0.95,
            cluster_rows = NULL,
            cluster_columns = NULL,
            colData = NULL,
            colorAnnot = NULL,
            showGrid = NULL,
            gparGrid = gpar(col = "black"),
            showValues = FALSE,
            Nsignif = 3,
            column_dend_reorder = FALSE,
            row_dend_reorder = FALSE,
            squareHt = NULL,
            row_split = NULL,
            column_split = NULL,
            sce_assay = "logcounts",
            ...) {

    if (inherits(matrix, "SummarizedExperiment")) {
        if (!is.null(colData)) {
            if(inherits(colData, "character")){
            	if (sum(!colData %in% colnames(colData(matrix))==0)) {
            		colData <- data.frame(colData(matrix)[colData])
            	} else {
            		stop(
            			setdiff(colData, colnames(colData(matrix))),
            			"not found in colData of SingleCellExperiment"
            		)
          	  }
            }
        }
        matrix <- assay(matrix, sce_assay)
    }
    args <- list()

    if (is.null(preSet)) {
        if (is.null(clustering_distance_rows))
            clustering_distance_rows <- covDist
        if (is.null(clustering_distance_columns))
            clustering_distance_columns <- covDist
        if (is.null(name))
            name <- "matrix"
        if (is.null(colorScale))
            colorScale <- c("#2E3672", "#4B9AD5", "white", "#FAB517", "#E5261D")
        if (is.null(additionnalRowNamesGpar))
            additionnalRowNamesGpar <- list()
        if (is.null(center))
            center <- TRUE
        if (is.null(scale))
            scale <- FALSE
    } else if (preSet == "expr") {
        if (is.null(clustering_distance_rows))
            clustering_distance_rows <- covDist
        if (is.null(clustering_distance_columns))
            clustering_distance_columns <- covDist
        if (is.null(name))
            name <- "centered log expression"
        if (is.null(colorScale))
            colorScale <- c("darkblue", "white", "red2")
        if (is.null(additionnalRowNamesGpar))
            additionnalRowNamesGpar <- list(fontface = "italic")
        if (is.null(center))
            center <- TRUE
        if (is.null(scale))
            scale <- FALSE
    } else if (preSet == "cor") {
        if (is.null(clustering_distance_rows))
            clustering_distance_rows <- "euclidean"
        if (is.null(clustering_distance_columns))
            clustering_distance_columns <- "euclidean"
        if (is.null(name))
            name <- "Pearson\ncorrelation"
        if (is.null(colorScale))
            colorScale <- c("darkblue", "white", "#FFAA00")
        if (is.null(additionnalRowNamesGpar))
            additionnalRowNamesGpar <- list()
        if (is.null(center))
            center <- FALSE
        if (is.null(scale))
            scale <- FALSE
    } else if (preSet == "dist") {
        if (is.null(clustering_distance_rows))
            clustering_distance_rows <- "euclidean"
        if (is.null(clustering_distance_columns))
            clustering_distance_columns <- "euclidean"
        if (is.null(name))
            name <- "Euclidean\ndistance"
        if (is.null(colorScale))
            colorScale <- c("white", "yellow", "red", "purple")
        if (is.null(additionnalRowNamesGpar))
            additionnalRowNamesGpar <- list()
        if (is.null(center))
            center <- FALSE
        if (is.null(scale))
            scale <- FALSE
    } else if (preSet == "vanilla") {
        if (is.null(clustering_distance_rows))
            clustering_distance_rows <- "euclidean"
        if (is.null(clustering_distance_columns))
            clustering_distance_columns <- "euclidean"
        if (is.null(name))
            name <- "matrix"
        if (is.null(colorScale))
            colorScale <- c("#2E3672", "#4B9AD5", "white", "#FAB517", "#E5261D")
        if (is.null(additionnalRowNamesGpar))
            additionnalRowNamesGpar <- list()
        if (is.null(center))
            center <- FALSE
        if (is.null(scale))
            scale <- FALSE
    } else{
        stop("preSet must equal to one of this value: NULL, ",
            "'expr', 'cor', 'dist', 'vanilla'")
    }
    matrix <- as.matrix(matrix)

    if (min(apply(matrix, 1, sd, na.rm = TRUE)) == 0 &
        (scale |
        identical(corrDist, clustering_distance_rows))) {
        warning(
            "some row have a 0 sd. sd-based method ",
            "(correlation distance, scaling) ",
            "will be deactivated or switched."
        )
        scale <- FALSE
        if (identical(corrDist, clustering_distance_rows)) {
            args$clustering_distance_rows <- "euclidean"
        }
    }
    if (scale |
        center)
        matrix <-
        rowScale(matrix, scaled = scale, center = center)
    if (is.null(midColorIs0)) {
        if (min(matrix, na.rm = TRUE) < 0 & max(matrix, na.rm = TRUE) > 0) {
            midColorIs0 <- TRUE
        } else{
            midColorIs0 <- FALSE
        }
    }
    if (is.null(squareHt)) {
        if (nrow(matrix) == ncol(matrix) &
            identical(colnames(matrix),rownames(matrix))) {
            squareHt <- TRUE
            warning("colnames and rownames are identical, ",
                    "squareHt is set to TRUE")
        } else{
            squareHt <- FALSE
        }
    }
    if (squareHt) {
        if (is.null(cluster_columns)) {
            cluster_columns <-
                hierarchicalClustering(
                    matrix,
                    transpose = FALSE,
                    method.dist = clustering_distance_columns,
                    method.hclust = clustering_method_columns
                )
        }
        args$cluster_rows <- cluster_columns
        args$cluster_columns <- cluster_columns
    } else{
        if (is.null(cluster_rows)) {
            args$clustering_method_rows <- clustering_method_rows
            args$clustering_distance_rows <-
                clustering_distance_rows
        } else{
            args$cluster_rows <- cluster_rows
        }
        if (is.null(cluster_columns)) {
            args$clustering_method_columns <- clustering_method_columns
            args$clustering_distance_columns <-
                clustering_distance_columns
        } else{
            args$cluster_columns <- cluster_columns
        }
    }

    if (is.null(colorScaleFun)) {
        colorScaleFun <-
            computeColorScaleFun(
                colors = colorScale,
                values = unlist(matrix),
                useProb = useProb,
                probs = probs,
                minProb = minProb,
                maxProb = maxProb,
                midColorIs0 = midColorIs0,
                returnColorFun = TRUE
            )
    }
    args$col <- colorScaleFun
    if (is.null(showGrid)) {
        if (nrow(matrix) * ncol(matrix) < 500) {
            showGrid <- TRUE
        } else{
            showGrid <- FALSE
        }
    }
    if (showGrid) {
        args$rect_gp <- gparGrid
    }
    if (showValues) {
        args$cell_fun <- function(j, i, x, y, w, h, col) {
            #dark or light background .
            if (colSums(col2rgb(col)) < 382.5)
                col <- "white"
            else
                col <- "black"
            grid.text(
                as.character(
                    signif(matrix[i, j], Nsignif)
                ), x, y, gp = gpar(col = col)
            )
        }
    }
    if (autoFontSizeRow)
        args$row_names_gp <- do.call("autoGparFontSizeMatrix",
            c(list(nrow(matrix)), additionnalRowNamesGpar))
    if (autoFontSizeColumn)
        args$column_names_gp <- do.call("autoGparFontSizeMatrix",
            c(list(ncol(matrix)), additionnalColNamesGpar))

    if (!is.null(colData)) {
        args$top_annotation <- genTopAnnot(colData, colorAnnot)
    }

    args$column_dend_reorder <-
        column_dend_reorder
    args$row_dend_reorder <- row_dend_reorder
    args$row_split <- row_split
    args$column_split <- column_split
    args$matrix <- matrix
    args$name <- name
    args$border <- border
    args <- c(args, list(...))

    ht <- do.call("Heatmap", args)
    if (returnHeatmap) {
        return(ht)
    } else{
        print(ht)
    }
}
