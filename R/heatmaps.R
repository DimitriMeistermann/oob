
#' Complex Heatmap wrapper optimized for RNA-Seq analyses...
#'
#' @inheritParams ComplexHeatmap::Heatmap
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
#' library(ComplexHeatmap)
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
            ...) {
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

    if (min(apply(matrix, 1, sd)) == 0 &
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
        if (min(matrix) < 0 & max(matrix) > 0) {
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




#' Plot best marker per group on a tidy Heatmap
#'
#' @param countTable A matrix of numeric with samples as columns (in the RNA-Seq
#'   context, log counts)
#' @param group A feature of factor/character, same length as number of sample.
#'   Describe group of each sample (for example clusters).
#' @param markerData A matrix describing marker scores for each group.
#' @param topn Number of marker to plot, rankes from the best
#' @param show_column_names Whether show column names.
#' @param ... Arguments passed to `heatmap.DM`
#'
#' @inheritParams heatmap.DM
#' @inheritParams drawSamplePerGroup
#'
#' @details Draw the same number of observation from each condition/group, and
#'   take the top n marker per group. Heatmap is sliced by group for each gene
#'   and
#'
#' @return A Heatmap object if `returnHeatmap` or print the Heatmap in the
#'   current graphical device.
#' @export
#'
#' @seealso [heatmap.DM()]
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media)
#' htMarker(bulkLogCounts,  group=sampleAnnot$culture_media,
#'     markerData=extractFeatureMarkerData(markerData),
#'     colData=sampleAnnot[c("line","passage")])

htMarker <-
    function(countTable,
            group,
            markerData,
            colData = NULL,
            topn = 5,
            maxDrawSize = NULL,
            minDrawSize = NULL,
            replace = FALSE,
            returnHeatmap = FALSE,
            show_column_names = FALSE,
            ...)
    {
        if (length(group) != ncol(countTable))
            stop("Length of group should be the same as number of columns in ",
                "countTable")
        group <- as.factor(make.names(group))
        if (sum(!sort(colnames(markerData)) == sort(levels(group))) >
            0)
            stop("colnames of markerData must ",
                "correspond to the levels of group")
        if(is.null(colnames(countTable)))
            colnames(countTable) <- make.names(
                paste0("X",seq_len(ncol(countTable))))
        names(group) <- colnames(countTable)
        if(!all(rownames(markerData) %in% rownames(countTable))) {
            stop("rownames of markerData must ",
                "exist in the rownames of countTable")
        }

        drawPerGroup<-drawSamplePerGroup(names(group),group,
            minDrawSize = maxDrawSize,
            maxDrawSize = maxDrawSize, replace = replace)
        drawCells <- names(drawPerGroup)

        topMarker <- apply(markerData, 2, function(x) {
            names(x)[order(x, decreasing = TRUE)][seq_len(topn)]
        })
        topMarker <- as.list(data.frame(topMarker))
        heatmap.DM(
            countTable[unlist(topMarker), drawCells],
            colData = colData[drawCells,],
            column_split = drawPerGroup,
            row_split = VectorListToFactor(topMarker),
            cluster_row_slices = FALSE,
            cluster_column_slices = FALSE,
            returnHeatmap = returnHeatmap,
            show_column_names = show_column_names,
            ...
        )
    }


#' Generate a top annotation for ComplexHeatmap
#'
#' @param annot A vector of factor, character, numeric or logical. Or, a
#'   dataframe of any of these type of value. The annotation that will be
#'   displayed on the heatmap.
#' @param colorScales List or NULL. Precomputed color scales. Color scales will
#'   be only generated for the features not described. Must be in the format of
#'   a list named by columns of `annots`. Each element contains the colors at
#'   breaks for continuous values. In the case of factors, the colors are named
#'   to their corresponding level or in the order of the levels.
#' @param border Logical. Whether draw border. The value can be logical or a
#'   string of color.
#' @param ... Other parameters passed to `genColorsForAnnots`.
#'
#' @return A HeatmapAnnotation object. Can be used for example in the
#'   `top_annotation` argument of `Heatmap`.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("DEgenesPrime_Naive")
#'
#' bestDE <-
#'     rownames(DEgenesPrime_Naive)[
#'         whichTop(DEgenesPrime_Naive$pvalue,
#'         decreasing = FALSE,
#'         top = 50)
#'     ]
#' ComplexHeatmap::Heatmap(rowScale(bulkLogCounts[bestDE, ]), top_annotation =
#'     genTopAnnot(sampleAnnot[, c("culture_media", "line")]))
genTopAnnot <- function(annot,
                        colorScales = NULL,
                        border = TRUE,
                        ...) {
    if (is.factor(annot) | is.data.frame(annot))
        annot <- droplevels(annot)
    if (is.vector(annot) | is.factor(annot)) {
        annot <- data.frame(Annotation = annot)
        if (is.list(colorScales))
            colnames(annot) <- names(colorScales)[1]
        if ((!is.null(colorScales)) &
            !is.list(colorScales))
            colorScales <- list("Annotation" = colorScales)
    }
    colorScales <-
        genColorsForAnnots(
            annots = annot,
            colorScales = colorScales,
            returnContinuousFun = TRUE,
            ...
        )
    HeatmapAnnotation(df = annot,
                    col = colorScales,
                    border = border)
}


#' Generate a row annotation for ComplexHeatmap
#'
#' @param annot A vector of factor, character, numeric or logical. Or, a
#'   dataframe of any of these type of value. The annotation that will be
#'   displayed on the heatmap.
#' @param colorScales List or NULL. Precomputed color scales. Color scales will
#'   be only generated for the features not described. Must be in the format of
#'   a list named by columns of `annots`. Each element contains the colors at
#'   breaks for continuous values. In the case of factors, the colors are named
#'   to their corresponding level or in the order of the levels.
#' @param border Logical. Whether draw border. The value can be logical or a
#'   string of color.
#' @param ... Other parameters passed to `genColorsForAnnots`.
#'
#' @return A HeatmapAnnotation object. Can be used for example in the
#'   `top_annotation` argument of `Heatmap`.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' data("DEgenesPrime_Naive")
#'
#' library(ComplexHeatmap)
#'
#' bestDE <- rownames(DEgenesPrime_Naive)[whichTop(DEgenesPrime_Naive$pvalue,
#'                                           decreasing = FALSE,
#'                                           top = 50)]
#' Heatmap(rowScale(bulkLogCounts[bestDE, ]) |> t(),
#'         right_annotation  = genRowAnnot(
#'           sampleAnnot[, c("culture_media", "line")])
#'         )
genRowAnnot <-
    function(annot,
            colorScales = NULL,
            border = TRUE,
            ...) {
        if (is.factor(annot) | is.data.frame(annot))
            annot <- droplevels(annot)
        if (is.vector(annot) | is.factor(annot)) {
            annot <- data.frame(Annotation = annot)
            if (is.list(colorScales))
                colnames(annot) <- names(colorScales)[1]
            if ((!is.null(colorScales)) & !is.list(colorScales))
                colorScales <- list(Annotation = colorScales)
        }
        colorScales <-
            genColorsForAnnots(
                annots = annot,
                colorScales = colorScales,
                returnContinuousFun = TRUE,
                ...
            )
        ComplexHeatmap::rowAnnotation(df = annot,
            col = colorScales,
            border = border)
    }

