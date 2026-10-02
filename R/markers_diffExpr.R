
#' Compute a matrix of coef describing best gene marker per group of samples
#' from a GLM net regression
#'
#' @param data A matrix of numeric with rows as features (in the
#'   RNA-Seq context, log data).
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param group A feature of factor/character, same length as number of sample.
#'   Describe group of each sample.
#' @param transpose If TRUE, the input data is transposed before processing.
#'   Default is TRUE (feature as rows, samples as columns).
#' @param sce_assay Integer or character, if `data` is a
#'   `SummarizedExperiment` related object, the assay name to use.
#'
#' @return A matrix containing each gene coefficient for each group.
#'   If `data` is a `SummarizedExperiment` related object
#'   the function will add the marker metrics to `rowData`.
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' res <- getMarkerGLMnet(bulkLogCounts,sampleAnnot$culture_media)
#' sce <- SingleCellExperiment(assays = list(counts = bulkLogCounts))
#' res <- getMarkerGLMnet(sce,sampleAnnot$culture_media)
getMarkerGLMnet <- function(data, group, transpose = TRUE, sce_assay = 1) {
    sce_obj <-NULL
    if (inherits(data, "SummarizedExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }
    data <- as.matrix(data)
    if (transpose) data <- t(data)
    fit <- glmnet(
        data,
        group,
        family = "multinomial",
        alpha = .5,
        lambda = cv.glmnet(data, group, family = "multinomial")$lambda.1se
    )

    cf <- vapply(coef(fit), function(x)
        x[seq(2,length(x))], numeric(ncol(data)))
    rownames(cf) <- colnames(data)
        colnames(cf) <- paste0("coef_",colnames(cf))

    if(is.null(sce_obj)){
        return(cf)
    }else{
        rowData(sce_obj) <- cbind(rowData(sce_obj),cf)
        return(sce_obj)
    }
}

#' Compute over dispersion values for each gene.
#'
#' @param data Normalized count table with genes as rows.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param minCount Minimum average expression to not be filtered out.
#' @param plot Logical. Show the overdispersion plot.
#' @param returnPlot Logical, if `plot` return it as a ggplot object instead of
#'   printing it.
#' @param sce_assay Integer or character, if `data` is a
#'   `SummarizedExperiment` related object, the assay name to use.
#'
#' @return A ggplot graph if `returnPlot`, otherwise a dataframe with the
#'   following columns:
#' - mu: average expression
#' - var: variance
#' - cv2: squared coefficient of variation. Used as a dispersion value.
#' - residuals: y-distance from teh regression. Can be used as an
#'    overdispersion value.
#' - residuals2: squared residuals
#' - fitted: theoretical dispersion for the gene average (y value of the curve).
#'
#'   If `data` is a `SummarizedExperiment` related object
#'   the function will add the gene metrics to `rowData`.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' normCount<-2^(bulkLogCounts-1)
#' dispData<-getMostVariableGenes(normCount,minCount=1)
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts = normCount))
#' sce <- getMostVariableGenes(sce,minCount=1)
#' rowData(sce) |> chead()
getMostVariableGenes <-
    function(data,
            minCount = 0.01,
            plot = TRUE,
            returnPlot = FALSE,
            sce_assay = 1) {

    sce_obj <-NULL
    if (inherits(data, "SummarizedExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }
    data <- data[rowMeans(data) > minCount, ]
    dispTable <-
        data.frame(
            mu = rowMeans(data),
            var = apply(data, 1, var),
            row.names = rownames(data)
        )
    dispTable$cv2 <- dispTable$var / dispTable$mu ^ 2
    sumNullvariance <- sum(dispTable$cv2 <= 0)
    if (sumNullvariance > 0) {
        warning(sumNullvariance, " have null variance and will be removed")
        dispTable <- dispTable[dispTable$cv2 > 0, ]
    }
    fit <-
        loess(as.formula("cv2 ~ mu"),
            data = log10(dispTable[, c("mu", "cv2")]))
    dispTable$residuals <- fit$residuals
    dispTable$residuals2 <- dispTable$residuals ^ 2
    dispTable$fitted <- 10 ^ fit$fitted
    if (plot) {
        g <-
            ggplot(dispTable,
                aes(
                    x = .data$mu,
                    y = .data$cv2,
                    label = rownames(dispTable),
                    fill = .data$residuals
                )) +
            geom_point(stroke = 1 / 8,
                        colour = "black",
                        shape = 21) +
            geom_line(aes(y = fitted), color = "red", size = 1.5) +
            scale_x_log10() + scale_y_log10()
        if (returnPlot) {
            return(g)
        } else{
            print(g)
        }
    }
    if(is.null(sce_obj)){
        return(dispTable)
    }else{
        sce_obj <- sce_obj[rownames(dispTable),]
        rowData(sce_obj) <- cbind(rowData(sce_obj),dispTable)
        return(sce_obj)
    }
}
