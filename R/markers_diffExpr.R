#' Compute a dataframe with marker metrics describing best gene marker per group
#' of samples.
#'
#' @param BPPARAM A BPPARAM object as return by [BiocParallel::bpparam()]. Used
#'   for multi-threading.
#' @param returnAsList Return a list where each element is dataframe containing
#'   the marker metrics of a group.
#' @param data A dataframe with genes as rows and samples as columns.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param groups A vector of group names, same size as the number of columns in
#'   `data`.
#' @param transpose If TRUE, the input data is transposed before processing.
#'   Default is TRUE (feature as rows, samples as columns).
#' @param sce_assay Integer or character, if `data` is a
#'   `SummarizedExperiment` related object, the assay name to use.
#'
#' @return A dataframe containing four column per group: Log2(Fold-Change),
#'   AUROC, marker score (see details), p-value and BH adjusted p-value.
#'   If `data` is a `SummarizedExperiment` related object and `returnAsList` is
#'   `FALSE`, the function will add the marker metrics to `rowData`.
#' @details LogFC and pvalues are computed from a linear modelling of the data.
#'
#' Score is consisting of the geometrical mean of absolute LogFC, absolute(auroc
#' - 0.5), and -log10(pval), then signed by the logFC: score = sign(logFC) ×
#' gmean( abs(logFC), abs(aurocRes-0.5), -log10(pval) )
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' # SnowParam(1) : monocore execution
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media,
#'     BPPARAM=BiocParallel::SnowParam(1))
#' sce <- SingleCellExperiment(assays = list(counts = bulkLogCounts))
#' sce <- getMarkers(sce, sampleAnnot$culture_media, sce_assay = "counts")
#' rowData(sce) |> chead()
getMarkers <- function (data,
                        groups,
                        transpose = TRUE,
                        BPPARAM = NULL,
                        returnAsList = FALSE,
                        sce_assay = 1) {
    if (is.null(BPPARAM))
        BPPARAM <- BiocParallel::bpparam()

    sce_obj <-NULL
    if (inherits(data, "SummarizedExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }

    if(!transpose) data <- t(data)
    groups <- as.factor(make.names(groups))
    grpLvl <- levels(groups)
    resDF <- BiocParallel::bplapply(grpLvl,BPPARAM = BPPARAM,
        FUN = function(group) {
            logicGroup <- groups == group
            designMat <- cbind(1, logicGroup)
            n1 <- as.integer(sum(!logicGroup))
            n2 <- as.integer(sum(logicGroup))
            resDFgroup <- apply(data, 1,
                                getGeneMarkerStat, designMat, n1, n2) |>
                t() |> data.frame()
            colnames(resDFgroup) <-
                c("lfc", "auroc", "score", "pval")
            resDFgroup$padj <-
                p.adjust(resDFgroup$pval, method = "BH")
            resDFgroup
        })
    names(resDF) <- grpLvl

    if (returnAsList) return(resDF)
        resDF<-do.call("cbind", resDF)
    if(is.null(sce_obj)){
        return(resDF)
    }else{
        rowData(sce_obj) <- cbind(rowData(sce_obj),resDF)
        return(sce_obj)
    }

}


getGeneMarkerStat<-function(gene, designMat, n1, n2){
    aurocRes <- aurocCPP(
        score = gene,
        boolVect = designMat[, 2],
        n1 = n1,
        n2 = n2
    )
    lmOut <- .lm.fit(designMat, gene)
    pval <- pvalLmFit(
        lmOut$residuals,
        lmOut$coefficients,
        p = lmOut$rank,
        qr = lmOut$qr
    )[2]
    coef <- lmOut$coefficients[2]
    score <-
        sign(coef) * prod(c(abs(coef), abs(aurocRes - 0.5),
                            min(-log10(pval), 324))) ^ (1 / 3)
    #min(-log10(pval),324): avoid Inf
    return(c(coef, aurocRes, score, pval))
}


#' Approximation of area under ROC curve
#'
#' @description By Miron Kursa https://mbq.me
#'
#' @param score a vector of numeric representing the measure of a feature in a
#'   set of samples.
#' @param boolVect a vector of logical, same size as score. Is the sample in the
#'   target group?
#'
#' @return A numeric value. Close to 1 = perfect marker, around 0.5 = as good as
#'   random values, close to 0 = perfect anti-marker.
#' @export
#'
#' @examples
#' data(iris)
#' auroc(iris$Sepal.Length,iris$Species=="virginica")
auroc <- function(score, boolVect) {
    n1 <- sum(!boolVect)
    n2 <- sum(boolVect)
    U    <- sum(rank(score)[!boolVect]) - n1 * (n1 + 1) / 2
    return(1 - U / n1 / n2)
}

#' Quick approximation of area under ROC curve
#'
#' @description By Miron Kursa https://mbq.me
#'
#' @param score a vector of numeric representing the measure of a feature in a
#'   set of samples.
#' @param boolVect a vector of logical, same size as score. Is the sample in the
#'   target group?
#'
#' @return A numeric value. Close to 1 = perfect marker, around 0.5 = as good as
#'   random values, close to 0 = perfect anti-marker.
#' @export
#'
#' @examples
#' data(iris)
#' qauroc(iris$Sepal.Length,iris$Species=="virginica")
qauroc <- function(score, boolVect) {
    n1 <- sum(!boolVect)
    n2 <- sum(boolVect)
    aurocCPP(as.numeric(score),
            as.logical(boolVect),
            as.integer(n1),
            as.integer(n2))
}

#' Extract a specific feature/metric (pval, logFC...) from a marker result
#' dataframe
#'
#' @param markerData A dataframe returned by [oob::getMarkers].
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object where
#'   were [oob::getMarkers] has been performed.
#' @param feature The name of the feature that has to be extracted.
#'
#' @return
#' A matrix containing only the wanted feature where each column is a group.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media,
#'     BPPARAM=BiocParallel::SnowParam(1))
#' extractFeatureMarkerData(markerData) |> chead()
#' sc <- SingleCellExperiment(assays = list(counts = bulkLogCounts),
#'    colData = sampleAnnot,
#'    rowData = markerData)
#' extractFeatureMarkerData(sc) |> chead()
extractFeatureMarkerData <-
    function(markerData, feature = "score") {
        if (inherits(markerData, "SummarizedExperiment")) {
            markerData <- rowData(markerData)
        }
        featureStrLen <- nchar(feature) + 1
        columns2Pick <-
            grep(paste0("^.*\\.", feature, "$"),
                colnames(markerData),
                value = TRUE)
        markerData <- markerData[, columns2Pick, drop = FALSE] |> as.matrix()
        grpName <- colnames(markerData)
        colnames(markerData) <-
            substr(grpName, 1, nchar(grpName) - featureStrLen)
        return(markerData)
    }

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
