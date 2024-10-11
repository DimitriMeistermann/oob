#' Give QC metrics per sample
#'
#' @param x A matrix of numeric with samples as columns.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param uncenter uncenter the matrix first so the minimum is 0.
#' @param sce_assay Integer or character, if `data` is a `SummarizedExperiment`
#'    related object, the assay name to use.
#'
#' @return A data frame with each row as a sample and the following columns:
#' * `mean`: average expression in the sample
#' * `sd`: standard deviation
#' * `TotalGenEx`: number of expressed gene (count>0) in the sample
#' * `TotalCount`: sum of counts for the sample
#'
#' If `x` is a `SummarizedExperiment` related object, the function will
#' add the columns to the `colData` slot.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' countMat <-
#'     sapply(vector("numeric", length = length(geneLengthGRCh38)),
#'         function(x) {
#'             MASS::rnegbin(10, theta = abs(rnorm(1, mean = 10, sd = 20)),
#'             mu = abs(rnorm(1, mean = 10, sd = 20)))
#'     }) |> t()
#' rownames(countMat) <- names(geneLengthGRCh38)
#' colnames(countMat) <- letters[seq_len(ncol(countMat))]
#' computeQCmetricSamples(countMat)
#' sce <- SingleCellExperiment(assays = list(counts = countMat))
#' sce <- computeQCmetricSamples(sce)
#' colData(sce) |> head()
computeQCmetricSamples <- function(x, uncenter = FALSE, sce_assay = 1) {
    sce_obj <-NULL
    if (inherits(x, "SummarizedExperiment")) {
        sce_obj <- x
        x <- assay(sce_obj, sce_assay)
    }
    if (uncenter) {
        x <- x - min(x)
        zero <- 0
    } else{
        zero <- min(x)
    }
    mean <- apply(x, 2, mean)
    sd <- apply(x, 2, sd)
    count <- colSums(x)
    CV <- apply(x, 2, cv)
    noGenEx <- rep(0, ncol(x))
    for (i in seq_len(ncol(x)) )
        noGenEx[i] <- length(which(x[, i] > zero))

    resDt<-data.frame(
        mean = mean,
        sd = sd,
        CV = CV,
        TotalGenEx = noGenEx,
        TotalCount = count
    )
    if(is.null(sce_obj)){
        return(resDt)
    }else{
        colData(sce_obj) <- cbind(colData(sce_obj), resDt)
        return(sce_obj)
    }
}



merge0dist <- function(disMat) {
    mat <- as.matrix(disMat)
    merged <- list()
    found <- TRUE
    while (found == TRUE) {
        found <- FALSE
        for (i in 2:nrow(mat)) {
            for (j in seq_len(i - 1)) {
                if (mat[i, j] == 0) {
                    newNames <- rownames(mat)
                    newNames <- newNames[-i]
                    newMat <- mat[-i, -i]
                    colnames(newMat) <- rownames(newMat) <- newNames
                    merged[[rownames(mat)[j]]] <-
                        c(merged[[rownames(mat)[j]]], rownames(mat)[i])
                    mat <- newMat
                    found <- TRUE
                    break
                }
            }
            if (found)
                break
        }
    }
    return(list(distMat = as.dist(mat), merged = merged))
}

#' Count Per Million normalization
#'
#' @param data A raw count table with rows as genes
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param sce_assay Integer or character,
#'   if `data` is a `SummarizedExperiment` related object,
#'   the assay name to use.
#' @return A normalized count table.
#'   If `data` is a `SummarizedExperiment` related object, 
#'   the function will add the normalized 
#'   count to the `SummarizedExperiment` object `cpm`.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' countMat <-
#'     sapply(vector("numeric", length = length(geneLengthGRCh38)), function(x){
#'         MASS::rnegbin(10, theta = abs(rnorm(1, mean = 10, sd = 20)),
#'                      mu = abs(rnorm(1, mean = 10, sd = 20)))
#'     }) |> t()
#'     rownames(countMat)<-names(geneLengthGRCh38)
#' head(CPM(countMat))
#' sce <- SingleCellExperiment(assays = list(counts = countMat))
#' head(assay(CPM(sce),"cpm"))
CPM <- function(data, sce_assay = 1) {
    sce_obj <-NULL
    if (inherits(data, "SummarizedExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }
    #Normalisation CPM
    data.CPM <- sweep(data, 2, colSums(data), `/`)
    data.CPM <- data.CPM * 1000000

    if(is.null(sce_obj)){
        return(data.CPM)
    }else{
        assay(sce_obj,"cpm") <- data.CPM
        return(sce_obj)
    }
}


#' Transcript per million (TPM) normalization for full-length transcript, short
#' read sequencing.
#'
#' @param data A raw count table with rows as genes.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param gene.length A vector of numeric corresponding to gene length in base
#'   pair. Must be named by genes.
#' @param sce_assay Integer or character,
#'   if `data` is a `SummarizedExperiment` related object,
#'   the assay name to use.
#' @return A normalized count table where 1 count is equal in theory to one
#'   transcript per million.
#'   If `data` is a `SummarizedExperiment` related object, the function will
#'   add the normalized count to the `SummarizedExperiment` object as an
#'   assay named `tpm`.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' countMat <-
#'     sapply(vector("numeric", length = length(geneLengthGRCh38)), function(x){
#'         MASS::rnegbin(10, theta = abs(rnorm(1, mean = 10, sd = 20)),
#'                      mu = abs(rnorm(1, mean = 10, sd = 20)))
#'     }) |> t()
#'     rownames(countMat)<-names(geneLengthGRCh38)
#' TPMfullLength(countMat,geneLengthGRCh38) |> chead()
#'
#' sce <- SingleCellExperiment(assays = list(counts = countMat))
#' assay(TPMfullLength(sce,geneLengthGRCh38),"tpm") |> head()
TPMfullLength <- function(data, gene.length, sce_assay = 1) {
    sce_obj <-NULL
    if (inherits(data, "SummarizedExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }
    gene.length.kb <- gene.length[rownames(data)] / 1000
    data <- sweep(data, 1, gene.length.kb, `/`)
    res<-CPM(data)
    if(is.null(sce_obj)){
        return(res)
    }else{
        assay(sce_obj,"tpm") <- res
        return(sce_obj)
    }
}

#' Reads Per Kilobase per Million (RPKM) normalization for full-length
#' transcript, short read sequencing.
#'
#' @param data A raw count table with rows as genes.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment`  object.
#' @param gene.length A vector of numeric corresponding to gene length in base
#'   pair. Must be named by genes.
#' @param sce_assay Integer or character, if `data` is a
#'   `SummarizedExperiment` related object, the assay name to use.
#' @return A normalized count table of RPKM.
#'   If `data` is a `SummarizedExperiment` object, the function will add the
#'   normalized count to the `SummarizedExperiment` object as an assay named
#'   `rpkm`.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' countMat <-
#'     sapply(vector("numeric", length = length(geneLengthGRCh38)), function(x){
#'         MASS::rnegbin(10, theta = abs(rnorm(1, mean = 10, sd = 20)),
#'                      mu = abs(rnorm(1, mean = 10, sd = 20)))
#'     }) |> t()
#' rownames(countMat)<-names(geneLengthGRCh38)
#' chead(RPKM(countMat,geneLengthGRCh38))
#'
#' sce <- SingleCellExperiment(assays = list(counts = countMat))
#' assay(RPKM(sce,geneLengthGRCh38),"rpkm") |> head()
RPKM <- function(data, gene.length, sce_assay = 1) {
    sce_obj <-NULL
    if (inherits(data, "SummarizedExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }

    gene.length.kb <- gene.length[rn(data)] / 1000
    data <- CPM(data)
    res <- sweep(data, 1, gene.length.kb, `/`)
    if(is.null(sce_obj)){
        return(res)
    }else{
        assay(sce_obj,"rpkm") <- res
        return(sce_obj)
    }
}


#' Quick DESeq2 normalization
#'
#' @param countMatrix A raw count table with rows as genes.
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment`  object.
#' @param sce_assay Integer or character,
#'   if `data` is a `SummarizedExperiment` related object,
#'   the assay name to use.
#' @return A normalized count table.
#' If `data` is a `SummarizedExperiment` related object, the function will
#'   add the normalized count to
#'   the `SummarizedExperiment` object as an assay named `deseq2norm`.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' countMat <-
#'     sapply(vector("numeric", length = length(geneLengthGRCh38)), function(x){
#'         MASS::rnegbin(10, theta = abs(rnorm(1, mean = 10, sd = 20)),
#'         mu = abs(rnorm(1, mean = 10, sd = 20)))
#'     }) |> t()
#' rownames(countMat)<-names(geneLengthGRCh38)
#' normDeseq(countMat)
#' sce <- SingleCellExperiment(assays = list(counts = countMat))
#' assay(normDeseq(sce),"deseq2norm") |> head()
normDeseq <-function(countMatrix, sce_assay = 1) {
    sce_obj <-NULL
    if (inherits(countMatrix, "SummarizedExperiment")) {
        sce_obj <- countMatrix
        countMatrix <- assay(sce_obj, sce_assay)
    }
    #matrix where genes are rows and samples are columns
    # PS = pseudo reference sample
    PS <-
        apply(countMatrix, 1, gmean, keepZero = TRUE)
    # get a vector which consist of the geometrical mean
    # of each genes across all samples
    keptRow <-
        PS > 0 #get rid of genes containing one zero ore more
    PS <- PS[keptRow]
    ratioMat <- sweep(countMatrix[keptRow,], 1, PS, "/")
    # get the ratio matrix (expression/expression from PS)
    normFactors <-
        apply(ratioMat, 2, median)
    # get the median of the ratios for each sample
    # to get the normalization factors
    res<-sweep(countMatrix, 2, normFactors, "/")
    #divide each sample by the corresponding normalization factor
    if(is.null(sce_obj)){
        return(res)
    }else{
        assay(sce_obj,"deseq2norm") <- res
        return(sce_obj)
    }

}

# look at clusters <- quickCluster(sce)
#' Quick single cell normalization (scran method).
#'
#' @param rawCounts RNA-Seq raw counts with rows as genes.
#' @param returnLog Logical. Return log counts.
#' @param sizeFactors NULL or a vector of numeric containing precomputed size
#'   factor, same size as number of cells in `rawCounts`.
#' @param ... Other parameter passed to `computeSumFactors`.
#'
#' @return A normalized count table.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' countMat <-
#'     sapply(vector("numeric", length = length(geneLengthGRCh38)), function(x){
#'         MASS::rnegbin(10, theta = abs(rnorm(1, mean = 10, sd = 20)),
#'             mu = abs(rnorm(1, mean = 10, sd = 20)))
#'     }) |> t()
#' rownames(countMat)<-names(geneLengthGRCh38)
#' quickSCnorm(countMat, returnLog=FALSE)
#' quickSCnorm(countMat, returnLog=TRUE)
quickSCnorm <-
    function(rawCounts,
            returnLog = TRUE,
            sizeFactors = NULL,
            ...) {
        sce <-
            SingleCellExperiment(assays = list(counts = rawCounts))
        if (!is.null(sizeFactors)) {
            BiocGenerics::sizeFactors(sce) <- sizeFactors
        } else{
            sce <- scran::computeSumFactors(sce, ...)
        }
        if (returnLog) {
            scater::normalizeCounts(
                sce,
                transform = "log",
                pseudo_count = 1,
                size_factors = sizeFactors(sce)
            )
        } else{
            scater::normalizeCounts(
                sce,
                transform = "none",
                pseudo_count = 0,
                size_factors = sizeFactors(sce)
            )
        }
    }


#' Correlation from one gene to all others.
#'
#' @param gene A single charachter. The gene (or feature) name that will be used
#'   for correlating to all others.
#' @param expression A matrix of numeric with rows as features (in the RNA-Seq
#'   context, log counts).
#' @param corFun A function to compute a correlation between two feature.
#' @param ... Parameters passed to `corFun`.
#'
#' @return A vector of numeric. Correlations values named by their corresponding
#'   gene.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' corGeneToOthers("NANOG",bulkLogCounts) |> head()
corGeneToOthers <- function(gene, expression, corFun = cor, ...) {
    expression <- as.matrix(expression)
    t(corFun(expression[gene, ], t(expression), ...))[, 1]
}

#' Execute a fastMNN and rescale the counts.
#'
#' @description Use fastMNN fron the package batchelor, but return a matrix.
#' Rescale the counts so the range of each gene remains the same after the
#' transformation.
#' @param logCounts A matrix of numeric (in the RNA-Seq context, log counts).
#'   Can also be a `SummarizedExperiment` or `SingleCellExperiment` object.
#' @param batch A vector or factor specifying the batch of origin for all cells
#' @param k An integer scalar specifying the number of nearest neighbors to
#'   consider when identifying MNNs.
#' @param returnRescale Logical. Use reScale on the output so the dynamic range
#'   after batch correction of genes is the same than before.
#' @param sce_assay Integer or character, if `data`
#'   is a `SummarizedExperiment` related object, the assay name to use.
#' @param ... Other parameters passed to fastMNN.
#'
#' @return A matrix of corrected count table.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' countMat<-sapply(vector("numeric",length = 100),function(x){
#'      c(MASS::rnegbin(10,mu = 50,theta = 5),
#'        MASS::rnegbin(10,mu = 10,theta = 5))
#' }) |> t(); countMat<-log2(countMat+1)
#' heatmap.DM(countMat)
#' #warning because of the small matrix
#' correctedMat<-oobFastMNN(countMat,batch = c(rep(1,10),rep(2,10)), k=5)
#' heatmap.DM(correctedMat)
oobFastMNN <- function(logCounts, batch, k,
    returnRescale = TRUE, sce_assay = 1, ...) {

    if (inherits(data, "SummarizedExperiment")) {
        data <- assay(data, sce_assay)
    }
    scObj <- batchelor::fastMNN(logCounts, batch = batch, k = k, ...)
    res <-
        SummarizedExperiment::assay(scObj, "reconstructed") |> as.matrix()
    if (returnRescale)
        res <- reScale(res, logCounts)
    res
}
