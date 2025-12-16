
#' Principal Component Analysis
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#'      Can also be a SingleCellExperiment object.
#' @param transpose Logical. If `transpose`, samples are columns and features
#'   are rows.
#' @param scale Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param sce_assay Integer or character,
#'   if `data` is a `SingleCellExperiment` object, the assay name to use.
#' @return
#' A list with the following element:
#' - sdev: the standard deviations of the principal components (i.e.,
#'    the square roots of the eigenvalues of the covariance/correlation
#'    matrix, though the calculation is actually done with the singular
#'    values of the data matrix).
#' - rotation: the matrix of variable loadings
#'    (i.e., a matrix whose columns contain the eigenvectors).
#'    The function princomp returns this in the element loadings.
#' - center, scale: the centering and scaling used, or FALSE.
#' - x: the value of the rotated data (PC coordinates)
#' - n.obs: Number of samples
#' - propExplVar: proportion of explained variance by PCs on total variance.
#' - transform: a list containing the scaling (sdeviations) and centering
#'    factors (means) for each principal component.
#' - isFastPCA: Logical. If TRUE, the PCA was computed using the `fastPCA`
#'
#' If `data` is a `SingleCellExperiment.` returns a `SingleCellExperiment`
#' with the PCA computed in the `reducedDims` slot. The complete PCA object
#' is stored in the `metadata` slot at the index `PCA`.
#' @export
#'
#' @examples
#' data(iris)
#' pca<-PCA(iris[,seq_len(4)],transpose = FALSE,scale = TRUE,center = TRUE)
#'
#' iris_sce <- SingleCellExperiment(assays = list(counts = iris[,seq_len(4)]))
#' iris_sce<-PCA(iris_sce)
#' reducedDims(iris_sce)
PCA <- function(data,
                transpose = TRUE,
                scale = FALSE,
                center = TRUE,
                sce_assay = "logcounts") {
        sce_obj <-NULL
        if (inherits(data, "SingleCellExperiment")) {
            sce_obj <- data
            data <- assay(sce_obj, sce_assay)
        }
    if (transpose)
        data <- t(data)

    means <- 0
    sdeviations <- 1
    if (center) {
        means <- apply(data, 2, mean)
        data <- sweep(data, 2, means, "-")
    }
    if (scale) {
        sdeviations <- apply(data, 2, sd)
        data <- sweep(data, 2, sdeviations, "/")
    }
    resacp <- prcomp(
        x = data,
        retx = TRUE,
        center = FALSE,
        scale = FALSE
    )

    resacp$n.obs <- dim(data)[1]

    resacp$propExplVar <- resacp$sdev ^ 2 / sum(resacp$sdev ^ 2)
    resacp$scale <- scale
    resacp$center <- center
    resacp$transform <-
        list(sdeviations = sdeviations, means = means)
    resacp$isFastPCA <- FALSE

    if(is.null(sce_obj)){
        return(resacp)
    }else{
        reducedDim(sce_obj,"PCA") <- resacp$x
        metadata(sce_obj)$PCA <- resacp
        return(sce_obj)
    }

}


#' Faster Principal Component Analysis
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#'   Can also be a SingleCellExperiment object.
#' @param transpose Logical. If `transpose`, samples are columns and features
#'   are rows.
#' @param scale Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param nPC Integer. Number of Principal Component to be computed.
#' @param weight.by.var Logical. Weight PC by variables. If TRUE return a
#'   regular PCA.
#' @param sce_assay Integer or character,
#'   if `data` is a `SingleCellExperiment` object, the assay name to use.
#' @param ... Other parameters passed to `irlba`.
#'
#' @return
#' - sdev: the standard deviations of the principal components
#'    (i.e., the square roots of the eigenvalues of the covariance/correlation
#'    matrix, though the calculation is actually done with the singular
#'    values of the data matrix).
#' - rotation: the matrix of variable loadings
#'    (i.e., a matrix whose columns contain the eigenvectors).
#'    The function princomp returns this in the element loadings.
#' - center, scale: the centering and scaling used, or FALSE.
#' - x: the value of the rotated data (PC coordinates)
#' - n.obs: Number of samples
#' - propExplVar: proportion of explained variance by PCs on on variance of
#'    computed axes.
#' - transform: a list containing the scaling (sdeviations) and centering
#'    factors (means) for each principal component.
#' - isFastPCA: Logical. If TRUE, the PCA was computed using the `fastPCA`
#'
#' If `data` is a `SingleCellExperiment`, returns a `SingleCellExperiment`
#' with the PCA computed in the `reducedDims` slot. The complete PCA object
#' is stored in the `metadata` slot at the index `PCA`.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' pca<-fastPCA(bulkLogCounts,nPC=10)
#' pca2d(pca)
#'
#' sce <- SingleCellExperiment(assays = list(counts = bulkLogCounts))
#' sce<-fastPCA(sce,nPC=10)
#' reducedDims(sce)
fastPCA <-
    function(data,
            transpose = TRUE,
            scale = FALSE,
            center = TRUE,
            nPC = min(ncol(data) - 1, nrow(data) - 1, 30),
            weight.by.var = TRUE,
            sce_assay = "logcounts",
            ...) {

        sce_obj <-NULL
        if (inherits(data, "SingleCellExperiment")) {
                sce_obj <- data
                data <- assay(sce_obj, sce_assay)
        }else{
                data <- as.matrix(data)
        }
        if (transpose)
            data <- t(data)


        means <- 0
        sdeviations <- 1
        if (center | scale)
            data <- scale(data, scale = scale, center = center)
        if (center)
            means <- attr(data, "scaled:center")
        if (scale)
            sdeviations <- attr(data, "scaled:scale")

        resacp <- list()
        resacp$n.obs <- dim(data)[1]

        resacp$scale <- scale
        resacp$center <- center
        resacp$transform <-
            list(sdeviations = sdeviations, means = means)

        irlbaResults <- irlba::irlba(A = data, nv = nPC, ...)
        rotation <- irlbaResults$v
        resacp$sdev <- irlbaResults$d / sqrt(max(1, nrow(data) - 1))
        if (weight.by.var) {
            if (nPC > 1) {
                reducedSpace <- irlbaResults$u %*% diag(irlbaResults$d)
            } else {
                reducedSpace <- irlbaResults$u %*% irlbaResults$d
            }
        } else {
            reducedSpace <- irlbaResults$u
        }
        rownames(rotation) <- colnames(data)
        colnames(rotation) <- paste0("PC", seq_len(nPC))
        rownames(reducedSpace) <- rownames(data)
        colnames(reducedSpace) <- colnames(rotation)
        resacp$x <- reducedSpace
        resacp$rotation <- rotation
        resacp$propExplVar <- resacp$sdev ^ 2 / sum(resacp$sdev ^ 2)
        resacp$isFastPCA <- TRUE

        if(is.null(sce_obj)){
            return(resacp)
        }else{
            reducedDim(sce_obj,"PCA") <- resacp$x
            metadata(sce_obj)$PCA <- resacp
            return(sce_obj)
        }
    }

#' Add new samples to an existing PCA.
#'
#' @param pca The existing PCA as an object returned by `PCA` or `fastPCA`.
#'   Can also be a `SingleCellExperiment` object. Where one of these functions
#'   have been performed.
#' @param newSamplesMatrix The matrix of numeric of new samples. It must have
#'   the same features than the original PCA.
#' @param transpose Logical. If `transpose`, samples are columns and features
#'   are rows.
#' @param returnPCA Logical. Return the updated PCA. If false, just returned the
#'   coordinates of new samples on PCs.
#'
#' @return A PCA list if `returnPCA` or a matrix of coordinates of samples on
#'   PCs.
#'
#' @export
#'
#' @examples
#' data("iris")
#' iris1<-iris[1:75,]
#' iris2<-iris[76:150,]
#' pca<-PCA(iris1[,seq_len(4)],transpose = FALSE,scale = TRUE,center = TRUE)
#' pca2d(pca,colorBy = iris1$Species)
#' pcaUpdated<-pcaAddSamples(pca,iris2[,seq_len(4)],transpose = FALSE)
#' pca2d(pcaUpdated,colorBy = iris$Species)
#'
#' iris_sce <- SingleCellExperiment(
#'     assays = list(counts = t(iris1[,seq_len(4)])))
#' iris_sce<-PCA(iris_sce, scale = TRUE, center = TRUE)
#' pcaUpdated<-pcaAddSamples(iris_sce,iris2[,seq_len(4)],transpose = FALSE)
#' pca2d(pcaUpdated,colorBy = iris$Species)
pcaAddSamples <-
    function(pca,
            newSamplesMatrix,
            transpose = TRUE,
            returnPCA = TRUE) {
        sce_obj <-NULL
        if (inherits(pca, "SingleCellExperiment")) {
                sce_obj <- pca
                if (is.null(metadata(sce_obj)$PCA))
                        stop("PCA not found in metadata,",
                            " please run oob::PCA or oob::fastPCA first.")
                pca <- metadata(sce_obj)$PCA
        }
        if (transpose)
            newSamplesMatrix <- t(newSamplesMatrix)
        newSamplesMatrix <-
            newSamplesMatrix[, rownames(pca$rotation), drop = FALSE]

        if (pca$center)
            newSamplesMatrix <-
            sweep(newSamplesMatrix, 2, pca$transform$means, "-")
        if (pca$scale)
            newSamplesMatrix <-
            sweep(newSamplesMatrix, 2, pca$transform$sdeviations, "/")

        newSamplesCoord <-
            t(apply(newSamplesMatrix, 1, function(x) {
                colSums(x * pca$rotation)
            }))
        if (returnPCA) {
            pca$x <- rbind(pca$x, newSamplesCoord)
            return(pca)
        } else{
            return(newSamplesCoord)
        }
    }

#' Principal Component analysis of variance (Type II ANOVA)
#'
#' @param pca A PCA as an object returned by `PCA` or `fastPCA`.
#' @param colData A dataframe of feature (numeric or factor) with rows as samples.
#' @param nComponent Integer, number of PC used for the analysis.
#'
#' @return A dataframe with the PC, (partial) percentage of sum of squares and the
#'   corresponding p-value (Type II tests).
#' @export
PCaov <- function(pca, colData, nComponent = 10) {
	stopifnot(is.data.frame(colData))
	if (is.null(pca$x)) stop("`pca` must have a score matrix in `pca$x`.")
	if (!requireNamespace("car", quietly = TRUE)) {
		stop("Type II ANOVA requires the 'car' package. Install with: install.packages('car')")
	}

	pcs <- colnames(pca$x)
	if (is.null(pcs)) pcs <- paste0("PC", seq_len(ncol(pca$x)))
	nComponent <- min(nComponent, length(pcs))

	# Type II/III ANOVA is typically used with sum-to-zero contrasts
	old_contr <- options("contrasts")
	on.exit(options(old_contr), add = TRUE)
	options(contrasts = c("contr.sum", "contr.poly"))

	features <- colnames(colData)
	rhs <- paste(sprintf("`%s`", features), collapse = " + ")
	lm_formula <- stats::as.formula(paste0(".PC ~ ", rhs))

	anovPerPCdt <- do.call(
		"rbind",
		lapply(pcs[seq_len(nComponent)], function(PC) {
			md <- data.frame(colData, .PC = pca$x[, PC], check.names = FALSE)
			fit <- stats::lm(lm_formula, data = md, na.action = stats::na.omit)

			a2 <- car::Anova(fit, type = 2)

			# Extract in the same order as input features
			rn <- rownames(a2)
			idx <- match(features, rn)
			if (anyNA(idx)) {
				# fallback if matching fails for non-syntactic names
				tl <- attr(stats::terms(fit), "term.labels")
				idx <- match(tl, rn)
				feat_out <- tl
			} else {
				feat_out <- features
			}

			ss <- a2[idx, "Sum Sq"]
			pvals <- a2[idx, "Pr(>F)"]

			# For Type II SS, use partial SS (%): SS_term / (SS_term + SS_resid)
			rss <- sum(stats::residuals(fit)^2)
			sumSqPercent <- (ss / (ss + rss)) * 100

			data.frame(
				feature = feat_out,
				PC = PC,
				sumSqPercent = as.numeric(sumSqPercent),
				pval = as.numeric(pvals),
				row.names = NULL
			)
		})
	)

	anovPerPCdt$feature <- factor(anovPerPCdt$feature, levels = features)
	anovPerPCdt$padj <- p.adjust(anovPerPCdt$pval)

	# Keep your original PC formatting logic when possible
	pc_num <- suppressWarnings(as.integer(sub("^PC", "", anovPerPCdt$PC)))
	if (!all(is.na(pc_num))) {
		anovPerPCdt$PC <- paste0("PC", formatNumber2Character(as.character(pc_num)))
	}

	anovPerPCdt
}



#' UMAP projection
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#'   Can also be a SingleCellExperiment object.
#' @param nDimPCA Integer. If not NULL compute a PCA first and take n first
#'   principal components.
#' @param transpose Logical. If `transpose`, samples are columns and features
#'   are rows.
#' @param n_neighbors The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation. Larger values
#'   result in more global views of the manifold, while smaller values result in
#'   more local data being preserved. In general values should be in the range 2
#'   to 100.
#' @param n_components Integer. Dimensions of the embedded space.
#' @param min_dist The effective minimum distance between embedded points.
#'   Smaller values will result in a more clustered/clumped embedding where
#'   nearby points on the manifold are drawn closer together, while larger
#'   values will result on a more even dispersal of points. The value should be
#'   set relative to the spread value, which determines the scale at which
#'   embedded points will be spread out.
#' @param init Type of initialization for the coordinates (see `uwot` for more
#'   details).
#' @param metric Type of distance metric to use to find nearest neighbors (see
#'   `uwot` for more details).
#' @param ret_model If TRUE, then return extra data that can be used to add new
#'   data to an existing embedding via umap_transform. The embedded coordinates
#'   are returned as the list item embedding. If FALSE, just return the
#'   coordinates. This parameter can be used in conjunction with ret_nn and
#'   ret_extra. Note that some settings are incompatible with the production of
#'   a UMAP model: external neighbor data (passed via a list to nn_method), and
#'   factor columns that were included via the metric parameter. In the latter
#'   case, the model produced is based only on the numeric data. A
#'   transformation using new data is possible, but the factor columns in the
#'   new data are ignored
#' @param ret_nn If TRUE, then in addition to the embedding, also return nearest
#'   neighbor data that can be used as input to nn_method to avoid the overhead
#'   of repeatedly calculating the nearest neighbors when manipulating unrelated
#'   parameters (e.g. min_dist, n_epochs, init). See the "Value" section for the
#'   names of the list items. If FALSE, just return the coordinates. Note that
#'   the nearest neighbors could be sensitive to data scaling, so be wary of
#'   reusing nearest neighbor data if modifying the scale parameter. This
#'   parameter can be used in conjunction with ret_model and ret_extra.
#' @param sce_assay Integer or character,
#'   if `data` is a `SingleCellExperiment` object, the assay name to use.
#' @param ... Other parameters passed to uwot.
#'
#' @return A matrix of coordinates with samples as rows, or:
#' - if ret_model = TRUE (or ret_extra contains "model"),
#'    returns a list containing extra information that can be used to add new
#'    data to an existing embedding via umap_transform.
#'    In this case, the coordinates are available in the list item embedding.
#'    NOTE: The contents of the model list should not be considered stable or
#'    part of the public API, and are purposely left undocumented.
#' - if ret_nn = TRUE (or ret_extra contains "nn"),
#' returns the nearest neighbor data as a list called nn.
#' This contains one list for each metric calculated,
#' itself containing a matrix idx with the integer ids of the neighbors;
#' and a matrix dist with the distances. The nn list (or a sub-list)
#' can be used as input to the nn_method parameter.
#' - if data is a SingleCellExperiment object, returns a SingleCellExperiment
#' object with the UMAP coordinates in the `reducedDims` slot, and the eventual
#' complete UMAP model and nearest neighbor in the metadata slot at the index
#' `UMAP`.
#' @export
#'
#' @examples
#' data(iris)
#' irisUMAP<-UMAP(iris[,seq_len(4)],transpose = FALSE)
#' proj2d(irisUMAP,colorBy = iris$Species)
#' irisUMAP<-UMAP(rowScale(iris[,seq_len(4)],center = TRUE,scaled = TRUE),
#'     transpose = FALSE,n_neighbors = nrow(iris),
#'     ret_nn = TRUE, ret_model = TRUE)
#' proj2d(irisUMAP$embedding,colorBy = iris$Species,
#'     nnMatrix = irisUMAP$nn$euclidean$idx[,c(1, 2, 3)],fixedCoord = TRUE)
#'
#' irisSCE <- SingleCellExperiment::SingleCellExperiment(
#'         assays = list(counts = t(iris[, -5])),
#'         colData = iris[5]
#' )
#' irisSCE<-UMAP(irisSCE)
UMAP <-
    function(data,
            nDimPCA = NULL,
            transpose = TRUE,
            n_neighbors = 20,
            n_components = 2,
            min_dist = 0.01,
            init = "laplacian",
            metric = "euclidean",
            ret_model = FALSE,
            ret_nn = FALSE,
            sce_assay = "logcounts",
            ...) {
    sce_obj <-NULL
    if (inherits(data, "SingleCellExperiment")) {
            sce_obj <- data
            data <- assay(sce_obj, sce_assay)
    }
    if (transpose)
        data <- t(data)
    if (nrow(data) < n_neighbors) {
        n_neighbors <- nrow(data)
        warning(
            "n_neighbors must not exceed number of samples, ",
            "adjusting n_neighbors to number of samples (",
            n_neighbors,
            ")"
        )
    }
    if (is.null(n_neighbors))
        n_neighbors <- nrow(data)
    if (!is.null(nDimPCA)) {
        data <- fastPCA(data,
                        transpose = FALSE,
                        scale = FALSE,
                        nPC = nDimPCA)$x
    }
    res <-
        uwot::umap(
            as.matrix(data),
            n_neighbors = n_neighbors,
            n_components = n_components,
            min_dist = min_dist,
            init = init,
            metric = metric,
            ret_model = ret_model,
            ret_nn = ret_nn,
            ...
        )
    if (!ret_model & !ret_nn)
        rownames(res) <- rownames(data)
    if(is.null(sce_obj)){
            return(res)
    }else{
        if(ret_model | ret_nn){
            reducedDim(sce_obj,"UMAP") <- res$embedding
            metadata(sce_obj)$UMAP <- res
        }else{
            reducedDim(sce_obj,"UMAP") <- res
        }
        return(sce_obj)
}
}


#' TriMap dimension reduction.
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#'    Can also be a SingleCellExperiment object.
#' @param n_dims Integer. Dimensions of the embedded space.
#' @param transpose Logical. If `transpose`, samples are columns and features
#'   are rows.
#' @param n_inliers Number of nearest neighbors for forming the nearest neighbor
#'   triplets.
#' @param n_outliers Number of outliers for forming the nearest neighbor
#'   triplets.
#' @param apply_pca Reduce the number of dimensions of the data to 100 if
#'   necessary before applying the nearest-neighbor search.
#' @param n_iters Number of iterations.
#' @param knn_tuple Use the precomputed nearest-neighbors information in form of
#'   a tuple (knn_nbrs, knn_distances).
#' @param sce_assay Integer or character,
#'   if `data` is a `SingleCellExperiment` object, the assay name to use.
#'
#' @return A matrix of coordinates with samples as rows.
#'
#' If `data` is a `SingleCellExperiment` returns a `SingleCellExperiment`
#' with the TriMap computed in the `reducedDims` slot at the index `TriMap`.
#' @export
#'
#' @examples
#' data(iris)
#' irisTrimap<-TRIMAP(iris[,seq_len(4)],transpose=FALSE)
#' proj2d(irisTrimap,colorBy = iris$Species)
#'
#' irisSCE <- SingleCellExperiment::SingleCellExperiment(
#'         assays = list(counts = t(iris[, -5])),
#'         colData = iris[5]
#' )
#' irisSCE<-TRIMAP(irisSCE)
TRIMAP <-
    function(data,
            n_dims = 2,
            transpose = TRUE,
            n_inliers = 10,
            n_outliers = 5,
            apply_pca = TRUE,
            n_iters = 400,
            knn_tuple = NULL,
            sce_assay = "logcounts") {
    sce_obj <-NULL
    if (inherits(data, "SingleCellExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }
    if (transpose)
        data <- t(data)

    proc <- basiliskStart(pythonEnv)
    on.exit(basiliskStop(proc))

    trimapOut<-basiliskRun(proc,
        fun=function(n_dims, n_inliers,n_outliers,
            apply_pca, n_iters, knn_tuple, data) {

        trimap_module <-
            reticulate::import(module = "trimap", delay_load = TRUE)
        trimap <- trimap_module$TRIMAP(
            n_dims = as.integer(n_dims),
            n_inliers = as.integer(n_inliers),
            n_outliers = as.integer(n_outliers),
            apply_pca = apply_pca,
            n_iters = as.integer(n_iters),
            knn_tuple = knn_tuple
        )
        res <- trimap$fit_transform(as.matrix(data))
        rownames(res) <- rownames(data)
        res
    }, n_dims, n_inliers, n_outliers, apply_pca, n_iters, knn_tuple, data)
    if(is.null(sce_obj)){
        return(trimapOut)
    }else{
        reducedDim(sce_obj,"TriMap") <- trimapOut
        return(sce_obj)
    }
}


#' Non-metric Multidimensional Scaling (dimension reduction)
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#'   Can also be a SingleCellExperiment object.
#' @param transpose Logical. If `transpose`, samples are columns and features
#'   are rows.
#' @param scale  Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param metric A function that return a object of class "dist".
#' @param ndim Integer. Dimensions of the embedded space.
#' @param maxit The maximum number of iterations.
#' @param sce_assay Integer or character,
#'   if `data` is a `SingleCellExperiment` object, the assay name to use.
#' @return A matrix of coordinates with samples as rows, or a
#'   `SingleCellExperiment` object if `data` is a `SingleCellExperiment.` with
#'   the NMDS stored in a `reducedDim` slot.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' NMDSproj<-NMDS(bulkLogCounts)
#' proj2d(NMDSproj)
#' sce <- SingleCellExperiment(assays = list(logcounts = bulkLogCounts))
#' sce<-NMDS(sce)
#' proj2d(sce)
NMDS <-
    function(data,
        transpose = TRUE,
        scale = FALSE,
        center = FALSE,
        metric = dist,
        ndim = 2,
        maxit = 100,
        sce_assay = "logcounts") {
    sce_obj <-NULL
    if (inherits(data, "SingleCellExperiment")) {
        sce_obj <- data
        data <- assay(sce_obj, sce_assay)
    }
    merged <- FALSE
    if (transpose)
        data <- t(data)
    data <- metric(data)  # euclidean distances between the rows
    if (min(data, na.rm = TRUE) == 0) {
        merged <- TRUE
        md <- merge0dist(data)
        data <- md$distMat
        mergedSample <- md$merged
    }
    fit <-
        MASS::isoMDS(data, k = ndim, maxit = maxit) # k is the number of dim
    fit$coord <- fit$points
    fit$points <- NULL
    if (merged) {
        for (sple in names(mergedSample)) {
            values <-
                matrix(
                    rep(fit$coord[sple, ], length(mergedSample[[sple]])),
                    nrow = length(mergedSample[[sple]]),
                    byrow = TRUE
                )
            rownames(values) <- mergedSample[[sple]]
            fit$coord <- rbind(fit$coord, values)
        }
    }
    if(is.null(sce_obj)){
        return(fit$coord)
    }else{
        reducedDim(sce_obj,"NMDS") <- fit$coord
        return(sce_obj)
    }
}

