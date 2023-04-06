#' Give QC metrics per sample
#'
#' @param x A matrix of numeric with samples as columns.
#' @param uncenter uncenter the matrix first so the minimum is 0.
#'
#' @return A data frame with each row as a sample and the following columns:
#' * `mean`: average expression in the sample
#' * `sd`: standard deviation
#' * `TotalGenEx`: number of expressed gene (count>0) in the sample
#' * `TotalCount`: sum of counts for the sample
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' library(MASS)
#' countMat<-t(sapply(vector("numeric",length = length(geneLengthGRCh38)),function(x){
#' 	rnegbin(10,theta = abs(rnorm(1,mean = 10,sd = 20)),mu = abs(rnorm(1,mean = 10,sd = 20)))
#' }));rownames(countMat)<-names(geneLengthGRCh38)
#' colnames(countMat)<-letters[1:ncol(countMat)]
#' computeQCmetricSamples(countMat)
computeQCmetricSamples <- function(x, uncenter = FALSE) {
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
	for (i in 1:ncol(x))
		noGenEx[i] <- length(which(x[, i] > zero))

	return(data.frame(
		mean = mean,
		sd = sd,
		CV = CV,
		TotalGenEx = noGenEx,
		TotalCount = count
	))
}


#' Principal Component Analysis
#'
#' @param d A matrix of numeric (in the RNA-Seq context, log counts).
#' @param transpose Logical. If `transpose`, samples are columns and features are rows.
#' @param scale Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#'
#' @return
#' A list with the following element:
#' - sdev: the standard deviations of the principal components (i.e., the square roots of the eigenvalues of the covariance/correlation matrix, though the calculation is actually done with the singular values of the data matrix).
#' - rotation: the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors). The function princomp returns this in the element loadings.
#' - x: if retx is true the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned. Hence, cov(x) is the diagonal matrix diag(sdev^2). For the formula method, napredict() is applied to handle the treatment of values omitted by the na.action.
#' - center, scale: the centering and scaling used, or FALSE.
#' - n.obs: Number of samples
#' - propExplVar: proportion of explained variance by principals components on total variance.
#' - transform: a list containing the scaling (sdeviations) and centering factors (means) for each principal component.
#' @export
#'
#' @examples
#' data(iris)
#' pca<-PCA(iris[,1:4],transpose = FALSE,scale = TRUE,center = TRUE)
PCA <- function(d,
								transpose = T,
								scale = F,
								center = T) {
	if (transpose)
		d <- t(d)

	means <- 0
	sdeviations <- 1
	if (center) {
		means <- apply(d, 2, mean)
		d <- sweep(d, 2, means, "-")
	}
	if (scale) {
		sdeviations <- apply(d, 2, sd)
		d <- sweep(d, 2, sdeviations, "/")
	}
	resacp <- prcomp(
		x = d,
		retx = T,
		center = FALSE,
		scale = FALSE
	)

	resacp$n.obs <- dim(d)[1]

	resacp$propExplVar <- resacp$sdev ^ 2 / sum(resacp$sdev ^ 2)
	resacp$scale <- scale
	resacp$center <- center
	resacp$transform <- list(sdeviations = sdeviations, means = means)
	resacp$isFastPCA<-FALSE
	return(resacp)

}


#' Faster Principal Component Analysis
#'
#' @param d A matrix of numeric (in the RNA-Seq context, log counts).
#' @param transpose Logical. If `transpose`, samples are columns and features are rows.
#' @param scale Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param nPC Integer. Number of Principal Component to be computed.
#' @param weight.by.var Logical. Weight PC by variables. If TRUE return a regular PCA.
#' @param ... Other parameters passed to `irlba`.
#'
#' @return
#' A list with the following element:
#' - sdev: the standard deviations of the principal components (i.e., the square roots of the eigenvalues of the covariance/correlation matrix, though the calculation is actually done with the singular values of the data matrix).
#' - rotation: the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors). The function princomp returns this in the element loadings.
#' - x: if retx is true the value of the rotated data (the centred (and scaled if requested) data multiplied by the rotation matrix) is returned. Hence, cov(x) is the diagonal matrix diag(sdev^2). For the formula method, napredict() is applied to handle the treatment of values omitted by the na.action.
#' - center, scale: the centering and scaling used, or FALSE.
#' - n.obs: Number of samples
#' - propExplVar: proportion of explained variance by principals components on total variance.
#' - transform: a list containing the scaling (sdeviations) and centering factors (means) for each principal component.
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' pca<-fastPCA(bulkLogCounts,nPC=10)
#' pca2d(pca)
fastPCA <-
	function(d,
					 transpose = TRUE,
					 scale = FALSE,
					 center = TRUE,
					 nPC = min(ncol(d) - 1, nrow(d) - 1, 30),
					 weight.by.var = TRUE,
					 ...) {
		if (transpose)
			d <- t(d)

		d <- as.matrix(d)
		means <- 0
		sdeviations <- 1
		if (center | scale)
			d <- scale(d, scale = scale, center = center)
		if (center)
			means <- attr(d, "scaled:center")
		if (scale)
			sdeviations <- attr(d, "scaled:scale")

		resacp <- list()
		resacp$n.obs <- dim(d)[1]

		resacp$scale <- scale
		resacp$center <- center
		resacp$transform <- list(sdeviations = sdeviations, means = means)

		irlbaResults <- irlba::irlba(A = d, nv = nPC, ...)
		rotation <- irlbaResults$v
		resacp$sdev <- irlbaResults$d / sqrt(max(1, nrow(d) - 1))
		if (weight.by.var) {
			if (nPC > 1) {
				reducedSpace <- irlbaResults$u %*% diag(irlbaResults$d)
			} else {
				reducedSpace <- irlbaResults$u %*% irlbaResults$d
			}
		} else {
			reducedSpace <- irlbaResults$u
		}
		rownames(rotation) <- colnames(d)
		colnames(rotation) <- paste0("PC", 1:nPC)
		rownames(reducedSpace) <- rownames(d)
		colnames(reducedSpace) <- colnames(rotation)
		resacp$x <- reducedSpace
		resacp$rotation <- rotation
		resacp$propExplVar <- resacp$sdev ^ 2 / sum(resacp$sdev ^ 2)
		resacp$isFastPCA<-TRUE
		resacp
	}

#' Add new samples to an existing PCA.
#'
#' @param pca The existing PCA as an object returned by `PCA` or `fastPCA`.
#' @param newSamplesMatrix The matrix of numeric of new samples. It must have the same features than the original PCA.
#' @param transpose Logical. If `transpose`, samples are columns and features are rows.
#' @param returnPCA Logical. Return the updated PCA. If false, just returned the coordinates of new samples on PCs.
#'
#' @return A PCA list if `returnPCA` or a matrix of coordinates of samples on PCs.
#' @export
#'
#' @examples
#' data("iris")
#' iris1<-iris[1:75,]
#' iris2<-iris[76:150,]
#' pca<-PCA(iris1[,1:4],transpose = FALSE,scale = TRUE,center = TRUE)
#' pca2d(pca,colorBy = iris1$Species)
#' pcaUpdated<-pcaAddSamples(pca,iris2[,1:4],transpose = FALSE)
#' pca2d(pcaUpdated,colorBy = iris$Species)

pcaAddSamples <-
	function(pca,
					 newSamplesMatrix,
					 transpose = TRUE,
					 returnPCA = TRUE) {
		if (transpose)
			newSamplesMatrix <- t(newSamplesMatrix)
		newSamplesMatrix <-
			newSamplesMatrix[, rownames(pca$rotation), drop = FALSE]

		if (pca$center)
			newSamplesMatrix <- sweep(newSamplesMatrix, 2, pca$transform$means, "-")
		if (pca$scale)
			newSamplesMatrix <-
			sweep(newSamplesMatrix, 2, pca$transform$sdeviations, "/")

		newSamplesCoord <- t(apply(newSamplesMatrix, 1, function(x) {
			colSums(x * pca$rotation)
		}))
		if (returnPCA) {
			pca$x <- rbind(pca$x, newSamplesCoord)
			return(pca)
		} else{
			return(newSamplesCoord)
		}
	}


#' Principal Component Regression
#'
#' @description PCR is a way to link experimental variable to principal components.
#'
#' @param pca A PCA as an object returned by `PCA` or `fastPCA`.
#' @param annotationDF A dataframe of feature (numeric or factor) with rows as samples. Must have the same number of samples than the PCA.
#' @param nComponent Integer, number of PC used in the regression.
#'
#' @return A dataframe with the PC, annotation and the corresponding R squared.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#'
#' pca<-fastPCA(bulkLogCounts,nPC=10)
#' PCR(pca,annotationDF=sampleAnnot[,c("culture_media","line","passage")])
PCR <- function(pca, annotationDF, nComponent = 10) {
	annots <- colnames(annotationDF)
	rSquaredMat <- matrixFromDimnames(cn(pca$x), annots, value = NA)
	for (x in cn(pca$x)) {
		for (annot in annots) {
			rSquaredMat[x, annot] <-
				summary(lm(formula(paste0(x, "~", annot)), data = data.frame(annotationDF[, annot, drop =
																																										F], pca$x[, x, drop = F])))$r.squared
		}
	}
	if (ncol(annotationDF) < 2) {
		retDt <-
			data.frame(
				Rsquared = rSquaredMat[, 1],
				PC = rownames(rSquaredMat),
				Annotation = colnames(annotationDF)
			)
	} else{
		retDt <-
			reshape2::melt(rSquaredMat[1:nComponent, ],
										 value.name = "Rsquared",
										 varnames = c("PC", "Annotation"))
	}

	retDt$PC <-
		factor(retDt$PC, levels = cn(pca$x)[1:nComponent]) #so the levels of PCs are well ordered
	retDt
}

merge0dist <- function(disMat) {
	mat <- as.matrix(disMat)
	merged <- list()
	found <- TRUE
	while (found == TRUE) {
		found <- FALSE
		for (i in 2:nrow(mat)) {
			for (j in 1:(i - 1)) {
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

#' Non-metric Multidimensional Scaling (dimension reduction)
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#' @param transpose Logical. If `transpose`, samples are columns and features are rows.
#' @param scale  Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param metric A function that return a object of class "dist".
#' @param ndim Integer. Dimensions of the embedded space.
#' @param maxit The maximum number of iterations.
#'
#' @return A matrix of coordinates with samples as rows.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' NMDSproj<-NMDS(bulkLogCounts)
#' proj2d(NMDSproj)
NMDS <-
	function(data,
					 transpose = TRUE,
					 scale = FALSE,
					 center = FALSE,
					 metric = dist,
					 ndim = 2,
					 maxit = 100) {
		merged <- FALSE
		if (transpose)
			data <- t(data)
		d <- metric(data)  # euclidean distances between the rows
		if (min(d, na.rm = TRUE) == 0) {
			merged <- TRUE
			md <- merge0dist(d)
			d <- md$distMat
			mergedSample <- md$merged
		}
		fit <-
			MASS::isoMDS(d, k = ndim, maxit = maxit) # k is the number of dim
		fit$coord <- fit$points
		fit$points <- NULL
		if (merged) {
			for (sple in names(mergedSample)) {
				values <-
					matrix(rep(fit$coord[sple, ], length(mergedSample[[sple]])),
								 nrow = length(mergedSample[[sple]]),
								 byrow = TRUE)
				rownames(values) <- mergedSample[[sple]]
				fit$coord <- rbind(fit$coord, values)
			}
		}
		return(fit$coord)
	}

#' Count Per Million normalization
#'
#' @param data A raw count table with rows as genes
#'
#' @return A normalized count table.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' library(MASS)
#' countMat<-t(sapply(vector("numeric",length = length(geneLengthGRCh38)),function(x){
#' 	rnegbin(10,theta = abs(rnorm(1,mean = 10,sd = 20)),mu = abs(rnorm(1,mean = 10,sd = 20)))
#' }));rownames(countMat)<-names(geneLengthGRCh38)
#' CPM(countMat)
CPM <- function(data) {
	#Normalisation CPM
	data.CPM <- sweep(data, 2, colSums(data), `/`)
	data.CPM <- data.CPM * 1000000
	return(data.CPM)
}


#' Transcript per milion (TPM) normalization for full-length transcript, short read sequencing.
#'
#' @param data A raw count table with rows as genes
#' @param gene.length A vector of numeric corresponding to gene length in base pair. Must be named by genes.
#'
#' @return A normalized count table where 1 count is equal in theory to one transcript per milion.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' library(MASS)
#' countMat<-t(sapply(vector("numeric",length = length(geneLengthGRCh38)),function(x){
#' 	rnegbin(10,theta = abs(rnorm(1,mean = 10,sd = 20)),mu = abs(rnorm(1,mean = 10,sd = 20)))
#' }));rownames(countMat)<-names(geneLengthGRCh38)
#' TPMfullLength(countMat,geneLengthGRCh38)
TPMfullLength <- function(data, gene.length) {
	gene.length.kb <- gene.length[rownames(data)] / 1000
	data <- sweep(data, 1, gene.length.kb, `/`)
	return(CPM(data))
}

#' Reads Per Kilobase per Million (RPKM) normalization for full-length transcript, short read sequencing.
#'
#' @param data A raw count table with rows as genes
#' @param gene.length A vector of numeric corresponding to gene length in base pair. Must be named by genes.
#'
#' @return A normalized count table of RPKM.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' library(MASS)
#' countMat<-t(sapply(vector("numeric",length = length(geneLengthGRCh38)),function(x){
#' 	rnegbin(10,theta = abs(rnorm(1,mean = 10,sd = 20)),mu = abs(rnorm(1,mean = 10,sd = 20)))
#' }));rownames(countMat)<-names(geneLengthGRCh38)
#' RPKM(countMat,geneLengthGRCh38)
RPKM <- function(data, gene.length) {
	gene.length.kb <- gene.length[rn(data)] / 1000
	data <- CPM(data)
	sweep(data, 1, gene.length.kb, `/`)
}


#' Quick DESeq2 normalization
#'
#' @param countMatrix A raw count table with rows as genes
#'
#' @return A normalized count table.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' library(MASS)
#' countMat<-t(sapply(vector("numeric",length = length(geneLengthGRCh38)),function(x){
#' 	rnegbin(10,theta = abs(rnorm(1,mean = 10,sd = 20)),mu = abs(rnorm(1,mean = 10,sd = 20)))
#' }));rownames(countMat)<-names(geneLengthGRCh38)
#' normDeseq(countMat)
normDeseq <-
	function(countMatrix) {
		#matrix where genes are rows and samples are columns
		# PS = pseudo reference sample
		PS <-
			apply(countMatrix, 1, gmean, keepZero = TRUE) #get a vector which consist of the geometrical mean of each genes across all samples
		keptRow <- PS > 0 #get rid of genes containing one zero ore more
		PS <- PS[keptRow]
		ratioMat <-
			sweep(countMatrix[keptRow, ], 1, PS, "/") #get the ratio matrix (expression/expression from PS)
		normFactors <-
			apply(ratioMat, 2, median) #get the median of the ratios for each sample to get the normalization factors
		sweep(countMatrix, 2, normFactors, "/") #divide each sample by the corresponding normalization factor
	}


#' Quick single cell normalization (scran method).
#'
#' @param rawCounts RNA-Seq raw counts with rows as genes.
#' @param returnLog Logical. Return log counts.
#' @param sizeFactors NULL or a vector of numeric containing precomputed size factor, same size as number of cells in `rawCounts`.
#' @param ... Other parameter passed to `computeSumFactors`.
#'
#' @return A normalized count table.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' library(MASS)
#' countMat<-t(sapply(vector("numeric",length = length(geneLengthGRCh38)),function(x){
#' 	rnegbin(10,theta = abs(rnorm(1,mean = 10,sd = 20)),mu = abs(rnorm(1,mean = 10,sd = 20)))
#' }));rownames(countMat)<-names(geneLengthGRCh38)
#' quickSCnorm(countMat, returnLog=FALSE)
#' quickSCnorm(countMat, returnLog=TRUE)
quickSCnorm <-
	function(rawCounts,
					 returnLog = TRUE,
					 sizeFactors = NULL,
					 ...) {
		sce <-
			SingleCellExperiment::SingleCellExperiment(assays = list(counts = rawCounts))
		if (!is.null(sizeFactors)) {
			sizeFactors(sce) <- sizeFactors
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


#' Determine the best partition in a hierarchichal clustering
#'
#' @param hc A hclust object.
#' @param min Minimum number of class in the partition.
#' @param max Maximum number of class in the partition
#' @param loss Logical. Return the list of computed partition with their derivative loss.
#' @param graph Logical. Plot a graph of computed partition with their derivative loss.
#'
#' @return A single integer (best partition) or print a graph if `graph` or a vector of numeric if `loss`.
#' @export
#'
#' @examples
#' data(iris)
#' resClust<-hierarchicalClustering(iris[,1:3],transpose = FALSE)
#' best.cutree(resClust,graph=TRUE)
#' best.cutree(resClust)
#' best.cutree(resClust, loss=TRUE)
#' cutree(resClust,k = best.cutree(resClust))
best.cutree <- function(hc,
												min = 2,
												max = 20,
												loss = FALSE,
												graph = FALSE) {
	if (class(hc) != "hclust")
		hc <- as.hclust(hc)
	max <- min(max, length(hc$height) - 1)
	inert.gain <- rev(hc$height)
	intra <- rev(cumsum(rev(inert.gain)))
	relative.loss = intra[min:(max + 1)] / intra[(min - 1):(max)]
	derivative.loss = relative.loss[2:length(relative.loss)] - relative.loss[1:(length(relative.loss) -
																																								1)]
	names(derivative.loss) <- min:max
	if (graph) {
		print(
			ggplot(
				data.frame(partition = min:max, derivative.loss = derivative.loss),
				aes(x = partition, y = derivative.loss)
			) +
				geom_point() +
				scale_x_continuous(breaks = min:max, minor_breaks = NULL)
		)
	} else {
		if (loss)
			derivative.loss
		else
			as.numeric(names(which.max(derivative.loss)))
	}
}


#' Perform a hierarchical clustering from a matrix/df of observation × features
#'
#' @param x A matrix or dataframe of numeric.
#' @param transpose Logical. If `transpose`, samples are columns and features are rows.
#' @param method.dist A method from the "dist" function. Can be also "pearson" or "bicor" for correlation distance.
#' @param method.hclust the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param bootstrap Logical. Use bootstrapping for determining the best clustering.
#' @param nboot The number of bootstrap replications.
#' @param PCAfirst Compute a PCA before computing teh clustering, if x contains a lot of features, it can reduce computation time.
#' @param nDimPCA Integer. If `PCAfirst`, compute a PCA first and take n first principal components.
#'
#' @return A hclust object.
#' @export
#'
#' @examples
#' data(iris)
#' resClust<-hierarchicalClustering(iris[,1:3],transpose = FALSE)
#' plot(resClust,hang=-1)
#' resClust<-hierarchicalClustering(iris[,1:3],transpose = TRUE,bootstrap = TRUE,nboot = 20)
#' plot(resClust,hang=-1)
hierarchicalClustering <-
	function(x,
					 transpose = TRUE,
					 method.dist = "euclidean",
					 method.hclust = "ward.D2",
					 bootstrap = FALSE,
					 nboot = 10,
					 PCAfirst = FALSE,
					 nDimPCA = NULL) {
		if (transpose)
			x <- t(x)
		if (PCAfirst) {
			x <- PCA(x, transpose = FALSE, scale = FALSE)$x
			if (!is.null(nDimPCA)) {
				x <- x[, 1:nDimPCA]
			}
		}
		if (bootstrap) {
			resClust <-
				pvclust::pvclust(
					t(x),
					nboot = nboot,
					method.hclust = method.hclust,
					parallel = TRUE,
					method.dist = method.dist
				)$hclust
		} else{
			if (method.dist == "pearson") {
				resDist <- corrDist(x)
			} else if (method.dist == "bicor") {
				resDist <-
					as.dist((1 - suppressWarnings(WGCNA::bicor(Matrix::t(
						x
					)))) / 2)
			} else{
				resDist <- dist(x, method = method.dist)
			}
			resClust <- stats::hclust(resDist, method = method.hclust)
		}
		return(resClust)
	}

#' Test a linear model on each gene following an experimental design.
#'
#' @param exprData A matrix of numeric with rows as features (in the RNA-Seq context, log counts).
#' @param colData A dataframe of feature (numeric or factor) with rows as samples. Must have the same number of samples than exprData
#' @param contrast A vector of 3 character.
#' 1. Name of the experimental variable that have to be used for differential activation. Must be a column name of `colData`.
#' 2. Condition considered as the reference.
#' 3. Condition considered as the target group.
#'
#' @return
#' A dataframe with the following columns:
#' - baseMean: mean
#' - log2FoldChange: Log(Log Fold Change)  between the two tested groups.
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' res<-multiLinearModel(bulkLogCounts,colData = sampleAnnot,contrast = c("culture_media","T2iLGO","KSR+FGF2"))
multiLinearModel <- function(exprData, colData, designFormula ,contrast) {
	samples<-rownames(colData)[colData[,contrast[1]]%in%contrast[2:3]]
	data<-exprData[,samples]
	groups<-droplevels(as.factor(colData[samples,contrast[1]]))
	logicGroup<-rep(F,len(groups))
	logicGroup[groups==contrast[2]]<-T
	regTabList<-apply(data,1,function(x){
		data.frame(data=x,group=logicGroup)
	})
	resList<-lapply(regTabList,function(regTab){
		summary(lm(data ~ group,data=regTab))$coefficients[2,c(1,4)]
	})
	res<-data.frame(do.call("rbind",resList));colnames(res)<-c("log2FoldChange","pval")
	res<-cbind(data.frame(baseMean=apply(exprData[,samples],1,mean)),res)
	res$padj<-p.adjust(res$pval,method = "BH")
	return(res)
}


#' Compute over dispersion values for each gene.
#'
#' @param counts Normalized count table with genes as rows.
#' @param minCount Minimum average expression to not be filtered out.
#' @param plot Logical. Show the overdispersion plot.
#' @param returnPlot Logical, if `plot` return it as a ggplot object instead of printing it.
#'
#' @return A ggplot graph if `returnPlot`, otherwise a dataframe with the following columns:
#' - mu: average expression
#' - var: variance
#' - cv2: squared coefficient of variation. Used as a dispersion value.
#' - residuals: y-distance from teh regression. Can be used as an overdispersion value.
#' - residuals2: squared residuals
#' - fitted: theoretical dispersion for the gene average (y value of the curve).
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' normCount<-2^(bulkLogCounts-1)
#' dispData<-getMostVariableGenes(normCount,minCount=1)
getMostVariableGenes <-
	function(counts,
					 minCount = 0.01,
					 plot = TRUE,
					 returnPlot = FALSE) {
		counts <- counts[rowMeans(counts) > minCount, ]
		dispTable <-
			data.frame(
				mu = rowMeans(counts),
				var = apply(counts, 1, var),
				row.names = rownames(counts)
			)
		dispTable$cv2 <- dispTable$var / dispTable$mu ^ 2
		sumNullvariance <- sum(dispTable$cv2 <= 0)
		if (sumNullvariance > 0) {
			warning(paste0(sumNullvariance, " have null variance and will be removed"))
			dispTable <- dispTable[dispTable$cv2 > 0, ]
		}
		fit <- loess(cv2 ~ mu, data = log10(dispTable[, c("mu", "cv2")]))
		dispTable$residuals <- fit$residuals
		dispTable$residuals2 <- dispTable$residuals ^ 2
		dispTable$fitted <- 10 ^ fit$fitted
		if (plot) {
			g <-
				ggplot(dispTable,
							 aes(
							 	x = mu,
							 	y = cv2,
							 	label = rownames(dispTable),
							 	fill = residuals
							 )) +
				geom_point(stroke = 1 / 8,
									 colour = "black",
									 shape = 21) + geom_line(aes(y = fitted), color = "red", size = 1.5) +
				scale_x_log10() + scale_y_log10()
			if (returnPlot) {
				return(g)
			} else{
				print(g)
			}
		}
		dispTable
	}

#' Compute over dispersion values for each gene from log counts. Do not used, need to be fixed
#'
#' @param logCounts Normalized log count table with genes as rows.
#' @param minCount Minimum average expression to not be filtered out.
#' @param plot Logical. Show the overdispersion plot.
#' @param returnPlot Logical, if `plot` return it as a ggplot object instead of printing it.
#'
#' @return A ggplot graph if `returnPlot`, otherwise a dataframe with the following columns:
#' - mu: average expression
#' - var: variance
#' - cv2: squared coefficient of variation. Used as a dispersion value.
#' - residuals: y-distance from the regression. Can be used as an overdispersion value.
#' - residuals2: squared residuals
#' - fitted: theoretical dispersion for the gene average (y value of the curve).
#' not export
#'
#' @examples
#' data("bulkLogCounts")
#' dispDataLog<-getMostVariableGenes(bulkLogCounts,minCount=1)
getMostVariableGenesLogCount <-
	function(logCounts,
					 minCount = 0,
					 plot = TRUE,
					 returnPlot = FALSE) {
		if (minCount > 0)
			logCounts <- logCounts[rowMeans(logCounts) > minCount, ]
		dispTable <-
			data.frame(
				mu = rowMeans(logCounts),
				var = apply(logCounts, 1, var),
				row.names = rownames(logCounts)
			)

		fit <- loess(var ~ mu, data = dispTable)
		if (length(rownames(fit$x)) < nrow(dispTable)) {
			warning(nrow(dispTable) - length(rownames(fit$x)),
							" genes have null variance and will be removed")
		}
		dispTable <- dispTable[rownames(fit$x), ]
		dispTable$residuals <- fit$residuals
		dispTable$fitted <- fit$fitted

		if (plot) {
			g <-
				ggplot(dispTable,
							 aes(
							 	x = mu,
							 	y = var,
							 	label = rownames(dispTable),
							 	fill = residuals
							 )) +
				geom_point(stroke = 1 / 8,
									 colour = "black",
									 shape = 21) + geom_line(aes(y = fitted), color = "red", size = 1.5)
			if (returnPlot) {
				return(g)
			} else{
				print(g)
			}
		}
		dispTable
	}


#' Approximation of area under ROC curve
#'
#' @description By Miron Kursa https://mbq.me
#'
#' @param score a vector of numeric representing the measure of a feature in a set of samples.
#' @param boolVect a vector of logical, same size as score. Is the sample in the target group?
#'
#' @return A numeric value. Close to 1 = perfect marker, around 0.5 = as good as random values, close to 0 = perfect anti-marker.
#' @export
#'
#' @examples
#' data(iris)
#' auroc(iris$Sepal.Length,iris$Species=="virginica")
auroc <- function(score, boolVect) {
	n1 <- sum(!boolVect)
	n2 <- sum(boolVect)
	U	<- sum(rank(score)[!boolVect]) - n1 * (n1 + 1) / 2
	return(1 - U / n1 / n2)
}

#' Quick approximation of area under ROC curve
#'
#' @description By Miron Kursa https://mbq.me
#'
#' @param score a vector of numeric representing the measure of a feature in a set of samples.
#' @param boolVect a vector of logical, same size as score. Is the sample in the target group?
#' @param n1 a single integer, number of observation not in the target group.
#' @param n2 a single integer, number of observation in the target group.
#'
#' @return A numeric value. Close to 1 = perfect marker, around 0.5 = as good as random values, close to 0 = perfect anti-marker.
#' @export
#'
#' @examples
#' data(iris)
#' qauroc(iris$Sepal.Length,iris$Species=="virginica", sum(!iris$Species=="virginica"),sum(iris$Species=="virginica"))
qauroc <- function(score, boolVect) {
	n1 <- sum(!boolVect)
	n2 <- sum(boolVect)
	aurocCPP(as.numeric(score),
					 as.logical(boolVect),
					 as.integer(n1),
					 as.integer(n2))
}

#bug, if group contain only 2 col

#' Compute a dataframe with marker metrics describing best gene marker per group of samples.
#'
#' @param expressionMatrix A matrix of numeric with rows as features (in the RNA-Seq context, log counts).
#' @param group A feature of factor/character, same length as number of sample. Describe group of each sample.
#' @param useBiocParallel A logical describing if the function should be launched in multi-threading mode.
#' @param BPPARAM A BPPARAM object as return by `BiocParallel::bpparam()`. Use for multi-threading.
#' @param returnAsList Return a list where each element is dataframe containing the marker metrics of a group.
#'
#' @return A dataframe containing four column per group: Log2(Fold-Change), AUROC, marker score (see details), p-value and BH adjusted p-value.
#'
#' @details
#' LogFC and pvalues are computed from a linear modelling of the data.
#'
#' Score is consisting of the geometrical mean of absolute LogFC, absolute(auroc - 0.5), and -log10(pval), then signed by the logFC:
#' score = sign(logFC) × gmean( abs(logFC), abs(aurocRes-0.5), -log10(pval) )
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media)
#'
getMarkers <-
	function (exprData,
						groups,
						multiThread = FALSE,
						BPPARAM = NULL,
						returnAsList = FALSE) {

		if (is.null(BPPARAM))
			BPPARAM = BiocParallel::bpparam()
		BiocParallel::register(BPPARAM)
		if (multiThread) {
			doParallel::registerDoParallel(cores = BPPARAM$workers)
		}else{
			doParallel::registerDoParallel(cores = 1)
		}

		groups <- as.factor(make.names(groups))
		grpLvl <- levels(groups)

		resDF <- foreach(group = grpLvl) %dopar% {
			logicGroup <- groups == group
			designMat <- cbind(1, logicGroup)
			n1 = as.integer(sum(!logicGroup))
			n2 = as.integer(sum(logicGroup))
			resDFgroup <- apply(exprData, 1, function(gene) {
				aurocRes <- aurocCPP(
					score = gene,
					boolVect = designMat[,
															 2],
					n1 = n1,
					n2 = n2
				)
				lmOut <- .lm.fit(designMat, gene)
				pval <- pvalLmFit(lmOut$residuals,
													lmOut$coefficients,
													p = lmOut$rank,
													qr = lmOut$qr)[2]
				coef <- lmOut$coefficients[2]
				score <- sign(coef)*prod(c(abs(coef),abs(aurocRes-0.5), min(-log10(pval),324)))^(1/3) #min(-log10(pval),324): avoid Inf
				return(c(coef, aurocRes, score, pval))
			}) |> t() |> data.frame()
			colnames(resDFgroup) <- c("lfc", "auroc", "score", "pval")
			resDFgroup$padj <- p.adjust(resDFgroup$pval, method = "BH")
			resDFgroup
		}
		names(resDF) <- grpLvl

		if (returnAsList) {
			return(resDF)
		} else{
			return(do.call("cbind", resDF))
		}

	}

#' Extract a specific feature/metric (pval, logFC...) from a marker result dataframe
#'
#' @param markerData A dataframe returned by `getMarkers`.
#' @param feature The name of the feature that has to be extracted.
#'
#' @return
#' A matrix containing only the wanted feature where each column is a group.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' markerData <- getMarkers(bulkLogCounts,sampleAnnot$culture_media)
#' markerData <- extractFeatureMarkerData(markerData)
extractFeatureMarkerData<-function(markerData, feature="score"){
	featureStrLen<-nchar(feature)+1
	columns2Pick<-grep(paste0("^.*\\.",feature,"$"),colnames(markerData),value = TRUE)
	markerData<-markerData[,columns2Pick] |> as.matrix()
	grpName<-colnames(markerData)
	colnames(markerData)<-substr(grpName,1,nchar(grpName)-featureStrLen)
	return(markerData)
}

#' Compute a matrix of coef describing best gene marker per group of samples from a GLM net regression
#'
#' @param expressionMatrix A matrix of numeric with rows as features (in the RNA-Seq context, log counts).
#' @param group A feature of factor/character, same length as number of sample. Describe group of each sample.
#'
#' @return A matrix containing each gene coefficient for each group.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#' res <- getMarkerGLMnet(bulkLogCounts,sampleAnnot$culture_media)
getMarkerGLMnet <- function(expressionMatrix, group) {
	geneMatrix <- as.matrix(expressionMatrix) |> t()
	fit <- glmnet(
		geneMatrix,
		group,
		family = "multinomial",
		alpha = .5,
		lambda = cv.glmnet(geneMatrix, group, family = "multinomial")$lambda.1se
	)

	cf <- sapply(coef(fit), function(x)
		x[2:length(x)])
	rownames(cf) <- colnames(geneMatrix)
	cf
}

#' Correlation from one gene to all others.
#'
#' @param gene A single charachter. The gene (or feature) name that will be used for correlating to all others.
#' @param expression A matrix of numeric with rows as features (in the RNA-Seq context, log counts).
#' @param corFun A function to compute a correlation between two feature.
#' @param ... Parameters passed to `corFun`.
#'
#' @return A vector of numeric. Correlations values named by their corresponding gene.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' corGeneToOthers("NANOG",bulkLogCounts)
corGeneToOthers <- function(gene, expression, corFun = cor, ...) {
	expression <- as.matrix(expression)
	t(corFun(expression[gene, ], t(expression), ...))[, 1]
}


#' Compute the activation score of a gene set from 1st component of its PCA
#'
#' @param exprMatrix A matrix of numeric with rows as features (in the RNA-Seq context, log counts).
#' @param genes A character vector. The gene set where the activation score has to be computed. Must be a subset of `exprMatrix` row names.
#' @param scale  Logical. Divide features by their standard deviation.
#' @param center Logical. Subtract features by their average.
#' @param returnContribution Logical. Return list with activation score and contribution of genes to the activation score.
#'
#' @return A vector of numeric corresponding to activation scores, named by genes. If `returnContribution` return a list with activation scores and contributions of genes.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' keggData<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSet<-keggData$kegg$`hsa00190 Oxidative phosphorylation`
#' geneSet<-intersect(geneSet,rownames(bulkLogCounts))
#' activScorePCA(bulkLogCounts,genes = geneSet)
#' activScorePCA(bulkLogCounts,genes = geneSet,returnContribution = TRUE)
activScorePCA <-
	function(exprMatrix,
					 genes,
					 scale = FALSE,
					 center = TRUE,
					 returnContribution = FALSE) {
		pca <-
			fastPCA(exprMatrix[genes, ],
							center = center,
							scale = scale,
							nPC = 1)
		activScore <- pca$x[, 1]
		contribution <- pca$rotation[, 1]
		if (cor(colMeans(exprMatrix[genes, ]), activScore) < 0) {
			activScore <- -activScore
			contribution <- -contribution
		}
		if (returnContribution) {
			list(activScore = activScore, contribution = contribution)
		} else{
			activScore
		}
	}

#' Compute PCA activation from a gene list.
#'
#' @param exprMatrix A matrix of numeric with rows as features (in the RNA-Seq context, log counts).
#' @param geneList A named list containing character vectors of genes.
#' @param scale A character vector. The gene set where the activation score has to be computed. Must be a subset of `exprMatrix` row names.
#' @param center Logical. Subtract features by their average.
#'
#' @return A list of 2 element: `activScoreMat`: The activation score matrix. `contributionList`: A list of named numeric containing contributions of genes.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' keggData<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' res<-activeScorePCAlist(bulkLogCounts,geneList = keggData$kegg[1:3])
#' activeScoreMat<-res$activScoreMat
#' contributionList<-res$contributionList
activeScorePCAlist <-
	function(exprMatrix,
					 geneList,
					 scale = FALSE,
					 center = TRUE) {
		res <-
			lapply(listGenePerRegulon, function(genesOfReg)
				activScorePCA(
					exprMatrix,
					genesOfReg,
					returnContribution = TRUE,
					scale = scale,
					center = center
				))
		list(
			activScoreMat = sapply(res, function(x)
				x$activScore),
			contributionList = lapply(res, function(x)
				x$contribution)
		)
	}

#' UMAP projection
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#' @param nDimPCA Integer. If not NULL compute a PCA first and take n first princpal components.
#' @param transpose Logical. If `transpose`, samples are columns and features are rows.
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
#' Larger values result in more global views of the manifold, while smaller values result in more local data being preserved.
#' In general values should be in the range 2 to 100.
#' @param n_components Integer. Dimensions of the embedded space.
#' @param min_dist The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out.
#' @param init Type of initialization for the coordinates (see `uwot` for more details).
#' @param metric Type of distance metric to use to find nearest neighbors (see `uwot` for more details).
#' @param ret_model If TRUE, then return extra data that can be used to add new data to an existing embedding via umap_transform. The embedded coordinates are returned as the list item embedding. If FALSE, just return the coordinates. This parameter can be used in conjunction with ret_nn and ret_extra. Note that some settings are incompatible with the production of a UMAP model: external neighbor data (passed via a list to nn_method), and factor columns that were included via the metric parameter. In the latter case, the model produced is based only on the numeric data. A transformation using new data is possible, but the factor columns in the new data are ignored
#' @param ret_nn If TRUE, then in addition to the embedding, also return nearest neighbor data that can be used as input to nn_method to avoid the overhead of repeatedly calculating the nearest neighbors when manipulating unrelated parameters (e.g. min_dist, n_epochs, init). See the "Value" section for the names of the list items. If FALSE, just return the coordinates. Note that the nearest neighbors could be sensitive to data scaling, so be wary of reusing nearest neighbor data if modifying the scale parameter. This parameter can be used in conjunction with ret_model and ret_extra.
#' @param ... Other parameters passed to uwot.
#'
#' @return A matrix of coordinates with samples as rows, or:
#' - if ret_model = TRUE (or ret_extra contains "model"), returns a list containing extra information that can be used to add new data to an existing embedding via umap_transform. In this case, the coordinates are available in the list item embedding. NOTE: The contents of the model list should not be considered stable or part of the public API, and are purposely left undocumented.
#' - if ret_nn = TRUE (or ret_extra contains "nn"), returns the nearest neighbor data as a list called nn. This contains one list for each metric calculated, itself containing a matrix idx with the integer ids of the neighbors; and a matrix dist with the distances. The nn list (or a sub-list) can be used as input to the nn_method parameter.
#' @export
#'
#' @examples
#' data(iris)
#' irisUMAP<-UMAP(iris[,1:4],transpose = FALSE)
#' proj2d(irisUMAP,colorBy = iris$Species)
#' irisUMAP<-UMAP(rowScale(iris[,1:4],center = TRUE,scaled = TRUE),transpose = FALSE,n_neighbors = nrow(iris),ret_nn = TRUE)
#' proj2d(irisUMAP$embedding,colorBy = iris$Species,nnMatrix = irisUMAP$nn$euclidean$idx[,1:3],fixedCoord = TRUE)
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
					 ...) {
		if (transpose)
			data <- t(data)
		if (nrow(data) < n_neighbors) {
			n_neighbors <- nrow(data)
			warning(
				"n_neighbors must not exceed number of samples, adjusting n_neighbors to number of samples (",
				n_neighbors,
				")"
			)
		}
		if (is.null(n_neighbors))
			n_neighbors = nrow(data)
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
		res
	}


#' TriMap dimension reduction.
#'
#' @param data A matrix of numeric (in the RNA-Seq context, log counts).
#' @param n_dims Integer. Dimensions of the embedded space.
#' @param transpose Logical. If `transpose`, samples are columns and features are rows.
#' @param n_inliers Number of nearest neighbors for forming the nearest neighbor triplets.
#' @param n_outliers Number of outliers for forming the nearest neighbor triplets.
#' @param apply_pca Reduce the number of dimensions of the data to 100 if necessary before applying the nearest-neighbor search.
#' @param n_iters Number of iterations.
#' @param knn_tuple Use the precomputed nearest-neighbors information in form of a tuple (knn_nbrs, knn_distances).
#'
#' @return A matrix of coordinates with samples as rows.
#' @export
#'
#' @examples
#' data(iris)
#' irisTrimap<-TRIMAP(iris[,1:4],transpose=FALSE)
#' proj2d(irisTrimap,colorBy = iris$Species)
TRIMAP <-
	function(data,
					 n_dims = 2,
					 transpose = TRUE,
					 n_inliers = 10,
					 n_outliers = 5,
					 apply_pca = TRUE,
					 n_iters = 400,
					 knn_tuple = NULL) {
		if (transpose)
			data <- t(data)
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
	}

#' Compute a Leiden clustering from a UMAP model.
#'
#' @param umapWithNN A list with the UMAP coordinates and the nearest neighbor data as returned by `UMAP` if `ret_nn = TRUE`.
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for the Leiden clustering.
#' @param metric Character. One of the metric used for computing the UMAP model. If NULL take the first available one.
#' @param partition_type 	Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details)
#' @param returnAsFactor If TRUE return the clusters attributions as a factor vector and not as characters.
#' @param n_iterations Number of iterations. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
#' @param resolution_parameter A parameter controlling the coarseness of the clusters.
#' @param seed Seed for the random number generator. By default uses a random seed if nothing is specified.
#' @param laplacian_init  Derive edge weights from the Laplacian matrix, otherwise weights are derived from the distance between cells. Not used for the moment
#' @param initial_membership Initial membership for the partition. If MULL then defaults to a singleton partition.
#' @param max_comm_size Maximal total size of nodes in a community. If zero (the default), then communities can be of any size.
#' @param ... Other parameters passed to `leidenFromPygraph`.
#'
#' @return A vector of character or factor if `returnAsFactor`, same length as number of samples in the UMAP. Cluster attribution of samples.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' umapWithNN<-UMAP(bulkLogCounts,ret_nn = TRUE)
#' proj2d(umapWithNN$embedding,colorBy = leidenFromUMAP(umapWithNN))
leidenFromUMAP <- function(umapWithNN,
													 n_neighbors = 10,
													 metric = NULL,
													 partition_type = c(
													 	"RBConfigurationVertexPartition",
													 	"ModularityVertexPartition",
													 	"RBERVertexPartition",
													 	"CPMVertexPartition",
													 	"MutableVertexPartition",
													 	"SignificanceVertexPartition",
													 	"SurpriseVertexPartition",
													 	"ModularityVertexPartition.Bipartite",
													 	"CPMVertexPartition.Bipartite"
													 ),
													 returnAsFactor = FALSE,
													 n_iterations = -1,
													 resolution_parameter = .5,
													 seed = 666,
													 laplacian_init = TRUE,
													 initial_membership = NULL,
													 max_comm_size = 0L,
													 ...) {
	if (is.null(metric)) {
		metric <- names(umapWithNN$nn)[1]
	} else{
		if (!metric %in% names(umapWithNN$nn))
			stop(metric, " distance metric is not computed in the provided model")
	}

	partition_type <- match.arg(partition_type)
	pygraph <-
		adjMat2Pygraph(getAdjMatfromUMAPWithNN(umapWithNN, n_neighbors = n_neighbors, metric =
																					 	metric))

	leidenFromPygraph(
		pygraph,
		returnAsFactor = returnAsFactor,
		n_iterations = n_iterations,
		resolution_parameter = resolution_parameter,
		seed = seed,
		laplacian_init = laplacian_init,
		initial_membership = initial_membership,
		max_comm_size = max_comm_size,
		partition_type = partition_type,
		...
	)
}


#' Python igraph object from an adjacency matrix.
#'
#' @param adjMat An adjacency matrix of a graph.
#' @param mode A character, the mode to be used. Possible values are:
#' -"directed" - the graph will be directed and a matrix element gives the number of edges between two vertex.
#' -"undirected" - alias to "max" for convenience.
#' -"max" - undirected graph will be created and the number of edges between vertex i and j is max(A(i,j),A(j,i))
#' -"min" - like "max", but with min(A(i,j),A(j,i))
#' -"plus" - like "max", but with A(i,j)+A(j,i)
#' -"upper" - undirected graph with the upper right triangle of the matrix (including the diagonal)
#' -"lower" - undirected graph with the lower left triangle of the matrix (including the diagonal)
#' @param ... Other parameters passed to the igraph python function `Weighted_Adjacency`.
#'
#' @return A python igraph object.
#' @export
#'
#' @examples
#' adjMat<-matrix(round(runif(25,min = 0,max = 1)),ncol = 5)
#' res<-adjMat2Pygraph(adjMat)
adjMat2Pygraph <- function(adjMat, mode = "directed", ...) {
	ig <- reticulate::import("igraph")
	graph <- ig$Graph$Weighted_Adjacency(matrix = adjMat, mode = mode, ...)
}

#' Compute a Leiden clustering from K Nearest Neighbors graph from python igraph.
#'
#' @param pygraph  A python igraph object as computed by `adjMat2Pygraph`.
#' @param returnAsFactor  If TRUE return the clusters attributions as a factor vector and not as characters.
#' @param n_iterations Number of iterations. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
#' @param resolution_parameter  A parameter controlling the coarseness of the clusters.
#' @param seed  Seed for the random number generator. By default uses a random seed if nothing is specified.
#' @param laplacian_init  Derive edge weights from the Laplacian matrix, otherwise weights are derived from the distance between cells. Not used for the moment.
#' @param initial_membership  Initial membership for the partition. If MULL then defaults to a singleton partition.
#' @param max_comm_size  Maximal total size of nodes in a community. If zero (the default), then communities can be of any size.
#' @param node_sizes Vector of numerics. The quality function takes into account the size of a community, which is defined as the sum over the sizes of each individual node. By default, the node sizes are set to 1, meaning that the size of a community equals the number of nodes of a community. If a node already represents an aggregation, this could be reflect in its node size.
#' @param weight_parameter Weight of the graph as a numeric value for each edge.  Not used for the moment.
#' @param partition_type 	Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details).
#'
#' @return  A vector of character or factor if `returnAsFactor`, same length as number of nodes in the graph. Cluster attribution of samples.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' umapWithNN<-UMAP(bulkLogCounts,ret_nn = TRUE)
#' pygraph<-adjMat2Pygraph(getAdjMatfromUMAPWithNN(umapWithNN))
#' proj2d(umapWithNN$embedding,colorBy = leidenFromPygraph(pygraph))
leidenFromPygraph <-
	function(pygraph,
					 returnAsFactor = FALSE,
					 n_iterations = -1,
					 resolution_parameter = .5,
					 seed = 666,
					 laplacian_init = TRUE,
					 initial_membership = NULL,
					 max_comm_size = 0L,
					 node_sizes = NULL,
					 weight_parameter = NULL,
					 partition_type = c(
					 	"RBConfigurationVertexPartition",
					 	"ModularityVertexPartition",
					 	"RBERVertexPartition",
					 	"CPMVertexPartition",
					 	"MutableVertexPartition",
					 	"SignificanceVertexPartition",
					 	"SurpriseVertexPartition",
					 	"ModularityVertexPartition.Bipartite",
					 	"CPMVertexPartition.Bipartite"
					 )) {
		ig <- import("igraph")
		numpy <- import("numpy", delay_load = TRUE)
		leidenalg <- import("leidenalg", delay_load = TRUE)
		#
		# if (laplacian_init && is.null(weight_parameter)) {
		# 	pygraph <- ig$Graph$simplify(pygraph, multiple = TRUE, loops = TRUE)
		# 	laplacian <- do.call("cbind",ig$Graph$laplacian(pygraph))
		# 	pygraph$es$set_attribute_values("weight",value = -as.matrix(laplacian)[as.matrix(laplacian) < 	0])
		# }
		graphIsWeighted <- ig$Graph$is_weighted(pygraph)

		w <- unlist(pygraph$es$get_attribute_values("weight"))

		partition_type <- match.arg(partition_type)
		if (partition_type == "ModularityVertexPartition.Bipartite")
			degree_as_node_size <- TRUE
		if (!is.null(seed))
			seed <- as.integer(seed)
		if (!is.integer(n_iterations))
			n_iterations <- as.integer(n_iterations)
		max_comm_size <- as.integer(max_comm_size)

		part <- switch(
			EXPR = partition_type,
			RBConfigurationVertexPartition = leidenalg$find_partition(
				pygraph,
				leidenalg$RBConfigurationVertexPartition,
				initial_membership = initial_membership,
				weights = w,
				seed = seed,
				n_iterations = n_iterations,
				max_comm_size = max_comm_size,
				resolution_parameter = resolution_parameter
			),

			ModularityVertexPartition = leidenalg$find_partition(
				pygraph,
				leidenalg$ModularityVertexPartition,
				initial_membership = initial_membership,
				weights = w,
				seed = seed,
				n_iterations = n_iterations,
				max_comm_size = max_comm_size
			),

			RBERVertexPartition = leidenalg$find_partition(
				pygraph,
				leidenalg$RBERVertexPartition,
				initial_membership = initial_membership,
				weights = w,
				seed = seed,
				n_iterations = n_iterations,
				max_comm_size = max_comm_size,
				node_sizes = node_sizes,
				resolution_parameter = resolution_parameter
			),

			CPMVertexPartition = leidenalg$find_partition(
				pygraph,
				leidenalg$CPMVertexPartition,
				initial_membership = initial_membership,
				weights = w,
				seed = seed,
				n_iterations = n_iterations,
				max_comm_size = max_comm_size,
				node_sizes = node_sizes,
				resolution_parameter = resolution_parameter
			),

			MutableVertexPartition = leidenalg$find_partition(
				pygraph,
				leidenalg$MutableVertexPartition,
				initial_membership = initial_membership,
				seed = seed,
				n_iterations = n_iterations,
				max_comm_size = max_comm_size
			),

			SignificanceVertexPartition = leidenalg$find_partition(
				pygraph,
				leidenalg$SignificanceVertexPartition,
				initial_membership = initial_membership,
				seed = seed,
				n_iterations = n_iterations,
				max_comm_size = max_comm_size,
				node_sizes = node_sizes,
				resolution_parameter = resolution_parameter
			),

			SurpriseVertexPartition = leidenalg$find_partition(
				pygraph,
				leidenalg$SurpriseVertexPartition,
				initial_membership = initial_membership,
				weights = w,
				seed = seed,
				n_iterations = n_iterations,
				max_comm_size = max_comm_size,
				node_sizes = node_sizes
			),

			ModularityVertexPartition.Bipartite = run_bipartite_partitioning(
				pygraph,
				initial_membership = initial_membership,
				weights = w,
				resolution_parameter_01 = resolution_parameter,
				resolution_parameter_0 = 0,
				resolution_parameter_1 = 0,
				degree_as_node_size = TRUE,
				types = "type",
				seed = seed,
				n_iterations = n_iterations
			),

			CPMVertexPartition.Bipartite = run_bipartite_partitioning(
				pygraph,
				initial_membership = initial_membership,
				weights = w,
				resolution_parameter_01 = resolution_parameter,
				resolution_parameter_0 = 0,
				resolution_parameter_1 = 0,
				degree_as_node_size = degree_as_node_size,
				types = "type",
				seed = seed,
				n_iterations = n_iterations
			),
			stop(
				"please specify a partition type as a string out of those documented"
			)
		)

		res <- paste0("k", formatNumber2Character(part$membership + 1))
		if (returnAsFactor)
			res <- as.factor(res)
		res
	}


#' Compute an adjacency matrix from the nearest neighbor matrix of a UMAP model
#'
#' @param umapWithNN A UMAP model as returned by UMAP if `ret_nn = TRUE`,
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for the Leiden clustering.
#' @param metric  Character. One of the metric used for computing the UMAP model. If NULL take the first available one.
#'
#' @return A sparse adjacency matrix from the class `dgTMatrix`.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' umapWithNN<-UMAP(bulkLogCounts,ret_nn = TRUE)
#' getAdjMatfromUMAPWithNN(umapWithNN)
getAdjMatfromUMAPWithNN <-
	function(umapWithNN,
					 n_neighbors = 10,
					 metric = NULL) {
		if (is.null(metric)) {
			metric <- names(umapWithNN$nn)[1]
		} else{
			if (!metric %in% names(umapWithNN$nn))
				stop(metric, " distance metric is not computed in the provided model")
		}
		if (ncol(umapWithNN$nn[[metric]]$idx) < n_neighbors) {
			stop(
				"The provided umap model contains the data for ",
				ncol(umapWithNN$nn[[metric]]$idx),
				" neighbors, please decrease the n_neighbors parameter or recompute the model on a higher number of neighbors"
			)
		}
		#from 2 --> rm diagonals
		knn_indices <- umapWithNN$nn[[metric]]$idx[, 2:n_neighbors]
		knn_dists <- umapWithNN$nn[[metric]]$dist[, 2:n_neighbors]

		n <- nrow(knn_indices)
		as(
			Matrix::sparseMatrix(
				i = as.integer(rep(1:n, each = ncol(knn_indices))),
				j = as.integer(as.vector(t(knn_indices))),
				x = as.numeric(as.vector(t(knn_dists))),
				dims = c(n, n)
			),
			"dgTMatrix"
		)
	}


#' Execute a fastMNN and rescale the counts.
#'
#' @description
#' Use fastMNN fron the package batchelor, but return a matrix. Rescale the counts so the range of each gene remains the same after the transformation.
#' @param logCounts A matrix of numeric (in the RNA-Seq context, log counts).
#' @param batch A vector or factor specifying the batch of origin for all cells
#' @param k An integer scalar specifying the number of nearest neighbors to consider when identifying MNNs.
#' @param returnRescale Logical. Use reScale on the output so the dymaic range after bacth correction of genes is the same than before.
#' @param ... Other parameters passed to fastMNN.
#'
#' @return A matrix of corrected count table.
#' @export
#'
#' @examples
#' data("geneLengthGRCh38")
#' library(MASS)
#' countMat<-sapply(vector("numeric",length = 100),function(x){
#'  	c(rnegbin(10,mu = 50,theta = 5), rnegbin(10,mu = 10,theta = 5))
#' }) |> t(); countMat<-log2(countMat+1)
#' heatmap.DM(countMat)
#' correctedMat<-oobFastMNN(countMat,batch = c(rep(1,10),rep(2,10)), k=5) #warning because of the small matrix
#' heatmap.DM(correctedMat)
oobFastMNN <- function(logCounts, batch, k, returnRescale = TRUE, ...) {
	scObj <- batchelor::fastMNN(logCounts, batch = batch, k = k, ...)
	res <- SummarizedExperiment::assay(scObj, "reconstructed") |> as.matrix()
	if (returnRescale)
		res <- reScale(res, logCounts)
	res
}
