
#' Title
#'
#' @param d
#' @param transpose
#' @param scale
#' @param center
#'
#' @return
#' @export
#'
#' @examples
PCA<-function(d,transpose=T,scale=F,center=T) {
	if(transpose) d<-t(d);
	means<-0;sdeviations<-1
	if(center){
		means<-apply(d,2,mean)
		d<-sweep(d,2,means,"-")
	}
	if(scale){
		sdeviations<-apply(d,2,sd)
		d<-sweep(d,2,sdeviations,"/")
	}
	resacp <-prcomp(x = d,retx = T,center = FALSE,scale = FALSE);
	resacp$n.obs<-dim(d)[1];
	resacp$percentVar<- resacp$sdev^2 / sum( resacp$sdev^2 )
	resacp$scale<-scale
	resacp$center<-center
	resacp$transform<-list(sdeviations=sdeviations,means=means)
	return(resacp);
}


#' Title
#'
#' @param d
#' @param transpose
#' @param scale
#' @param center
#' @param nPC
#' @param weight.by.var
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fastPCA <- function(d,transpose=TRUE,scale=FALSE,center=TRUE,nPC=NULL,
										weight.by.var = TRUE, ...) {
	require(irlba)
	if(transpose) d<-t(d);
	d<-as.matrix(d)
	means<-0;sdeviations<-1
	if(center | scale) d<-scale(d,scale = scale,center = center)
	if(center) means<-attr(d,"scaled:center")
	if(scale) sdeviations<-attr(d,"scaled:scale")

	if(is.null(nPC | nPC>=min(nrow(d),ncol(d)))) {
		nPC <- min(ncol(d) - 1, nrow(d) -1)
	}

	resacp<-list()
	resacp$n.obs<-dim(d)[1];
	resacp$scale<-scale
	resacp$center<-center
	resacp$transform<-list(sdeviations=sdeviations,means=means)

	irlbaResults <- irlba(A = d, nv = nPC, ...)
	rotation <- irlbaResults$v
	resacp$sdev <- irlbaResults$d/sqrt(max(1, nrow(d) - 1))
	if (weight.by.var) {
		reducedSpace <- irlbaResults$u %*% diag(irlbaResults$d)
	} else {
		reducedSpace <- irlbaResults$u
	}
	rownames(rotation) <- colnames(d)
	colnames(rotation) <- paste0("PC", 1:nPC)
	rownames(reducedSpace) <- rownames(d)
	colnames(reducedSpace) <- colnames(rotation)
	resacp$x<-reducedSpace
	resacp$rotation<-rotation
	resacp$percentVar<- resacp$sdev^2 / sum( resacp$sdev^2 )
	resacp
}

#' Title
#'
#' @param pca
#' @param newSamplesMatrix
#' @param transpose
#' @param combineMat
#'
#' @return
#' @export
#'
#' @examples
pcaAddSamples<-function(pca,newSamplesMatrix,transpose=TRUE,combineMat=TRUE){
	if(transpose) newSamplesMatrix<-t(newSamplesMatrix)
	newSamplesMatrix<-newSamplesMatrix[,rownames(pca$rotation),drop=FALSE]

	if(pca$center) newSamplesMatrix<-sweep(newSamplesMatrix,2,pca$transform$means,"-")
	if(pca$scale) newSamplesMatrix<-sweep(newSamplesMatrix,2,pca$transform$sdeviations,"/")

	newSamplesCoord<-t(apply(newSamplesMatrix,1,function(x){
		colSums(x*pca$rotation)
	}))
	if(combineMat){
		pca$x<-rbind(pca$x,newSamplesCoord)
		return(pca)
	}else{
		return(newSamplesCoord)
	}
}


#' Title
#'
#' @param pca
#' @param annotationDF
#' @param nComponent
#'
#' @return
#' @export
#'
#' @examples
PCR<-function(pca,annotationDF,nComponent=10){
	require(reshape2)
	annots<-colnames(annotationDF)
	rSquaredMat<-matrixFromDimnames(cn(pca$x),annots,value = NA)
	for(x in cn(pca$x)){
		for(annot in annots){
			rSquaredMat[x,annot]<-summary(lm(formula(paste0(x,"~",annot)),data = data.frame(annotationDF[,annot,drop=F],pca$x[,x,drop=F])))$r.squared
		}
	}
	if(ncol(annotationDF)<2){
		retDt<-data.frame(Rsquared=rSquaredMat[,1],PC=rownames(rSquaredMat),Annotation=colnames(annotationDF))
	}else{
		retDt<-melt(rSquaredMat[1:nComponent,], value.name = "Rsquared",varnames=c("PC","Annotation"))
	}

	retDt$PC<-factor(retDt$PC,levels=cn(pca$x)[1:nComponent]) #so the levels of PCs are well ordered
	retDt
}


#' Title
#'
#' @param data
#' @param transpose
#' @param scale
#' @param center
#' @param metric
#' @param ndim
#' @param maxit
#'
#' @return
#' @export
#'
#' @examples
NMDS<-function(data,transpose=TRUE,scale=FALSE,center=FALSE,metric=dist,ndim=2,maxit=100){
	merged<-FALSE
	require(MASS)
	if(transpose) data <- t(data)
	d <- metric(data)  # euclidean distances between the rows
	if(min(d,na.rm=TRUE)==0){
		merged<-TRUE
		md<-merge0dist(d)
		d<-md$distMat
		mergedSample<-md$merged
	}
	fit <- isoMDS(d, k=ndim, maxit=maxit) # k is the number of dim
	fit$coord<-fit$points
	fit$points<-NULL
	if(merged){
		for(sple in names(mergedSample)){
			values<-matrix(rep(fit$coord[sple,],length(mergedSample[[sple]])),nrow=length(mergedSample[[sple]]),byrow = TRUE)
			rownames(values)<-mergedSample[[sple]]
			fit$coord<-rbind(fit$coord,values)
		}
	}
	return(fit)
}

#' Title
#'
#' @param data
#' @param nDimPCA
#' @param transpose
#' @param n_neighbors
#' @param n_components
#' @param min_dist
#' @param init
#' @param metric
#' @param ret_model
#' @param ret_nn
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
make.umap<-function(data,nDimPCA=NULL,transpose=TRUE,n_neighbors=NULL, n_components = 2,min_dist=0.01,
										init = "laplacian", metric = "euclidean",ret_model=FALSE,ret_nn=FALSE,...){
	require(uwot)
	if(transpose) data<-t(data)
	if(is.null(n_neighbors)) n_neighbors=nrow(data)
	if(!is.null(nDimPCA)){
		data<-fastPCA(data,transpose = FALSE,scale = FALSE,nPC = nDimPCA)$x
	}
	res<-uwot::umap(as.matrix(data),n_neighbors = n_neighbors, n_components = n_components,
									min_dist=min_dist, init = init, metric = metric,ret_model=ret_model,ret_nn = ret_nn,...)
	if(!ret_model & !ret_nn) rownames(res)<-rownames(data)
	res
}


#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
CPM<-function(data){ #Normalisation CPM
	data.CPM <- sweep(data, 2, colSums(data),`/`)
	data.CPM <-data.CPM * 1000000
	return(data.CPM)
}


#' Title
#'
#' @param data
#' @param gene.length
#'
#' @return
#' @export
#'
#' @examples
TPMfullLength<-function(data, gene.length){
	gene.length.kb <- gene.length[rn(data)]/1000
	data<-sweep(data, 1, gene.length.kb,`/`)
	return(CPM(data))
}

#' Title
#'
#' @param data
#' @param gene.length
#'
#' @return
#' @export
#'
#' @examples
RPKM<-function(data, gene.length){
	gene.length.kb <- gene.length[rn(data)]/1000
	data<-CPM(data)
	sweep(data, 1, gene.length.kb,`/`)
}


#' Title
#'
#' @param hc
#' @param min
#' @param max
#' @param loss
#' @param graph
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
best.cutree <- function(hc, min=2, max=20, loss=FALSE, graph=FALSE, ...){
	if (class(hc)!="hclust") hc <- as.hclust(hc)
	max <- min(max, length(hc$height)-1)
	inert.gain <- rev(hc$height)
	intra <- rev(cumsum(rev(inert.gain)))
	relative.loss = intra[min:(max+1)]/intra[(min - 1):(max)]
	derivative.loss = relative.loss[2:length(relative.loss)]-relative.loss[1:(length(relative.loss)-1)]
	names(derivative.loss) <- min:max
	if (graph) {
		print(
			ggplot(data.frame(partition=min:max,derivative.loss=derivative.loss),aes(x=partition,y=derivative.loss))+
				geom_point()+
				scale_x_continuous(breaks=min:max)
		)
	} else {
		if (loss)
			derivative.loss
		else
			as.numeric(names(which.max(derivative.loss)))
	}
}


# x : matrix
#' Title
#'
#' @param x
#' @param transpose
#' @param method.dist
#' @param method.hclust
#' @param bootstrap
#' @param nboot
#' @param PCAfirst
#' @param nDimPCA
#'
#' @return
#' @export
#'
#' @examples
hierarchicalClustering<-function(x,transpose=TRUE,method.dist="euclidean",method.hclust="ward.D2",
																 bootstrap=FALSE,nboot=10,PCAfirst=FALSE,nDimPCA=NULL){
	if(transpose) x<-t(x)
	if(PCAfirst){
		x<-ACP(x,transpose = FALSE,scale = FALSE)$x
		if(!is.null(nDimPCA)){
			x<-x[,1:nDimPCA]
		}
	}
	if(bootstrap){
		require(pvclust)
		resClust<-pvclust(t(x),nboot=nboot,method.hclust = method.hclust,parallel = TRUE,method.dist = method.dist)$hclust
	}else{
		if(method.dist=="pearson"){
			resDist<-corrDist(x)
		}else if(method.dist=="bicor"){
			require("WGCNA")
			resDist<-as.dist((1 - suppressWarnings(bicor(Matrix::t(x))))/2)
		}else{
			resDist<-dist(x, method = method.dist)
		}
		resClust<-stats::hclust(resDist,method = method.hclust)
	}
	return(resClust)
}


#' Title
#'
#' @param countMatrix
#'
#' @return
#' @export
#'
#' @examples
normDeseq<-function(countMatrix){ #matrix where genes are rows and samples are columns
	# PS = pseudo reference sample
	PS<-apply(countMatrix,1,gmean,keepZero=TRUE) #get a vector which consist of the geometrical mean of each genes across all samples
	keptRow<-PS>0 #get rid of genes containing one zero ore more
	PS<-PS[keptRow]
	ratioMat<-sweep(countMatrix[keptRow,],1,PS,"/") #get the ratio matrix (expression/expression from PS)
	normFactors<-apply(ratioMat,2,median) #get the median of the ratios for each sample to get the normalization factors
	sweep(countMatrix,2,normFactors,"/") #divide each sample by the corresponding normalization factor
}

#' Title
#'
#' @param exprData
#' @param sampleData
#' @param contrast
#'
#' @return
#' @export
#'
#' @examples
multiLinearModel<-function(exprData,sampleData,contrast){
	samples<-rownames(sampleData)[sampleData[,contrast[1]]%in%contrast[2:3]]
	data<-exprData[,samples]
	groups<-droplevels(as.factor(sampleData[samples,contrast[1]]))
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


#' Title
#'
#' @param counts
#' @param minCount
#' @param plot
#' @param returnPlot
#'
#' @return
#' @export
#'
#' @examples
getMostVariableGenes<-function(counts,minCount=0.01,plot=TRUE,returnPlot=FALSE){
	counts<-counts[rowMeans(counts)>minCount,]
	dispTable<-data.frame(mu=rowMeans(counts),var=apply(counts,1,var),row.names =rownames(counts))
	dispTable$cv2<- dispTable$var / dispTable$mu^2
	sumNullvariance<-sum(dispTable$cv2 <= 0)
	if(sumNullvariance>0){
		warning(paste0(sumNullvariance, " have null variance and will be removed"))
		dispTable<-dispTable[dispTable$cv2 > 0,]
	}
	fit<-loess(cv2~mu,data = log10(dispTable[,c("mu","cv2")]))
	dispTable$residuals<-fit$residuals
	dispTable$residuals2<-dispTable$residuals^2
	dispTable$fitted<-10^fit$fitted
	if(plot){
		require(ggplot2)
		require(ggrepel)
		require(circlize)
		g<-ggplot(dispTable,aes(x=mu,y=cv2,label=rownames(dispTable),fill=residuals))+
			geom_point(stroke=1/8,colour = "black",shape=21)+geom_line(aes(y=fitted),color="red",size=1.5)+
			scale_x_log10()+scale_y_log10()
		if(returnPlot){
			return(g)
		}else{
			print(g)
		}
	}
	dispTable
}

#' Title
#'
#' @param logCounts
#' @param minCount
#' @param plot
#' @param returnPlot
#'
#' @return
#' @export
#'
#' @examples
getMostVariableGenesLogCount<-function(logCounts,minCount=0,plot=TRUE,returnPlot=FALSE){
	if(minCount>0) logCounts<-logCounts[rowMeans(logCounts)>minCount,]
	dispTable<-data.frame(mu=rowMeans(logCounts),var=apply(logCounts,1,var),row.names =rownames(logCounts))

	fit<-loess(var~mu,data = dispTable)
	if(length(rownames(fit$x)) < nrow(dispTable)){
		warning(nrow(dispTable)-length(rownames(fit$x))," genes have null variance and will be removed")
	}
	dispTable<-dispTable[rownames(fit$x),]
	dispTable$residuals<-fit$residuals
	dispTable$fitted<-fit$fitted

	if(plot){
		require(ggplot2)
		require(ggrepel)
		require(circlize)
		g<-ggplot(dispTable,aes(x=mu,y=var,label=rownames(dispTable),fill=residuals))+
			geom_point(stroke=1/8,colour = "black",shape=21)+geom_line(aes(y=fitted),color="red",size=1.5)
		if(returnPlot){
			return(g)
		}else{
			print(g)
		}
	}
	dispTable
}


# By Miron Kursa https://mbq.me
#' Title
#'
#' @param score
#' @param bool
#'
#' @return
#' @export
#'
#' @examples
auroc <- function(score, bool) {
	n1 <- sum(!bool)
	n2 <- sum(bool)
	U	<- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
	return(1 - U / n1 / n2)
}


#' Title
#'
#' @param expressionMatrix
#' @param group
#' @param BPPARAM
#'
#' @return
#' @export
#'
#' @examples
getMarkers<-function(expressionMatrix,group,BPPARAM=NULL){
	if(is.null(BPPARAM )) BPPARAM=bpparam()
	if(!is.matrix(expressionMatrix)) expressionMatrix<-as.matrix(expressionMatrix)
	if(length(group)!=ncol(expressionMatrix)) stop("group should be a vector with same length as number of column in expressionMatrix")

	group<-as.factor(group)
	group<-droplevels(group)
	binaryGroup<-lapply(levels(group),function(x){
		x==group & !is.na(group)
	});names(binaryGroup)<-levels(group)

	res<-as.matrix(data.frame(lapply(binaryGroup,function(labels){
		unlist(bplapply(seq_len(nrow(expressionMatrix)),function(i,expressionMatrix,labels,auroc){
			auroc(expressionMatrix[i,],labels)
		},labels=labels,expressionMatrix=expressionMatrix,auroc=auroc,BPPARAM=BPPARAM),recursive = FALSE)
	})))
	rownames(res)<-rownames(expressionMatrix)
	res
}


#' Title
#'
#' @param gene
#' @param expression
#' @param corFun
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
corGeneToOthers<-function(gene,expression,corFun=cor,...){
	expression<-as.matrix(expression)
	t(corFun(expression[gene,],t(expression),...))[,1]
}


#' Title
#'
#' @param exprMatrix
#' @param genes
#' @param scale
#' @param center
#' @param returnContribution
#'
#' @return
#' @export
#'
#' @examples
eigengenes<-function(exprMatrix,genes,scale=F,center=T,returnContribution=F){
	pca<-prcomp(x = t(exprMatrix[genes,]),retx = T,center = center,scale = scale)
	eigen<-pca$x[,1]
	contribution<-pca$rotation[,1]
	if(cor(colMeans(exprMatrix[genes,]),eigen)<0){
		eigen<- -eigen
		contribution<- -contribution
	}
	if(returnContribution){
		list(eigengenes=eigen,contribution=contribution)
	}else{
		eigen
	}
}


#' Title
#'
#' @param rawCounts
#' @param returnLog
#' @param sizeFactors
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
quickSCnorm<-function(rawCounts,returnLog=TRUE,sizeFactors=NULL,...){
	require(scran)
	sce <- SingleCellExperiment(assays=list(counts=rawCounts))
	if(!is.null(sizeFactors)){
		sizeFactors(sce)<-sizeFactors
	}else{
		sce <- computeSumFactors(sce,...)
	}
	scater::normalizeCounts(sce,log=returnLog,pseudo_count = 0,size_factors = sizeFactors(sce))
}


#' @title pacmap
#'
#' @description An R wrapper for the PaCMAP Python module found at
#' https://github.com/YingfanWang/PaCMAP
#'
#' PaCMAP (Pairwise Controlled Manifold Approximation) is a
#' dimensionality reduction method that can be used for
#' visualization, preserving both local and global structure
#' of the data in original space. PaCMAP optimizes the low
#' dimensional embedding using three kinds of pairs of points:
#' neighbor pairs (pair_neighbors), mid-near pair (pair_MN),
#' and further pairs (pair_FP).
#'
#' @param rdf A variable by observation data frame
#' @param n_components integer Dimensions of the embedded space. Default: 3
#' @param perplexity numeric The perplexity is related to the
#' number of nearest neighbors that is used in other manifold learning
#' algorithms. Larger datasets usually require a larger perplexity. Consider
#' selecting a value between 5 and 50. The choice is not extremely critical
#' since t-SNE is quite insensitive to this parameter. Default: 30
#' @param early_exaggeration numeric Controls how tight natural
#' clusters in the original space are in the embedded space and how much space
#' will be between them. For larger values, the space between natural clusters
#' will be larger in the embedded space. Again, the choice of this parameter
#' is not very critical. If the cost function increases during initial
#' optimization, the early exaggeration factor or the learning rate might be
#' too high. Default: 12.0
#' @param learning_rate numeric The learning rate for t-SNE is
#' usually in the range [10.0, 1000.0]. If the learning rate is too high, the
#' data may look like a ‘ball’ with any point approximately equidistant from
#' its nearest neighbours. If the learning rate is too low, most points may
#' look compressed in a dense cloud with few outliers. If the cost function
#' gets stuck in a bad local minimum increasing the learning rate may help. Default: 200.0
#' @param n_iter integer Maximum number of iterations for the
#' optimization. Should be at least 250. Default: 1000
#' @param n_iter_without_progress integer Maximum number of
#' iterations without progress before we abort the optimization, used after
#' 250 initial iterations with early exaggeration. Note that progress is only
#' checked every 50 iterations so this value is rounded to the next multiple
#' of 50. Default: 300
#' @param min_grad_norm numeric If the gradient norm is below
#' this threshold, the optimization will be stopped. Default: 1e-7
#' @param metric character or callable The metric to use when calculating distance
#' between instances in a feature array. If metric is a character, it must be one
#' of the options allowed by scipy.spatial.distance.pdist for its metric
#' parameter, or a metric listed in pairwise.PAIRWISE.DISTANCE.FUNCTIONS. If
#' metric is “precomputed”, X is assumed to be a distance matrix.
#' Alternatively, if metric is a callable function, it is called on each pair
#' of instances (rows) and the resulting value recorded. The callable should
#' take two arrays from X as input and return a value indicating the distance
#' between them. The default is “euclidean” which is interpreted as squared
#' euclidean distance.
#' @param init character or numpy array Initialization of
#' embedding. Possible options are ‘random’, ‘pca’, and a numpy array of shape
#' (n.samples, n.components). PCA initialization cannot be used with
#' precomputed distances and is usually more globally stable than random
#' initialization. Default: “random”
#' @param verbose integer Verbosity level. Default: 0
#' @param random_state int, RandomState instance or NULL If int,
#' random.state is the seed used by the random number generator; If
#' RandomState instance, random.state is the random number generator; If NULL,
#' the random number generator is the RandomState instance used by np.random.
#' Note that different initializations might result in different local minima
#' of the cost function. Default: NULL
#' @param method character By default the gradient
#' calculation algorithm uses Barnes-Hut approximation running in \eqn{O(N log N)}
#' time. method=’exact’ will run on the slower, but exact, algorithm in \eqn{O(N^2)}
#' time. The exact algorithm should be used when nearest-neighbor errors need
#' to be better than 3%. However, the exact method cannot scale to millions of
#' examples. Default: ‘barnes.hut’
#' @param angle numeric Only used if method=’barnes.hut’ This is
#' the trade-off between speed and accuracy for Barnes-Hut T-SNE. ‘angle’ is
#' the angular size (also referred to as theta) of a distant node as
#' measured from a point. If this size is below ‘angle’ then it is used as a
#' summary node of all points contained within it. This method is not very
#' sensitive to changes in this parameter in the range of 0.2 - 0.8. Angle
#' less than 0.2 has quickly increasing computation time and angle greater 0.8
#' has quickly increasing error.#' Default: 0.5
#' @param auto_iter boolean Should optimal parameters be determined?
#' If false, behaves like stock MulticoreTSNE Default: TRUE
#' @param auto_iter_end intNumber of iterations for parameter
#' optimization. Default: 5000
#' @param n_jobs Number of processors to use.  Default: all.
#'
#' @importFrom reticulate import py_module_available
#' @importFrom parallel detectCores
#'
#' @return data.frame with tSNE coordinates
#' @export
#'
pacmap <- function(rdf,
									 n_dims         = 2,
									 n_neighbors    = NULL,
									 MN_ratio       = 0.5,
									 FP_ratio       = 2.0,
									 pair_neighbors = NULL,
									 pair_MN        = NULL,
									 pair_FP        = NULL,
									 distance       = "euclidean",
									 lr             = 1.0,
									 num_iters      = 450,
									 verbose        = FALSE,
									 apply_pca      = TRUE,
									 intermediate   = FALSE,
									 init=NULL){

	require(reticulate)
	if (!py_module_available("pacmap")){
		stop("The pacmap module is unavailable.  Please activate the appropriate environment or install the module.")
	}

	pacmap_module <-import(
		module = "pacmap",
		delay_load = TRUE
	)

	pacmap <- pacmap_module$PaCMAP(
		n_dims         = as.integer(n_dims),
		n_neighbors    = as.integer(n_neighbors),
		MN_ratio       = as.numeric(MN_ratio),
		FP_ratio       = as.numeric(FP_ratio),
		pair_neighbors = pair_neighbors,
		pair_MN        = pair_MN,
		pair_FP        = pair_FP,
		distance       = distance,
		lr             = as.numeric(lr),
		num_iters      = as.integer(num_iters),
		verbose        = verbose,
		apply_pca      = apply_pca,
		intermediate   = intermediate
	)

	pacmap$fit_transform(rdf,init=init)
}


#' Title
#'
#' @param rdf
#' @param n_dims
#' @param n_inliers
#' @param n_outliers
#' @param apply_pca
#' @param n_iters
#' @param knn_tuple
#'
#' @return
#' @export
#'
#' @examples
trimap<-function(rdf,n_dims = 2,n_inliers = 10,n_outliers = 5,
								 apply_pca=TRUE,n_iters=400,knn_tuple=NULL){
	# n_dims: Number of dimensions of the embedding (default = 2)
	# n_inliers: Number of nearest neighbors for forming the nearest neighbor triplets (default = 10).
	# n_outliers: Number of outliers for forming the nearest neighbor triplets (default = 5).
	# n_random: Number of random triplets per point (default = 5).
	# distance: Distance measure ('euclidean' (default), 'manhattan', 'angular', 'hamming')
	# weight_adj: The value of gamma for the log-transformation (default = 500.0).
	# lr: Learning rate (default = 1000.0).
	# n_iters: Number of iterations (default = 400).
	#
	# The other parameters include:
	#
	# 	knn_tuple: Use the precomputed nearest-neighbors information in form of a tuple (knn_nbrs, knn_distances) (default = None)
	# use_dist_matrix: Use the precomputed pairwise distance matrix (default = False)
	# apply_pca: Reduce the number of dimensions of the data to 100 if necessary before applying the nearest-neighbor search (default = True).
	# opt_method: Optimization method {'sd' (steepest descent), 'momentum' (GD with momentum), 'dbd' (delta-bar-delta, default)}.
	# verbose: Print the progress report (default = True).
	# return_seq: Store the intermediate results and return the results in a tensor (default = False).
	#

	trimap_module <- reticulate::import( module = "trimap", delay_load = TRUE)

	trimap <- trimap_module$TRIMAP(
		n_dims = as.integer(n_dims),
		n_inliers = as.integer(n_inliers),
		n_outliers = as.integer(n_outliers),
		apply_pca=apply_pca,
		n_iters = as.integer(n_iters),
		knn_tuple=knn_tuple
	)

	trimap$fit_transform(rdf)

}

#' Title
#'
#' @param umapModel
#' @param n_neighbors
#' @param metric
#' @param partition_type
#' @param returnAsFactor
#' @param n_iterations
#' @param resolution_parameter
#' @param seed
#' @param laplacian_init
#' @param initial_membership
#' @param max_comm_size
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
leidenFromUMAP<-function(umapModel,n_neighbors=10,metric=NULL,
												 partition_type=c("RBConfigurationVertexPartition", "ModularityVertexPartition","RBERVertexPartition", "CPMVertexPartition", "MutableVertexPartition",
												 								 "SignificanceVertexPartition", "SurpriseVertexPartition","ModularityVertexPartition.Bipartite", "CPMVertexPartition.Bipartite"),
												 returnAsFactor=FALSE,n_iterations=-1,resolution_parameter = .5,seed=666,laplacian_init = TRUE,
												 initial_membership=NULL, max_comm_size = 0L,...){

	require(leiden)
	require(igraph)

	if(is.null(metric)){
		metric<-names(umapModel$nn)[1]
	}else{
		if(!metric %in% names(umapModel$nn)) stop(metric, " distance metric is not computed in the provided model")
	}

	partition_type <- match.arg(partition_type)
	pygraph<-adjMat2Pygraph(getAdjMatfromUMAPmodel(umapModel,n_neighbors=n_neighbors,metric=metric))

	leidenFromPygraph(pygraph,returnAsFactor=returnAsFactor,n_iterations=n_iterations,
										resolution_parameter = resolution_parameter,seed=seed,laplacian_init = laplacian_init,
										initial_membership=initial_membership, max_comm_size = max_comm_size,partition_type = partition_type,...)
}


adjMat2Pygraph<-function(adjMat,mode = "directed",...){
	ig<-reticulate::import("igraph")
	graph<-ig$Graph$Weighted_Adjacency(matrix = adjMat,mode = mode,...)
}

leidenFromPygraph<-function(pygraph,returnAsFactor=FALSE,n_iterations=-1,resolution_parameter = .5,seed=666,laplacian_init = TRUE,
														initial_membership=NULL, max_comm_size = 0L,node_sizes = NULL,weight_parameter=NULL,
														partition_type = c("RBConfigurationVertexPartition", "ModularityVertexPartition",
																							 "RBERVertexPartition", "CPMVertexPartition", "MutableVertexPartition",
																							 "SignificanceVertexPartition", "SurpriseVertexPartition",
																							 "ModularityVertexPartition.Bipartite", "CPMVertexPartition.Bipartite"),...){
	require(reticulate)
	ig<-import("igraph")
	numpy <- import("numpy", delay_load = TRUE)
	leidenalg <- import("leidenalg", delay_load = TRUE)
	#
	# if (laplacian_init && is.null(weight_parameter)) {
	# 	pygraph <- ig$Graph$simplify(pygraph, multiple = TRUE, loops = TRUE)
	# 	laplacian <- do.call("cbind",ig$Graph$laplacian(pygraph))
	# 	pygraph$es$set_attribute_values("weight",value = -as.matrix(laplacian)[as.matrix(laplacian) < 	0])
	# }
	graphIsWeighted<-ig$Graph$is_weighted(pygraph)

	w<-unlist(pygraph$es$get_attribute_values("weight"))

	partition_type <- match.arg(partition_type)
	if (partition_type == "ModularityVertexPartition.Bipartite") degree_as_node_size <- TRUE
	if (!is.null(seed))  seed <- as.integer(seed)
	if (!is.integer(n_iterations)) n_iterations <- as.integer(n_iterations)
	max_comm_size<-as.integer(max_comm_size)

	part <- switch(EXPR = partition_type,

								 RBConfigurationVertexPartition = leidenalg$find_partition(pygraph,
								 																													leidenalg$RBConfigurationVertexPartition, initial_membership = initial_membership,
								 																													weights = w, seed = seed, n_iterations = n_iterations,
								 																													max_comm_size = max_comm_size,
								 																													resolution_parameter = resolution_parameter),

								 ModularityVertexPartition = leidenalg$find_partition(pygraph,
								 																										 leidenalg$ModularityVertexPartition, initial_membership = initial_membership,
								 																										 weights = w, seed = seed, n_iterations = n_iterations,
								 																										 max_comm_size = max_comm_size),

								 RBERVertexPartition = leidenalg$find_partition(pygraph,
								 																							 leidenalg$RBERVertexPartition, initial_membership = initial_membership,
								 																							 weights = w, seed = seed, n_iterations = n_iterations,
								 																							 max_comm_size = max_comm_size,
								 																							 node_sizes = node_sizes, resolution_parameter = resolution_parameter),

								 CPMVertexPartition = leidenalg$find_partition(pygraph,
								 																							leidenalg$CPMVertexPartition, initial_membership = initial_membership,
								 																							weights = w, seed = seed, n_iterations = n_iterations,
								 																							max_comm_size = max_comm_size,
								 																							node_sizes = node_sizes, resolution_parameter = resolution_parameter),

								 MutableVertexPartition = leidenalg$find_partition(pygraph,
								 																									leidenalg$MutableVertexPartition, initial_membership = initial_membership,
								 																									seed = seed, n_iterations = n_iterations, max_comm_size = max_comm_size),

								 SignificanceVertexPartition = leidenalg$find_partition(pygraph,
								 																											 leidenalg$SignificanceVertexPartition, initial_membership = initial_membership,
								 																											 seed = seed, n_iterations = n_iterations, max_comm_size = max_comm_size,
								 																											 node_sizes = node_sizes, resolution_parameter = resolution_parameter),

								 SurpriseVertexPartition = leidenalg$find_partition(pygraph,
								 																									 leidenalg$SurpriseVertexPartition, initial_membership = initial_membership,
								 																									 weights = w, seed = seed, n_iterations = n_iterations,
								 																									 max_comm_size = max_comm_size, node_sizes = node_sizes),

								 ModularityVertexPartition.Bipartite = run_bipartite_partitioning(pygraph,
								 																																 initial_membership = initial_membership, weights = w,
								 																																 resolution_parameter_01 = resolution_parameter,
								 																																 resolution_parameter_0 = 0, resolution_parameter_1 = 0,
								 																																 degree_as_node_size = TRUE, types = "type", seed = seed,
								 																																 n_iterations = n_iterations),

								 CPMVertexPartition.Bipartite = run_bipartite_partitioning(pygraph,
								 																													initial_membership = initial_membership, weights = w,
								 																													resolution_parameter_01 = resolution_parameter,
								 																													resolution_parameter_0 = 0, resolution_parameter_1 = 0,
								 																													degree_as_node_size = degree_as_node_size, types = "type",
								 																													seed = seed, n_iterations = n_iterations), stop("please specify a partition type as a string out of those documented"))

	res <- paste0("k",formatNumber2Character(part$membership + 1))
	if(returnAsFactor) res<-as.factor(res)
	res
}



#' Title
#'
#' @param umapModel
#' @param n_neighbors
#' @param metric
#'
#' @return
#' @export
#'
#' @examples
getAdjMatfromUMAPmodel<-function(umapModel,n_neighbors=10,metric=NULL){
	if(ncol(umapModel$nn[[metric]]$idx) < n_neighbors){
		stop("The provided umap model contains the data for ",ncol(umapModel$nn[[metric]]$idx),
				 " neighbors, please decrease the n_neighbors parameter or recompute the model on a higher number of neighbors")
	}
	if(is.null(metric)){
		metric<-names(umapModel$nn)[1]
	}else{
		if(!metric %in% names(umapModel$nn)) stop(metric, " distance metric is not computed in the provided model")
	}
	#from 2 --> rm diagonals
	knn_indices<-umapModel$nn[[metric]]$idx[,2:n_neighbors]
	knn_dists<-umapModel$nn[[metric]]$dist[,2:n_neighbors]

	n<-nrow(knn_indices)
	as(Matrix::sparseMatrix(
		i = as.integer(rep(1:n,each=ncol(knn_indices))),
		j = as.integer(as.vector(t(knn_indices))),
		x = as.numeric(as.vector(t(knn_dists))),
		dims = c(n,n)),"dgTMatrix")
}

