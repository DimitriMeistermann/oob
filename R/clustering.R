
#' Compute a Leiden clustering from a UMAP model.
#'
#' @param umapWithNN A list with the UMAP coordinates and the nearest neighbor data as returned by `UMAP` if `ret_nn = TRUE`.
#' @param n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for the Leiden clustering.
#' @param metric Character. One of the metric used for computing the UMAP model. If NULL take the first available one.
#' @param partition_type     Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details)
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
			stop(metric,
					 " distance metric is not computed in the provided model")
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
	graph <-
		ig$Graph$Weighted_Adjacency(matrix = adjMat, mode = mode, ...)
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
#' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details).
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
					 	"SurpriseVertexPartition"
					 )) {
		ig <- import("igraph")
		numpy <- import("numpy", delay_load = TRUE)
		leidenalg <- import("leidenalg", delay_load = TRUE)
		#
		# if (laplacian_init && is.null(weight_parameter)) {
		#     pygraph <- ig$Graph$simplify(pygraph, multiple = TRUE, loops = TRUE)
		#     laplacian <- do.call("cbind",ig$Graph$laplacian(pygraph))
		#     pygraph$es$set_attribute_values("weight",value = -as.matrix(laplacian)[as.matrix(laplacian) <     0])
		# }
		graphIsWeighted <- ig$Graph$is_weighted(pygraph)

		w <- unlist(pygraph$es$get_attribute_values("weight"))

		partition_type <- match.arg(partition_type)
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
			stop(
				"please specify a partition type as a string out of those documented"
			)
		)

		res <-
			paste0("k", formatNumber2Character(part$membership + 1))
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
				stop(metric,
						 " distance metric is not computed in the provided model")
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
				i = as.integer(rep(seq_len(n), each = ncol(knn_indices))),
				j = as.integer(as.vector(t(knn_indices))),
				x = as.numeric(as.vector(t(knn_dists))),
				dims = c(n, n)
			),
			"TsparseMatrix"
		)
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
#' resClust<-hierarchicalClustering(iris[,c(1, 2, 3)],transpose = FALSE)
#' best.cutree(resClust,graph=TRUE)
#' best.cutree(resClust)
#' best.cutree(resClust, loss=TRUE)
#' cutree(resClust,k = best.cutree(resClust))
best.cutree <- function(hc,
												min = 2,
												max = 20,
												loss = FALSE,
												graph = FALSE) {
	if (is(hc, "hclust"))
		hc <- as.hclust(hc)
	max <- min(max, length(hc$height) - 1)
	inert.gain <- rev(hc$height)
	intra <- rev(cumsum(rev(inert.gain)))
	relative.loss = intra[min:(max + 1)] / intra[(min - 1):(max)]
	derivative.loss = relative.loss[seq(2, length(relative.loss), 1)] -
		relative.loss[seq(1, length(relative.loss) - 1, 1)]
	names(derivative.loss) <- min:max
	if (graph) {
		print(
			ggplot(
				data.frame(
					"partition" = min:max,
					"derivative.loss" = derivative.loss
				),
				aes_string(x = "partition", y = "derivative.loss")
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
#' resClust<-hierarchicalClustering(iris[,c(1, 2, 3)],transpose = FALSE)
#' plot(resClust,hang=-1)
#' resClust<-hierarchicalClustering(iris[,c(1, 2, 3)],transpose = TRUE,bootstrap = TRUE,nboot = 20)
#' plot(resClust,hang=-1)

hierarchicalClustering <-
	function (x,
						transpose = TRUE,
						method.dist = "euclidean",
						method.hclust = "ward.D2",
						bootstrap = FALSE,
						nboot = 10,
						PCAfirst = FALSE,
						nDimPCA = NULL)
	{
		if (transpose)
			x <- t(x)
		if (PCAfirst) {
			x <- PCA(x, transpose = FALSE, scale = FALSE)$x
			if (!is.null(nDimPCA)) {
				x <- x[, seq_len(nDimPCA)]
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
		}
		else {
			if (is.function(method.dist)) {
				resDist <- method.dist(x)
			}
			else if (method.dist == "pearson") {
				resDist <- corrDist(x)
			}
			else {
				resDist <- dist(x, method = method.dist)
			}
			resClust <-
				stats::hclust(resDist, method = method.hclust)
		}
		return(resClust)
	}

