#' Compute a Leiden clustering from a UMAP model.
#'
#' @param umapWithNN A list with the UMAP coordinates and the nearest neighbor
#'   data as returned by `UMAP` if `ret_nn = TRUE`.
#'   This can also be a `SingleCellExperiment` object with where the UMAP
#'   function has been run with `ret_nn = TRUE`.
#' @param n_neighbors The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for the Leiden clustering.
#' @param metric Character. One of the metric used for computing the UMAP model.
#'   If NULL take the first available one.
#' @param partition_type     Type of partition to use. Defaults to
#'   RBConfigurationVertexPartition. Options include: ModularityVertexPartition,
#'   RBERVertexPartition, CPMVertexPartition, MutableVertexPartition,
#'   SignificanceVertexPartition, SurpriseVertexPartition,
#'   ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the
#'   Leiden python module documentation for more details)
#' @param returnAsFactor If TRUE return the clusters attributions as a factor
#'   vector and not as characters.
#' @param n_iterations Number of iterations. If the number of iterations is
#'   negative, the Leiden algorithm is run until an iteration in which there was
#'   no improvement.
#' @param resolution_parameter A parameter controlling the coarseness of the
#'   clusters.
#' @param seed Seed for the random number generator. By default uses a random
#'   seed if nothing is specified.
#' @param initial_membership Initial membership for the partition. If `NULL`
#'   then defaults to a singleton partition.
#' @param max_comm_size Maximal total size of nodes in a community. If zero (the
#' @param prefix Prefix for the cluster
#'   default), then communities can be of any size.
#' @return A vector of character or factor if `returnAsFactor`, same length as
#'   number of samples in the UMAP. Cluster attribution of samples.
#'   If `umapWithNN` is a `SingleCellExperiment` object, the function will add
#'   the cluster attribution to the colData of the object and return the object.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' umapWithNN <- UMAP(bulkLogCounts, ret_nn = TRUE)
#' proj2d(umapWithNN$embedding, colorBy = leidenFromUMAP(umapWithNN))
#' sce <- SingleCellExperiment(assays = list(counts = bulkLogCounts))
#' sce <- UMAP(sce, ret_nn = TRUE)
#' sce<-leidenFromUMAP(sce)
#' proj2d(sce, colorBy = "leiden_cluster")
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
                    initial_membership = NULL,
                    max_comm_size = 0L,
										prefix = "k") {
    sce_obj <-NULL
    if (inherits(umapWithNN, "SingleCellExperiment")) {
        sce_obj <- umapWithNN
        if(!is.null(metadata(sce_obj)$UMAP)){
                umapWithNN <- metadata(sce_obj)$UMAP
        }else{
                stop("UMAP model not found in the provided",
                    "SingleCellExperiment object, please run oob::UMAP first.")
        }
    }
    if (is.null(metric)) {
        metric <- names(umapWithNN$nn)[1]
    } else {
        if (!metric %in% names(umapWithNN$nn)) {
            stop(
                metric,
                " distance metric is not computed in the provided model"
            )
        }
    }

    partition_type <- match.arg(partition_type)
    adjMat <- getAdjMatfromUMAPWithNN(
        umapWithNN,
        n_neighbors = n_neighbors,
        metric = metric
    )

    proc <- basiliskStart(pythonEnv)
    on.exit(basiliskStop(proc))
    res<-basiliskRun(proc, function(adjMat, mode,returnAsFactor, n_iterations,
        resolution_parameter, seed, initial_membership,
        max_comm_size, partition_type){
            ig <- import("igraph")

            pygraph <- ig$Graph$Weighted_Adjacency(matrix = adjMat, mode = mode)

            leidenFromPygraph(
                    pygraph,
                    returnAsFactor = returnAsFactor,
                    n_iterations = n_iterations,
                    resolution_parameter = resolution_parameter,
                    seed = seed,
                    initial_membership = initial_membership,
                    max_comm_size = max_comm_size,
                    partition_type = partition_type,
                    prefix = prefix)
        }, adjMat, mode = "undirected", returnAsFactor = returnAsFactor,
                n_iterations = n_iterations,
                resolution_parameter = resolution_parameter,
                seed = seed, initial_membership = initial_membership,
                max_comm_size = max_comm_size,
                partition_type = partition_type)

    if(!is.null(sce_obj)){
        colData(sce_obj)$leiden_cluster <- res
        return(sce_obj)
    }else{
        return(res)
    }
}


leidenFromPygraph <-
    function(pygraph,
            returnAsFactor = FALSE,
            n_iterations = -1,
            resolution_parameter = .5,
            seed = 666,
            initial_membership = NULL,
            max_comm_size = 0L,
            node_sizes = NULL,
    				prefix = "k",
            partition_type = c(
                "RBConfigurationVertexPartition",
                "ModularityVertexPartition",
                "RBERVertexPartition",
                "CPMVertexPartition",
                "MutableVertexPartition",
                "SignificanceVertexPartition",
                "SurpriseVertexPartition"
            )){

            ig <- import("igraph")
            numpy <- import("numpy")
            leidenalg <- import("leidenalg")
            graphIsWeighted <- ig$Graph$is_weighted(pygraph)

            w <- unlist(pygraph$es$get_attribute_values("weight"))

            partition_type <- match.arg(partition_type)
            if (!is.null(seed)) {
                seed <- as.integer(seed)
            }
            if (!is.integer(n_iterations)) {
                n_iterations <- as.integer(n_iterations)
            }
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
                    "please specify a partition type ",
                    "as a string out of those documented"
                )
            )
        res<-paste0(prefix, formatNumber2Character(part$membership + 1))
        if (returnAsFactor) {
            res <- as.factor(res)
        }
        res
    }


#' Compute an adjacency matrix from the nearest neighbor matrix of a UMAP model
#'
#' @param umapWithNN A UMAP model as returned by UMAP if `ret_nn = TRUE`,
#' @param n_neighbors The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for the Leiden clustering.
#' @param metric  Character. One of the metric used for computing the UMAP
#'   model. If NULL take the first available one.
#'
#' @return A sparse adjacency matrix from the class `dgTMatrix`.
getAdjMatfromUMAPWithNN<- function (umapWithNN, n_neighbors = 10, metric = NULL)
{
	if (is.null(metric)) {
		metric <- names(umapWithNN$nn)[1]
	}
	else {
		if (!metric %in% names(umapWithNN$nn)) {
			stop(metric, " distance metric is not computed in the provided model")
		}
	}
	full_idx <- umapWithNN$nn[[metric]]$idx
	full_dist <- umapWithNN$nn[[metric]]$dist
	if (ncol(full_idx) < n_neighbors) {
		stop("The provided umap model contains the data for ",
				 ncol(full_idx), " neighbors, please decrease the n_neighbors parameter ",
				 "or recompute the model on a higher number of neighbors")
	}
	n <- nrow(full_idx)
	idx <- full_idx[, seq_len(n_neighbors), drop = FALSE]
	dist <- full_dist[, seq_len(n_neighbors), drop = FALSE]

	# uwot::umap(..., ret_nn=TRUE) on a self-fit *usually* returns each point as its own
	# nearest neighbor at distance 0 (typically column 1, since neighbors are sorted by
	# ascending distance) -- but not reliably for every row: uwot's default approximate
	# nearest-neighbor search doesn't have perfect recall, and exact-duplicate rows tie at
	# distance 0 with no guaranteed ordering. An earlier version of this function assumed
	# column 1 was always self and hard-errored otherwise -- too strict in practice (hit on
	# real data where self wasn't found in column 1 for every row). Instead, self-matches
	# are dropped from wherever they actually fall within the requested n_neighbors columns
	# (same approach as HybridClustR::buildGraphFromNN, cross-checked against this code):
	# rows where self is found lose exactly that one column; rows where it isn't (or where
	# self happens to sit at a rank not covered by n_neighbors -- indistinguishable from
	# "not among the closest n_neighbors points", which is disqualifying either way) simply
	# keep all n_neighbors columns as their real neighbors. This also means fitting UMAP
	# with n_neighbors=k no longer requires requesting k+1 to leave headroom.
	keep <- idx != seq_len(n)   # n x n_neighbors logical, seq_len(n) recycled column-wise

	# Previous code (`knn_dists / knn_dists`) collapsed every real distance to a
	# meaningless constant 1 and produced NaN for any included self-loop -- those NaNs
	# landed directly in the weights passed to leidenalg::find_partition(). Self-loops are
	# excluded above instead of ever being weighted; inverse-distance similarity here keeps
	# edges reflecting actual UMAP distance (closer neighbors get a higher weight) rather
	# than a flat constant.
	weight <- 1 / (1 + dist[keep])
	as(Matrix::sparseMatrix(i = row(idx)[keep], j = idx[keep], x = weight,
													 dims = c(n, n)), "TsparseMatrix")
};assignInNamespace("getAdjMatfromUMAPWithNN",getAdjMatfromUMAPWithNN,"oob")



#' Determine the best partition in a hierarchical clustering
#'
#' @param hc A hclust object.
#' @param min Minimum number of class in the partition.
#' @param max Maximum number of class in the partition
#' @param loss Logical. Return the list of computed partition with their
#'   derivative loss.
#' @param graph Logical. Plot a graph of computed partition with their
#'   derivative loss.
#'
#' @details Based on the higher relative loss of inertia, fucntion modified from
#'   the [JLutils](https://rdrr.io/github/larmarange/JLutils/src/R/clustering.R)
#'   package.
#'
#' @return A single integer (best partition) or print a graph if `graph` or a
#'   vector of numeric if `loss`.
#' @export
#'
#' @examples
#' data(iris)
#' resClust <- hierarchicalClustering(iris[, c(1, 2, 3)], transpose = FALSE)
#' best.cutree(resClust, graph = TRUE)
#' best.cutree(resClust)
#' best.cutree(resClust, loss = TRUE)
#' cutree(resClust, k = best.cutree(resClust))
best.cutree <- function(hc,
                        min = 2,
                        max = 20,
                        loss = FALSE,
                        graph = FALSE) {
    if (is(hc, "hclust")) {
        hc <- as.hclust(hc)
    }
    max <- min(max, length(hc$height) - 1)
    inert.gain <- rev(hc$height)
    intra <- rev(cumsum(rev(inert.gain)))
    relative.loss <- intra[min:(max + 1)] / intra[(min - 1):(max)]
    derivative.loss <-
        relative.loss[seq(2, length(relative.loss), 1)] -
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
        if (loss) {
            derivative.loss
        } else {
            as.numeric(names(which.max(derivative.loss)))
        }
    }
}


#' Download transcription factors of a species from the JASPAR database
#' @param taxID NCBI Taxonomy ID
#'
#' @returns A vector of TF genes
#' @export
#'
#' @examples
#' TFsOfHuman <- dl_JAFAR_geneSyms()
#' head(TFsOfHuman)
dl_JASPAR_geneSyms <- function(taxID = 9606) {
	if (!(requireNamespace("JASPAR2024", quietly = TRUE))) {
		stop("Error, please install 'JASPAR2024' package, then retry")
	}
	JASPAR2024 <- JASPAR2024::JASPAR2024()
	JASPARConnect <- dbConnect(SQLite(), JASPAR2024::db(JASPAR2024))
	geneNames <- dbGetQuery(
		JASPARConnect,
		paste0(
			"SELECT DISTINCT m.name FROM matrix AS m
	  JOIN matrix_species AS ms ON ms.id = m.id
	  WHERE ms.tax_id = ",
			taxID,
			";"
		)
	)[, "NAME"]
	strsplit(geneNames, split = ":", fixed = T) |> unlist() |> unique()
}


#' Compute gene modules with UMAP + Leiden clustering
#'
#' Clusters genes into co-expression modules by (i) building a kNN graph with
#' `UMAP()` (with `ret_nn = TRUE`) and (ii) running Leiden community detection via
#' `leidenFromUMAP()`. For each retained module, a per-sample activation score is
#' computed with `activeScorePCAlist()`. Modules smaller than `min_module_size`
#' are collapsed into `"M0"`.
#'
#' If `TF_db` is provided, module labels can be renamed with `renameModulePerTF()`
#' (assumed to relabel factor levels in `gene_annot$Module`).
#'
#' @param data A numeric matrix-like object with **genes in rows** and
#'   **samples/cells in columns**, or a `SingleCellExperiment` object.
#' @param transpose Logical; if `TRUE`, `data` is assumed to be samples/cells x
#'   genes and will be transposed internally.
#' @param metric Character; distance metric forwarded to `UMAP()` (e.g. `"cosine"`,
#'   `"correlation"`).
#' @param resolution_parameter Numeric; Leiden resolution parameter forwarded to
#'   `leidenFromUMAP()`.
#' @param k Integer; number of nearest neighbors used by UMAP and Leiden.
#' @param TF_list Optional; forwarded to `renameModulePerTF()` to rename module labels.
#' @param sce_assay Character; assay name to extract when `data` is a
#'   `SingleCellExperiment`.
#' @param min_module_size Integer; minimum number of genes to keep a module
#'   (otherwise assigned to `"M0"`).
#' @param prefix Character; prefix used by `leidenFromUMAP()` when naming modules.
#' @param seed Optional integer; if not `NULL`, used for `set.seed()` before
#'   running UMAP/Leiden.
#' @param return_all Logical; if `TRUE`, returns a list of results; if `FALSE`,
#'   returns only `gene_annot`.
#'
#' @return If `return_all = TRUE`, a list with:
#' \describe{
#'   \item{gene_annot}{`data.frame` with per-gene module assignment (`Module`),
#'     membership (correlation with module score), and contribution.}
#'   \item{module_genes}{Named list of character vectors; genes per module
#'     (includes `"M0"` as empty).}
#'   \item{activation_score}{Numeric matrix; module activation scores (rows are
#'     samples/cells).}
#'   \item{membership}{Numeric matrix; gene-by-module correlations.}
#'   \item{params}{List of parameters used.}
#' }
#' If `return_all = FALSE`, only `gene_annot` is returned.
#'
#' @seealso `UMAP()`, `leidenFromUMAP()`, `activeScorePCAlist()`, `renameModulePerTF()`
# not exported yet #' @export
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(2000), nrow = 100, ncol = 20,
#'               dimnames = list(paste0("G", 1:100), paste0("C", 1:20)))
#' res <- computegenemodules(mat, k = 10, resolution_parameter = 1)
#' chead(res$gene_annot)
computegenemodules <- function(data,
															 transpose = FALSE,
															 metric = "cosine",
															 resolution_parameter = 1.5,
															 k = 20,
															 TF_list = NULL,
															 sce_assay = "logcounts",
															 min_module_size = 10,
															 prefix = "M",
															 seed = NULL,
															 return_all = TRUE) {
	# Extract assay if needed
	sce_obj <- NULL
	if (inherits(data, "SingleCellExperiment")) {
		sce_obj <- data
		# assay() is defined/exported by SummarizedExperiment -- SingleCellExperiment only
		# inherits it (extends SummarizedExperiment) without re-exporting it under its own
		# namespace, so SingleCellExperiment::assay(...) fails to resolve.
		data <- SummarizedExperiment::assay(sce_obj, sce_assay)
	}
	if(is.data.frame(data)) data <- as.matrix(data)
	# Basic checks
	if (!is.matrix(data) && !inherits(data, "Matrix"))
		stop("`data` must be a matrix-like object or a SingleCellExperiment.")
	# is.numeric() is a base-R check that doesn't recognize Matrix-package S4 objects (e.g.
	# dgCMatrix) as numeric even when they are -- contradicting the check just above, which
	# already explicitly allows them through. inherits(data, "dMatrix") is the Matrix
	# package's own equivalent for its numeric-storage matrix classes.
	if (!is.numeric(data) && !inherits(data, "dMatrix"))
		stop("`data` must be numeric.")
	if (!isTRUE(all.equal(k, as.integer(k))) || k <= 1)
		stop("`k` must be an integer > 1.")
	if (!is.numeric(resolution_parameter) || length(resolution_parameter) != 1L || resolution_parameter <= 0)
		stop("`resolution_parameter` must be a single positive number.")
	if (!isTRUE(all.equal(min_module_size, as.integer(min_module_size))) || min_module_size < 1)
		stop("`min_module_size` must be an integer >= 1.")

	# Standardize orientation: genes x samples/cells
	if (isTRUE(transpose)) data <- t(data)

	if (is.null(rownames(data))) rownames(data) <- paste0("gene_", seq_len(nrow(data)))

	if (!is.null(seed)) set.seed(seed)

	# Gene clustering
	gene_cluster <- UMAP(
		data,
		metric = metric,
		transpose = FALSE,
		ret_nn = TRUE,
		n_neighbors = k
	) |>
		leidenFromUMAP(
			n_neighbors = k,
			resolution_parameter = resolution_parameter,
			prefix = prefix
		)

	names(gene_cluster) <- rownames(data)
	# Keep only modules of sufficient size; others -> M0
	n_gene_per_module <- table(gene_cluster)
	kept_modules <- names(n_gene_per_module)[n_gene_per_module >= min_module_size]

	module_genes <- lapply(kept_modules, function(m) names(gene_cluster)[gene_cluster == m])
	names(module_genes) <- kept_modules

	# Activation scores (rows = samples/cells)
	# activScorePC1list is internal to GSDS (not exported, not in oob's NAMESPACE imports --
	# only the singular activScorePC1 is) -- reference it directly rather than via a name
	# that was never actually resolvable from here.
	activ_dt <- GSDS:::activScorePC1list(data, geneList = module_genes)
	activation_score <- activ_dt$activScoreMat

	# Per-gene annotation
	gene_annot <- data.frame(row.names = rownames(data))
	gene_annot$Module <- "M0"
	gene_annot[names(gene_cluster), "Module"] <- gene_cluster
	gene_annot[!gene_annot$Module %in% kept_modules, "Module"] <- "M0"
	# "M0" (the too-small-to-keep catch-all, assigned above) is never in kept_modules /
	# names(module_genes) -- include it explicitly so those genes don't become NA here.
	gene_annot$Module <- factor(gene_annot$Module, levels = c(names(module_genes), "M0"))



	# Membership (gene-module correlation)
	# activation_score (activScoreMat from GSDS:::activScorePC1list) is documented and
	# returned as modules x samples -- stats::cor(x, y) correlates columns pairwise and
	# needs nrow(x)==nrow(y), so it must be transposed to samples x modules to line up
	# with t(data) (samples x genes) here. Produces genes x modules, matching how
	# membership[genes, m] is indexed below.
	membership <- suppressWarnings(stats::cor(
		t(as.matrix(data)),
		t(activation_score),
		use = "pairwise.complete.obs"
	))

	gene_annot$Membership <- NA_real_
	gene_annot$Contribution <- NA_real_

	# Fill membership + contribution for retained modules
	for (m in kept_modules) {
		genes <- module_genes[[m]]
		gene_annot[genes, "Membership"] <- membership[genes, m]
		contrib <- activ_dt$contributionList[[m]]
		if (!is.null(contrib)) {
			gene_annot[genes, "Contribution"] <- unname(contrib[genes])
		}
	}

	# Optional renaming based on TF database (assumed to relabel factor levels)
	if (!is.null(TF_list)) {
		old_levels <- levels(gene_annot$Module)
		gene_annot <- renameLevelsByBestVal(gene_annot, TF_list)
		new_levels <- levels(gene_annot$Module_newName)

		if (length(new_levels) == length(old_levels)) {
			colnames(activation_score) <- new_levels
			colnames(membership) <- new_levels
			names(module_genes) <- new_levels
		}
	}

	if (!isTRUE(return_all)) return(gene_annot)

	list(
		gene_annot = gene_annot,
		module_genes = module_genes,
		activation_score = activation_score,
		membership = membership,
		params = list(
			metric = metric,
			resolution_parameter = resolution_parameter,
			k = k,
			sce_assay = sce_assay,
			min_module_size = min_module_size,
			prefix = prefix,
			seed = seed,
			transposed_input = transpose,
			input_was_sce = !is.null(sce_obj)
		)
	)
}

