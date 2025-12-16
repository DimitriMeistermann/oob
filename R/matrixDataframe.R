#' Scale per row
#'
#' @param data A matrix or dataframe of numerics.
#' @param center either a logical value or numeric-alike vector of length equal
#'   to the number of columns of x, where ‘numeric-alike’ means that
#'   as.numeric(.) will be applied successfully if is.numeric(.) is not true.
#' @param scaled either a logical value or a numeric-alike vector of length
#'   equal to the number of columns of x.
#'
#' @return For scale.default, the centered, scaled matrix. The numeric centering
#'   and scalings used (if any) are returned as attributes "scaled:center" and
#'   "scaled:scale"
#' @export
#'
#' @examples
#' rowScale(matrix(rnorm(100),ncol = 5))
rowScale <- function(data,
                    center = TRUE,
                    scaled = FALSE) {
    data <- t(data)
    data <- t(scale(data, center = center, scale = scaled))
    return(data)
}

#' Print first columns and rows of a matrix/d
#'
#' @param x A matrix or a dataframe
#' @param n number of columns/rows to print
#'
#' @return Print the result in the console.
#' @export
#'
#' @examples
#' data("sampleAnnot")
#' chead(sampleAnnot)
chead <- function(x, n = 5) {
    print(x[seq_len(min(n, nrow(x))), seq_len(min(n, ncol(x)))])
}

#' Compute correlation distance
#'
#' @param x A matrix of numeric. The function will compute the distance between
#'   the rows (same as `dist`).
#' @param method Character. Name of the method to compute correlation. See
#'   `method` argument from `cor`.
#' @param ... Other arguments to be passed to `cor`.
#'
#' @return An object of class "dist".
#' @export
#'
#' @examples
#' data("iris")
#' corrDist(t(iris[,seq_len(3)]))
corrDist <- function(x, method = "pearson", ...) {
    x <- Matrix::t(x)
    if(requireNamespace("WGCNA", quietly = TRUE)) {
        x <- WGCNA::cor(x, method = method)
    } else {
        x <- stats::cor(x, method = method)
        message("Consider installing WGCNA for better performance when ",
                "computing correlation.")
    }
    return(as.dist( (1 - x)/2 ))
}

#' Compute cosine distance
#'
#' @param x  A matrix of numeric. The function will compute the distance between
#'   the rows (same as `dist`).
#'
#' @return An object of class "dist".
#' @export
#'
#' @examples
#' data("iris")
#' cosineDist(t(iris[,seq_len(3)]))
cosineDist<-function(x){
    1 - lsa::cosine(
        Matrix::t(x)
    ) |> as.dist()
}


#' Compute covariance distance
#'
#' @param x  A matrix of numeric. The function will compute the distance between
#'   the rows (same as `dist`).
#'
#' @return  An object of class "dist".
#' @export
#'
#' @examples
#' data("iris")
#' covDist(t(iris[,seq_len(3)]))
covDist<-function(x){
    x<-cov(
        Matrix::t(x)
    )
    max(x)-x |> as.dist()
}


#' Convert a list of key from an ID to another.
#'
#' @param keyList A vector of character. The list of key from an ID to be
#'   converted (for example gene symbols).
#' @param tabKey  A dataframe or a matrix of character. The database that
#'   contains correspondence between each key from each ID.
#' @param colOldKey Integer. Column of `tabKey` that contains the old (same as
#'   keyList) IDs.
#' @param colNewKey Integer. Column of `tabKey` that contains the new IDs.
#'
#' @return A vector of character.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("humanGeneIDtable")
#' geneSym<-rownames(bulkLogCounts)
#' geneEntrez<-ConvertKey(geneSym,tabKey = humanGeneIDtable,
#'     colOldKey = "SYMBOL",colNewKey = "ENTREZID")
ConvertKey <- function(keyList,
                    tabKey,
                    colOldKey = 1,
                    colNewKey = 2) {
    hashCorr <- tabKey[, colNewKey]
    names(hashCorr) <- tabKey[, colOldKey]
    returned <- hashCorr[keyList]
    names(returned) <- NULL
    return(as.character(returned))
}

#' Convert rownames from a matrix/df from one ID type to another.
#'
#' @param tab A matrix or a dataframe. Rownames are corresponding to keys from
#'   an ID to be converted
#' @param tabKey A dataframe or a matrix of character. The database that
#'   contains correspondence between each key from each ID.
#' @param colOldKey Integer. Column of `tabKey` that contains the old (same as
#'   keyList) IDs.
#' @param colNewKey Integer. Column of `tabKey` that contains the new IDs.
#' @param first If for given the old ID has several correspondence, take the
#'   first one as the new ID.
#' @param dim Interger, 1 or 2. If 1 convert rows, else columns.
#' @param fun If `first=FALSE`, A function to compute the summary statistics
#'   which can be applied to all data that correspond to ther same ID.
#'
#' @return A matrix or dataframe. Same format as tab.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("humanGeneIDtable")
#' bulkLogCountsEntrez<-ConvertKeyMatrix(bulkLogCounts,
#'     tabKey = humanGeneIDtable, colOldKey = "SYMBOL",colNewKey = "ENTREZID")
ConvertKeyMatrix <- function (tab, tabKey, colOldKey = 1, colNewKey = 2, first = TRUE,
															dim = 1, fun) {
	if (dim == 2)
		tab <- t(tab)
	keyList <- rownames(tab)
	newkey <- ConvertKey(keyList, tabKey, colOldKey, colNewKey)
	names(newkey) <- as.character(seq_along(newkey))
	if (first) {
		newkey <- newkey[which(!is.na(newkey))]
		newkey <- takefirst(newkey)
		tab <- tab[as.numeric(names(newkey)), ]
		rownames(tab) <- newkey
	}
	else {
		tab <- tab[which(!is.na(newkey)), ]
		newkey <- newkey[which(!is.na(newkey))]
		tab <- aggregate(x = tab, by = list(newkey), FUN = fun)
		rownames(tab) <- tab[, 1]
		tab <- as.matrix(tab[, -1])
	}
	if (dim == 2)
		tab <- t(tab)
	return(tab)
}


#' Delete row/column in a matrix/df with NA names
#'
#' @param x A dataframe or matrix
#' @param side 1:row names are treated, 2: column names.
#'
#' @return A matrix or dataframe. Same format as x.
#' @export
#'
#' @examples
#' x<-matrix(1:25,ncol=5)
#' rownames(x)<-c(letters[seq_len(4)],NA)
#' supprNAnames(x)
supprNAnames <- function(x , side = 1) {
    #Vector case
    if (is.null(dim(x))) {
        return(x[which(!is.na(names(x)))])
    }
    #Col case
    if (side == 2) {
        return(x[, which(!is.na(colnames(x)))])
    }
    #Row case
    else{
        return(x[which(!is.na(rownames(x))), ])
    }
}

#' Return the row and column index (2D coordinate) from a 1D coordinate in a
#' matrix.
#'
#' @param x An integer. 1D coordinate to be transposed in 2D.
#' @param Mat A matrix.
#'
#' @return A vector of 2 integer. Index of row and column.
#' @export
#'
#' @examples
#' test<-matrix(1:25,ncol=5)
#' test[10]
#' matrixCoord1D_2D(10,test)
matrixCoord1D_2D <- function(x, Mat) {
    matDim <- dim(Mat)
    c((x - 1) %% matDim[1] + 1 , (x - 1) %/% matDim[1] + 1)
}


#' Create empty matrix from a column and rownames
#'
#' @param row A vector of character.
#' @param col A vector of character.
#' @param value Default value for each element of the matrix.
#'
#' @return A matrix. Size of the matrix will be `length(row)*length(col)`.
#' @export
#'
#' @examples
#' matrixFromDimnames(1:10,letters[1:5])
matrixFromDimnames <- function(row, col, value = 0) {
    matrix(
        value,
        ncol = length(col),
        nrow = length(row),
        dimnames = list(row, col)
    )
}


#' Put rows of a matrix on the same range than another.
#'
#' @param matToAdjust Matrix of numeric where the range of each row has to be
#'   adjusted.
#' @param matGoodRange Matrix of numeric, same dimension as `matToAdjust`.
#'
#' @return Return `matToAdjust` but the each row has now same minimum and
#'   maximum to its corresponding row in matGoodRange
#' @export
#'
#' @examples
#' m1<-matrix(rnorm(100),ncol=50)
#' qplotDensity(m1[1,])
#' m2<-matrix(rnorm(100,20),ncol=50)
#' qplotDensity(m2[1,])
#' m1Rescale<-reScale(m1,m2)
#' qplotDensity(m1Rescale[1,])
reScale <- function(matToAdjust, matGoodRange) {
    apply(matGoodRange, 1, function(x)
        max(x) - min(x)) /
        (apply(matToAdjust, 1, function(x)
            max(x) - min(x))) *
        (matToAdjust - apply(matToAdjust, 1, max)) + apply(matGoodRange, 1, max)
}


#' Format a dataframe following instructions from another.
#'
#' @param annotDataFrame A dataframe of features from any kind.
#' @param metaAnnot A dataframe with at least a column "Type" containing the
#'   feature variable type and "colorScale" containing eventually the mapping of
#'   colors to values. Each row is named by a feature (column) from
#'   `annotDataFrame`. If type of `annotDataFrame` feature and column `Type`
#'   mismatch, the feature will be automatically converted. `colorScale` Contain
#'   a character vector describing the colors that has to be mapped to values,
#'   in the form "color1,color2,...,colorN". Can be an empty string if no color
#'   scale has to be exported. If the feature is a factor, the colors must be
#'   named by the possible values:
#'   "level1=color1,level2=color2,...,levelN=colorN". In this case, the factor
#'   levels order will correspond to the level order of the color scale.
#'
#' @return `annotDataFrame`, but Variable type are eventually converted
#'   following Type coloumn from `metaAnnot`. Has also an attribute
#'   "colorScales" which contains a list named by the features that had a
#'   described color scale in `metaAnnot`. Each element contains the colors at
#'   breaks for continuous values. In the case of factors, the colors are named
#'   to their corresponding level.
#' @export
#'
#' @examples
#' df <- data.frame(a = 1:12, b = rep(c("1", "19", "2"), each = 4))
#' metaDf <-
#'     data.frame(
#'         Type = c("numeric", "factor"),
#'         colorScale = c("", "1=blue,2=white,19=red"),
#'         row.names = c("a", "b")
#'     )
#' res <- formatAnnotFromMeta(df, metaDf)
#' res$b
#' attr(res, "colorScale")
formatAnnotFromMeta <- function(annotDataFrame, metaAnnot) {
	colorScales <- list()
	for (feature in rownames(metaAnnot)) {
		if (feature %in% colnames(annotDataFrame)) {
			type <- metaAnnot[feature, "Type"]
			if (!(type == "" | is.na(type))) {
				annotDataFrame[, feature] <-
					do.call(paste0("as.", type), list(x = annotDataFrame[, feature]))
			}
			if (metaAnnot[feature, "colorScale"] != "") {
				colorScale <-
					strsplit(str_remove_all(metaAnnot[feature, "colorScale"], " "), split = ",")[[1]]
				splitted <-
					vapply(colorScale, function(x)
						strsplit(x, "=")[[1]], FUN.VALUE = character(2))
				colorScales[[feature]] <- splitted[2, ]
				names(colorScales[[feature]]) <- splitted[1, ]
				if (metaAnnot[feature, "Type"] == "factor") {
					annotDataFrame[, feature] <-
						factor(annotDataFrame[, feature], levels = names(colorScales[[feature]]))
				}
			}
		}
	}
	attr(annotDataFrame, "colorScales") <- colorScales
	annotDataFrame
}


#' Convert factors in a dataframe to strings
#'
#' @param dt A dataframe containing factor column
#'
#' @return A dataframe where all the factor column have been converted to srings
#' @export
#' @examples
#' a<-data.frame(x=factor(c("a","a","b")),y=seq_len(3))
#' b<-factorAsStrings(a)
#' b$x
factorAsStrings <- function(dt) {
    for (i in seq_len(ncol(dt))) {
        if (is.factor(dt[, i]))
            dt[, i] <- as.character(dt[, i])
    }
    return(dt)
}

#' Aggregate columns or rows of a matrix by a vector of factor
#'
#' @param x A matrix of numeric.
#' @param by A vector of factor. Same size as number of row in `x` if  `byRow`,
#'   otherwise same size as number of columns.
#' @param FUN a function to compute the summary statistics which can be applied
#'   to all data subsets.
#' @param byRow Logical. Do the aggregation by row or by column. If NULL,
#'   determines automatically based on the size of `by`.
#'
#' @return A matrix of numeric. Aggregated dimension names are changed to match
#'   the levels of `by`.
#' @export
#'
#' @examples
#' data(iris)
#' aggregMatPerVector(iris[,seq_len(4)],iris$Species)
aggregMatPerVector <- function(x, by, FUN = mean, byRow = NULL) {
	if(is.data.frame(x)) x <- as.matrix(x)
	stopifnot(is.matrix(x) || inherits(x, "Matrix"))
	nr <- nrow(x); nc <- ncol(x)

	if (is.null(byRow)) {
		if (length(by) == nr) byRow <- TRUE
		else if (length(by) == nc) byRow <- FALSE
		else stop("`by` length must match nrow(x) or ncol(x).")
	}

	# Prepare grouping as an integer factor (cheap & consistent)
	f <- factor(by, exclude = NULL)                # keep NAs separate if present
	g <- as.integer(f)
	labs <- levels(f)
	k <- nlevels(f)
	cnt <- tabulate(g, nbins = k)

	is_sparse <- inherits(x, "dgCMatrix")

	if (byRow) {
		if (is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean") ||
											identical(FUN, sum)  || identical(FUN, base::sum)  || identical(FUN, "sum"))) {
			# Sparse row-group sum via G %*% x
			G <- Matrix::sparseMatrix(i = g, j = seq_len(nr), x = 1L, dims = c(k, nr))
			out <- G %*% x                              # (k x nr) %*% (nr x nc) -> (k x nc)
			if (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean")) {
				out <- Matrix::Diagonal(x = 1 / pmax(cnt, 1L)) %*% out
			}
			rownames(out) <- labs
			return(as.matrix(out))
		} else if (!is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean"))) {
			out <- rowsum(x, group = g, reorder = FALSE)  # dense, fast C code
			out <- out / pmax(cnt, 1L)
			rownames(out) <- labs
			return(as.matrix(out))
		} else if (!is_sparse && (identical(FUN, sum) || identical(FUN, base::sum) || identical(FUN, "sum"))) {
			out <- rowsum(x, group = g, reorder = FALSE)
			rownames(out) <- labs
			return(as.matrix(out))
		} else {
			# Generic FUN (slower): apply per group without copying whole matrix
			idx <- split(seq_len(nr), g, drop = TRUE)
			out <- vapply(idx, function(ix) apply(x[ix, , drop = FALSE], 2, FUN), numeric(nc))
			out <- t(out)
			rownames(out) <- labs
			return(as.matrix(out))
		}
	} else { # group columns
		if (is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean") ||
											identical(FUN, sum)  || identical(FUN, base::sum)  || identical(FUN, "sum"))) {
			# Sparse col-group sum via x %*% t(G)
			Gt <- Matrix::sparseMatrix(i = seq_len(nc), j = g, x = 1L, dims = c(nc, k)) # this is t(G)
			out <- x %*% Gt                           # (nr x nc) %*% (nc x k) -> (nr x k)
			if (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean")) {
				out <- Matrix::t(Matrix::t(out) / pmax(cnt, 1L))  # sweep-like without dense copy
			}
			colnames(out) <- labs
			return(as.matrix(out))
		} else if (!is_sparse && (identical(FUN, mean) || identical(FUN, base::mean) || identical(FUN, "mean"))) {
			# Dense columns: avoid data.frames; one transpose (fast BLAS) is usually OK
			out <- t(rowsum(t(x), g, reorder = FALSE))   # sum per column-group
			out <- sweep(out, 2, pmax(cnt, 1L), "/")
			colnames(out) <- labs
			return(as.matrix(out))
		} else if (!is_sparse && (identical(FUN, sum) || identical(FUN, base::sum) || identical(FUN, "sum"))) {
			out <- t(rowsum(t(x), g, reorder = FALSE))
			colnames(out) <- labs
			return(as.matrix(out))
		} else {
			idx <- split(seq_len(nc), g, drop = TRUE)
			out <- vapply(idx, function(ix) apply(x[, ix, drop = FALSE], 1, FUN), numeric(nr))
			colnames(out) <- labs
			return(as.matrix(out))
		}
	}
}
#' Draw n samples from each population and return it has a named vector
#'
#' @param sampleNames A character vector of the name of the samples.
#' @param group A factor or a character vector of the same length as
#'   `sampleNames`. Describe the population of each sample.
#' @param maxDrawSize Maximum number of observation to draw per group.
#' @param minDrawSize Minimum number of observation to draw per group.
#' @param replace Logical. Should the sample be drawn with replacement.
#'
#' @return A named vector with `sampleNames` with the population of each sample.
#' @export
#'
#' @examples
#' sampleNames<-paste0("sample",1:20)
#' group<-c(rep("A",5),rep("B",12),rep("C",3))
#' drawSamplePerGroup(sampleNames,group)
#' drawSamplePerGroup(sampleNames,group,minDrawSize=6)
#' drawSamplePerGroup(sampleNames,group,maxDrawSize=2)
#'
#' # If minDrawSize > number of sample in a group and replace=TRUE, samples
#' # will be drawn with replacement in this group
#' drawSamplePerGroup(sampleNames,group,minDrawSize=6, replace=TRUE)
drawSamplePerGroup<-function(sampleNames, group, maxDrawSize = NULL,
                            minDrawSize = NULL,replace=FALSE){
    if (is.null(group))
        return(sampleNames)
    if (is.null(names(group)))
        names(group)<-sampleNames
    if(! (is.factor(group) | is.character(group)))
        stop("group must be a factor or a character vector")
    drawSize <- min(table(group))
    if (!is.null(maxDrawSize))
        drawSize <- min(drawSize, maxDrawSize)
    if (!is.null(minDrawSize))
        drawSize <- max(drawSize, minDrawSize)
    drawCells <- unlist(lapply(unique(group), function(lvl) {
        cells <- sampleNames[group == lvl]
        if(replace & length(cells)<drawSize){
            return(sample(cells, drawSize, replace = TRUE))
        } else {
            sample(cells, min(drawSize, length(cells)))
        }
    }))
    return(group[drawCells])
}





