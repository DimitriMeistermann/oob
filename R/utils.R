#' Compute geometrical mean
#'
#' @param x numeric vector
#' @param keepZero get rid of 0 before computing.
#' @return A single numeric value.
#' @export
#' @examples
#' gmean(c(1,2,3))
#' gmean(c(0,2,3),keepZero = TRUE)
#' gmean(c(0,2,3),keepZero = FALSE)
gmean <- function(x, keepZero = TRUE) {
    #geometrical mean
    if (sum(x) == 0)
        return(0)
    if (!keepZero) {
        x <- x[x != 0]
    } else{
        if (length(which(x == 0)) > 0)
            return(0)
    }
    return(exp(sum(log(x)) / length(x)))
}

#' Get the precise random seed state
#' @export
#' @return A vector of integers representing the random seed state
#' @export
#' @seealso setRandState
#' @examples
#' randSate<-getRandState()
#' rnorm(5)
#' setRandState(randSate)
#' rnorm(5)
getRandState <- function() {
    # Using `get0()` here to have `NULL` output in case object doesn't exist.
    # Also using `inherits = FALSE` to get value exactly from global environment
    # and not from one of its parent.
    get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}




#' Rename module levels by the best-scoring gene
#'
#' For each level in a module/cluster column, this function picks the gene with
#' the highest value in `interestCol` among `candidates` (restricted to genes
#' present in `rownames(df)`). If no candidate is present for a module, it falls
#' back to using all genes in that module. The chosen gene name (plus
#' `newnameSuffix`) becomes the new label for that module level.
#'
#' Typical usage: rename gene modules using a list of transcription factors.
#'
#' @param df A data.frame (or tibble) with gene identifiers in `rownames(df)`.
#' @param candidates Character vector of candidate gene names.
#' @param ModuleNameCol Name of the column in `df` containing module labels.
#'   Will be converted to a factor if needed.
#' @param interestCol Name of the numeric column used to choose the "best" gene
#'   (max value) within each module.
#' @param newNameCol Name of an additional column to store the renamed module
#'   labels. If equal to `ModuleNameCol`, only `ModuleNameCol` is updated.
#' @param newnameSuffix Suffix appended to the chosen best gene to form the new
#'   module label.
#'
#' @return `df` with:
#' \itemize{
#'   \item `oldModuleName`: the original module labels
#'   \item updated factor levels in `ModuleNameCol`
#'   \item `newNameCol` (if different from `ModuleNameCol`): the renamed labels
#' }
#' Additionally, an attribute `level_map` (data.frame) is attached, giving the
#' old level, new level, and best gene per module.
#'
#' @export
#'
#' @examples
#' df <- data.frame  (
#'   Module = rep(c("M1","M2"), each = 5), Contribution = seq_len(10),
#'   row.names = c("TFAP2A", "SOX10", "TYR", "DCT", "PMEL", "OCA2", "SLC24A5", "SLC45A2", "MC1R", "MITF")
#' )
#' candidates <- c("TFAP2A", "MITF", "SOX10")
#' print(renameLevelsByBestVal(df, candidates))
#'
renameLevelsByBestVal <- function(
		df,
		candidates,
		ModuleNameCol = "Module",
		interestCol   = "Contribution",
		newNameCol    = paste0(ModuleNameCol, "_newName"),
		newnameSuffix = ".Mod"
) {
	stopifnot(is.data.frame(df))

	if (is.null(rownames(df)) || anyNA(rownames(df)) || any(duplicated(rownames(df)))) {
		stop("`df` must have non-missing, unique `rownames(df)` (gene identifiers).")
	}
	if (!ModuleNameCol %in% colnames(df)) {
		stop("`ModuleNameCol` not found in `df`: ", ModuleNameCol)
	}
	if (!interestCol %in% colnames(df)) {
		stop("`interestCol` not found in `df`: ", interestCol)
	}

	# Ensure module labels are a factor (stable levels)
	module <- df[[ModuleNameCol]]
	if (!is.factor(module)) module <- as.factor(module)

	score <- df[[interestCol]]
	if (!is.numeric(score)) {
		stop("`interestCol` must be numeric (got class: ", paste(class(score), collapse = "/"), ").")
	}

	candidates <- unique(as.character(candidates))
	candidates <- intersect(candidates, rownames(df))

	modules <- levels(module)
	newLevels <- setNames(character(length(modules)), modules)
	bestGene  <- setNames(rep(NA_character_, length(modules)), modules)

	for (m in modules) {
		idx <- which(module == m)
		if (length(idx) == 0L) {
			newLevels[m] <- m
			next
		}

		genes <- rownames(df)[idx]
		ok <- !is.na(score[idx])

		genes_ok <- genes[ok]
		if (length(genes_ok) == 0L) {
			# No usable score in that module: keep original label
			newLevels[m] <- m
			next
		}

		pool <- intersect(genes_ok, candidates)
		if (length(pool) == 0L) pool <- genes_ok

		pool_idx <- match(pool, rownames(df))
		pool_scores <- score[pool_idx]

		best <- pool[which.max(pool_scores)]
		bestGene[m]  <- best
		newLevels[m] <- paste0(best, newnameSuffix)
	}


	# Apply renaming to ModuleNameCol
	renamed <- factor(module, levels = modules, labels = unname(newLevels))
	df[[newNameCol]] <- renamed
	df
}



#' Set the precise random seed state
#' @param state Object saved by getRandState
#' @export
#' @seealso getRandState
#' @return NULL, set the random seed state
#' @examples
#' randSate<-getRandState()
#' rnorm(5)
#' setRandState(randSate)
#' rnorm(5)
setRandState <- function(state) {
    # Assigning `NULL` state might lead to unwanted consequences
    if (!is.null(state)) {
        assign(".Random.seed",
            state,
            envir = .GlobalEnv,
            inherits = FALSE)
    }
}

#' Coefficient of variation
#'
#' @param x Numeric vector
#'
#' @return A single numeric value
#' @export
#' @examples
#' cv(c(1,2,3,4))
cv <- function(x) {
    return(sd(x) / mean(x))

}

#' Coefficient of variation of squared mean and sd
#'
#' @param x Numeric vector
#'
#' @return  A single numeric value
#' @export
#' @examples
#' cv2(c(1,2,3,4))
cv2 <- function(x) {
    return(sd(x) ^ 2 / mean(x) ^ 2)

}

#' Standard mean error
#'
#' @param x Numeric vector
#'
#' @return A single numeric value
#'
#' @export
#' @examples
#' se(c(1,2,3,4))
se <- function(x) {
    #
    return(sd(x) / sqrt(length(x)))

}


#' Transform vector to have no negative value (new minimum is 0)
#'
#' @param x Numeric vector
#'
#' @return Numeric vector
#' @export
#' @examples
#' uncenter(-5:5)
uncenter <- function(x) {
    #transform vector to have no negative value
    return(x + abs(min(x)))

}

#' Take first element of multiple values in a vector
#'
#' @description Similar to unique but conserve vector names or return index
#' where you can find each first value of multiple element.
#'
#' @param x Vector.
#' @param returnIndex Logical. Should the index of first elements or vector of
#'   first elements.
#'
#' @return Vector of first elements or numeric vector of index.
#' @export
#' @examples
#' a<-c(1,2,3,3,4,3)
#' names(a)<-c("a","b","c","d","e","f")
#' takefirst(a)
#' takefirst(a,returnIndex = TRUE)
takefirst <- function(x, returnIndex = FALSE) {
    uniqDat <- unique(x)
    caseUniq <- c()
    for (i in uniqDat)
        caseUniq <- c(caseUniq, which(i == x)[1])
    if (returnIndex) {
        return(caseUniq)
    } else{
        return(x[caseUniq])
    }
}


#' Compute the mode of a distribution.
#'
#' @source https://github.com/benmarwick/LaplacesDemon/blob/master/R/Mode.R
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value.
#' @export
#' @examples
#' Mode(c(seq_len(10),3))
#'
Mode <- function(x) {
    ### Initial Checks
    if (missing(x))
        stop("The x argument is required.")
    if (!is.vector(x))
        x <- as.vector(x)
    x <- x[is.finite(x)]
    ### Discrete
    if (all(x == round(x))) {
        Mode <- as.numeric(names(which.max(table(x))))
    } else {
        ### Continuous (using kernel density)
        x <- as.vector(as.numeric(as.character(x)))
        kde <- density(x)
        Mode <- kde$x[kde$y == max(kde$y)]
    }
    return(Mode)
}


#' Copy paste ready vector
#'
#' @param x A vector.
#'
#' @return A string ready to be copied and embedded as R code.
#' @export
#' @examples
#' copyReadyVector(seq_len(5))
copyReadyVector <- function(x) {
    paste0("c('", paste0(x, collapse = "','"), "')")
}


#' Better make.unique
#'
#' @description Similar to make.unique, but also add a sequence member for the
#' first encountered duplicated element.
#'
#' @param sample.names Character vector
#' @param sep A character string used to separate a duplicate name from its
#'   sequence number.
#'
#' @return A character vector of same length as names with duplicates changed.
#' @export
#' @examples
#' make.unique2(c("a", "a", "b"))
make.unique2 <- function(sample.names, sep = ".") {
    # Create a table of occurrences for each name
    tab <- table(sample.names)

    # Initialize a vector to store the new sample names
    newSampleNames <- character(length(sample.names))

    # Loop through each unique name only once
    for (name in names(tab)) {
        indices <- which(sample.names == name)
        # For names that occur more than once, append a suffix
        if (tab[name] > 1) {
            suffixes <- paste0(sep, seq_len(tab[name]))
            newSampleNames[indices] <- paste0(name, suffixes)
        } else {
            newSampleNames[indices] <- paste0(name,sep,1)
        }
    }

    newSampleNames
}

#' String split with chosen returned element
#'
#' @param x character vector, each element of which is to be split. Other
#'   inputs, including a factor, will give an error.
#' @param split character vector (or object which can be coerced to such)
#'   containing regular expression(s) (unless fixed = TRUE) to use for
#'   splitting. If empty matches occur, in particular if split has length 0, x
#'   is split into single characters. If split has length greater than 1, it is
#'   re-cycled along x.
#' @param n Single integer, the element index to be returned
#' @param fixed logical. If TRUE match split exactly, otherwise use regular
#'   expressions. Has priority over perl.
#' @param perl logical. Should Perl-compatible regexps be used?
#' @param useBytes logical. If TRUE the matching is done byte-by-byte rather
#'   than character-by-character, and inputs with marked encoding are not
#'   converted. This is forced (with a warning) if any input is found which is
#'   marked as "bytes" (see Encoding).
#'
#' @return A vector of the same length than x, with the n-th element for the
#'   split of each value.
#' @export
#' @examples
#' strsplitNth(c("ax1","bx2"), "x",1)
#' strsplitNth(c("ax1","bx2"), "x",2)
strsplitNth <-
    function(x,
            split,
            n = 1,
            fixed = FALSE,
            perl = FALSE,
            useBytes = FALSE) {
        res <- strsplit(x, split, fixed, perl, useBytes)
        vapply(res, function(el) {
            el[n]
        }, character(1))
    }


#' Convert numeric to string, add 0 to the number to respect lexicographical
#' order.
#'
#' @param x A numeric vector.
#' @param digit A single integer value. The maximum number of digits in the
#'   number sequence. It will determine the number of 0 to add.
#'
#' @return A charactervector.
#' @export
#' @examples
#' formatNumber2Character(seq_len(10))
#' formatNumber2Character(seq_len(10),digit = 4)
formatNumber2Character <-
    function(x, digit = max(nchar(as.character(x)))) {
        x <- format(x, scientific = FALSE, trim = TRUE)
        vapply(as.list(x), function(el) {
            paste0(paste0(rep("0", digit - nchar(el)), collapse = ""), el)
        }, character(1))
    }

#' Convert a named factor vector to a list
#'
#' @param factorValues A vector of factor. It has to be named if
#'   `factorNames=NULL`.
#' @param factorNames A character vector for providing the names separately.
#'
#' @return A list. Each element is named by a factor level of `factorValues`,
#'   and contains the provided names that had this level has a value.
#' @export
#' @examples
#' x<-factor(c("a","a","b","b","c","c","c"))
#' names(x)<-paste0("x",seq_len(7))
#' factorToVectorList(x)
#'
#' @seealso VectorListToFactor
factorToVectorList <- function(factorValues, factorNames = NULL) {
    if (is.null(factorNames))
        factorNames <- names(factorValues)
    if (is.character(factorValues))
        factorValues <- as.factor(factorValues)
    res <-
        lapply(levels(factorValues), function(x)
            factorNames[factorValues == x])
    names(res) <- levels(factorValues)
    res
}



#alias

#' colnames alias (getter)
#'
#' @param x A matrix-like object.
#' @param do.NULL logical. If `FALSE` and names are `NULL`, names are created.
#' @param prefix for created names
#' @return A character vector.
#' @export
#' @examples
#' data(iris)
#' cn(iris)
cn <-
    function(x, do.NULL = TRUE, prefix = "col"){
        colnames(x, do.NULL = TRUE, prefix = "col")
    }


#' colnames alias (setter)
#'
#' @param x A matrix-like object.
#' @param value Either NULL or a character vector equal of length equal to the
#'   appropriate dimension.
#'
#' @return The modified object.
#' @export
#'
#' @examples
#' data(iris)
#' cn(iris)<-c("S.l","S.w", "P.l", "P.w", "Sp")
#' cn(iris)
'cn<-' <- function(x, value) {
    colnames(x) <- value
    x
}

#' rownames alias (getter)
#'
#' @param x A matrix-like object.
#' @param do.NULL logical. If `FALSE` and names are `NULL`, names are created.
#' @param prefix for created names
#' @return A character vector.
#' @export
#' @examples
#' data(iris)
#' rn(iris)
rn <-
    function(x, do.NULL = TRUE, prefix = "row"){
        rownames(x, do.NULL = TRUE, prefix = "row")
    }

#' rownames alias (setter)
#'
#' @param x A matrix-like object.
#' @param value Either NULL or a character vector equal of length equal to the
#'   appropriate dimension.
#'
#' @return The modified object.
#' @export
#'
#' @examples
#' data(iris)
#' rn(iris)<-paste0("f",nrow(iris) |> seq_len())
#' rn(iris)
'rn<-' <- function(x, value) {
    rownames(x) <- value
    x
}

#' length alias
#'
#' @param x A vector-like object.
#' @return An integer.
#' @examples len(1:10)
#' @export
len <- function(x){
    length(x)
}


#' intersect alias
#'
#' @param x A vector-like object.
#' @param y A vector-like object.
#' @param ... Additional arguments passed to `intersect`.
#' @return A vector-like object.
#' @examples inter(1:10,5:15)
#'
#' @export
inter <- function(x, y, ...){
    intersect(x, y, ...)
}


#' Copy data to clipboard (only on Windows)
#'
#' @param x A dataframe / matrix / table / vector to copy to the clipboard
#' @param vectSep Character used for separating vector
#' @param ... currently not used
#'
#' @returns Nothing! but copy the input to the clipboard,
#' ready to be paste somewhere else
#' @export
#'
#' @examples
#' v <- seq_len(3)
#' names(v) <- c("A", "B", "C")
#' to_clipboard(V)
#'
to_clipboard <- function(x, vectSep="\n",...) {
	if(.Platform$OS.type != "windows") stop("only available in Windows")
	if(inherits(x, "data.frame") ){
		cols        <- lapply(x, format, ...)
		header_text <- paste(names(x), collapse = "\t")
		body_text   <- do.call(paste, c(cols, sep = "\t", collapse = "\n"))
		all_text    <- paste0(header_text, "\n", body_text)
	} else if(inherits(x, "matrix") | inherits(x, "table")){
		hasRowNames <- !is.null(row.names(x))
		hasColNames <- !is.null(colnames(x))
		if(hasRowNames & hasColNames){
			all_text <- "\t"
		} else {
			all_text <- ""
		}
		if (hasColNames) {
			all_text <- paste0(all_text, paste(colnames(x), collapse = "\t"),"\n")
		}
		rows <- apply(x, 1, function(x) paste(x, collapse = "\t"))
		if(hasRowNames) rows <- paste0(row.names(x), "\t", rows)
		all_text <- paste0(all_text, paste(rows, collapse = "\n"))
		paste0(all_text, paste(apply(x, 1, function(x) paste(x, collapse = "\t")), collapse = "\n"))
	} else {
		all_text <- paste0(x, collapse = vectSep)
	}
	# Format the input as a tab-separated table

	# Write the formatted text to the system clipboard
	utils::writeClipboard(all_text)

	# Return the input, invisibly
	invisible(all_text)
}
