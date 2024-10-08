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

#' Convert a list to a named factor vector
#'
#' @param listOfVector A named list. Each element must contain a character
#'   vector.
#'
#' @return A named factor vector.
#' @export
#'
#' @examples
#' VectorListToFactor(list(a=c("x1","x2"),b=c("x3","x4"),c=c("x5","x6","x7")))
#'
#' @seealso factorToVectorList
VectorListToFactor <- function(listOfVector) {
    res <-
        factor(unlist(lapply(seq_along(listOfVector), function(i)
            rep(names(listOfVector)[i], length(listOfVector[[i]])))),
            levels = names(listOfVector))
    names(res) <- unlist(listOfVector)
    res
}


#' Transform a range of value to another by a linear relationship.
#'
#' @description Similar to the javascript function `d3.scaleLinear()`.
#'
#' @param vals A numeric vector. Values to be transposed in the new range.
#' @param newRange A vector of two numeric values: the new minimum and maximum.
#' @param returnFunction Logical. Return the linear scale as a function instead
#'   of the transposed values in a new scale. If set to `TRUE`, `vals` argument
#'   can be also a vector of 2 numeric corresponding to the minimum and maximum
#'   of the old range.
#' @return A vector of value or a function if `returnFunction=TRUE`.
#' @export
#'
#' @examples
#' oldValues<-seq_len(10)
#' linearScale(oldValues,c(0,1),returnFunction = FALSE)
#' scaleFun<-linearScale(c(1,10),c(0,1),returnFunction = TRUE)
#' scaleFun(oldValues)
linearScale <- function(vals, newRange, returnFunction = TRUE) {
    if (!is.numeric(vals))
        stop("x should be a vector of numerics")
    if (length(newRange) != 2 |
        !is.numeric(newRange))
        stop("newRange should be a vector of 2 numerics")

    oldMin <- min(vals)
    oldMax <- max(vals)
    newMin <- newRange[1]
    newMax <- newRange[2]

    mfac <- (newMax - newMin) / (oldMax - oldMin)
    scaleFun <- function(x)
        newMin + (x - oldMin) * mfac

    if (returnFunction) {
        scaleFun
    } else{
        scaleFun(vals)
    }
}

#' Return ordered index of element with top n value.
#'
#' @param x Numeric vector.
#' @param top Integer. Number of element to be returned.
#' @param decreasing Logical, return top element by decreasing or increasing
#'   order.
#'
#' @return A vector of integer.
#' @export
#'
#' @examples
#' x<-c(1,5,6,10,5.2,3,8)
#' whichTop(x)
#' whichTop(x,top=3)
#' whichTop(x,decreasing=FALSE)
whichTop <- function(x, top = 5, decreasing = TRUE) {
    order(x, decreasing = decreasing)[seq_len(top)]
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

