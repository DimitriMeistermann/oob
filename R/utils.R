#' Compute geometrical mean
#'
#' @param x numeric vector
#' @param keepZero get rid of 0 before computing.
#'
#' @return A single numeric value.
#' @examples
#' gmean(c(1,2,3))
#' gmean(c(0,2,3),keepZero = TRUE)
#' gmean(c(0,2,3),keepZero = FALSE)
gmean<-function(x, keepZero=FALSE){ #geometrical mean
	if(sum(x)==0) return(0)
	if(!keepZero){
		x<-x[x!=0]
	}else{
		if(length(which(x==0))>0) return(0)
	}
	return( exp( sum(log(x))/length(x) ) )
}

#' Get the precise random seed state
getRandState <- function() {
	# Using `get0()` here to have `NULL` output in case object doesn't exist.
	# Also using `inherits = FALSE` to get value exactly from global environment
	# and not from one of its parent.
	get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

#' Set the precise random seed state
#' @param state Object saved by getRandState
setRandState <- function(state) {
	# Assigning `NULL` state might lead to unwanted consequences
	if (!is.null(state)) {
		assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
	}
}

#' Coefficient of variation
#'
#' @param x Numeric vector
#'
#' @return A single numeric value
#'
#' @examples
#' cv(c(1,2,3,4))
cv<-function(x){
	return(sd(x)/mean(x));
}

#' Coefficient of variation of squared mean and sd
#'
#' @param x Numeric vector
#'
#' @return  A single numeric value
#'
#' @examples
#' cv2(c(1,2,3,4))
cv2<-function(x){
	return(sd(x)^2/mean(x)^2);
}

#' Standard mean error
#'
#' @param x Numeric vector
#'
#' @return A single numeric value
#'
#'
#' @examples
#' se(c(1,2,3,4))
se<-function(x){ #
	return(sd(x)/sqrt(length(x)));
}


#' Transform vector to have no negative value (new minimum is 0)
#'
#' @param x Numeric vector
#'
#' @return Numeric vector
#' @examples
#' uncenter(-5:5)
uncenter<-function(x){
	#transform vector to have no negative value
	return(x+abs(min(x)));
}

#' Take first element of multiple values in a vector
#'
#' @description
#' Similar to unique but conserve vector names or return index where you can find each first value of multiple element.
#'
#' @param x Vector.
#' @param returnIndex Logical. Should the index of first elements or vector of first elements.
#'
#' @return Vector of first elements or numeric vector of index.
#'
#' @examples
#' a<-c(1,2,3,3,4,3)
#' names(a)<-c("a","b","c","d","e","f")
#' takefirst(a)
#' takefirst(a,returnIndex = TRUE)
takefirst<-function(x,returnIndex=FALSE){
	uniqDat<-unique(x)
	caseUniq<-c()
	for(i in uniqDat) caseUniq<-c(caseUniq,which(i==x)[1])
	if(returnIndex){
		return(caseUniq)
	}else{
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
#'
#' @examples Mode(c(1:10,3))
#'
Mode <- function(x) {
	### Initial Checks
	if(missing(x)) stop("The x argument is required.")
	if(!is.vector(x)) x <- as.vector(x)
	x <- x[is.finite(x)]
	### Discrete
	if(all(x == round(x))) {
		Mode <- as.numeric(names(which.max(table(x))))
	} else {
		### Continuous (using kernel density)
		x <- as.vector(as.numeric(as.character(x)))
		kde <- density(x)
		Mode <- kde$x[kde$y == max(kde$y)]}
	return(Mode)
}


#' Copy paste ready vector
#'
#' @param x A vector.
#'
#' @return A string ready to be copied and embedded as R code.
#'
#' @examples
#' copyReadyVector(1:5)
copyReadyVector<-function(x){
	paste0("c('",paste0(x,collapse = "','"),"')")
}


#' Better make.unique
#'
#' @description
#' Similar to make.unique, but also add a sequence member for the first encountered duplicated element.
#'
#' @param sample.name Character vector
#' @param sep A character string used to separate a duplicate name from its sequence number.
#'
#' @return A character vector of same length as names with duplicates changed.
#'
#' @examples
#' make.unique2(c("a", "a", "b"))
make.unique2<-function(sample.name,sep="."){
	counts=table(as.factor(sample.name))
	nmRep<-sapply(as.list(counts),function(x) 1:x)
	paste0(rep(names(nmRep),counts),sep,as.character(unlist(nmRep,use.names = F)))[order(sample.name)[order(sample.name)]]
}

#' String split with chosen returned element
#'
#' @param x character vector, each element of which is to be split. Other inputs, including a factor, will give an error.
#' @param split character vector (or object which can be coerced to such) containing regular expression(s) (unless fixed = TRUE) to use for splitting. If empty matches occur, in particular if split has length 0, x is split into single characters. If split has length greater than 1, it is re-cycled along x.
#' @param n Single integer, the element index to be returned
#' @param fixed logical. If TRUE match split exactly, otherwise use regular expressions. Has priority over perl.
#' @param perl logical. Should Perl-compatible regexps be used?
#' @param useBytes logical. If TRUE the matching is done byte-by-byte rather than character-by-character, and inputs with marked encodings are not converted. This is forced (with a warning) if any input is found which is marked as "bytes" (see Encoding).
#'
#' @return A vector of the same length than x, with the n-th element for the split of each value.
#'
#' @examples
#' strsplitNth(c("ax1","bx2"), "x",1)
#' strsplitNth(c("ax1","bx2"), "x",2)
strsplitNth<-function(x, split, n=1, fixed=FALSE, perl=FALSE, useBytes=FALSE){
	res<-strsplit(x, split, fixed, perl, useBytes)
	sapply(res,function(el){ el[n] })
}


#' Convert numeric to string, add 0 to the number to respect lexicographical order.
#'
#' @param x A numeric vector.
#' @param digit A single integer value. The maximum number of digits in the number sequence. It will determine the number of 0 to add.
#'
#' @return A charactervector.
#'
#' @examples
#' formatNumber2Character(1:10)
#' formatNumber2Character(1:10,digit = 4)
formatNumber2Character<-function(x,digit=max(nchar(as.character(x)))){
	x<-as.character(x)
	sapply(as.list(x),function(el){ paste0(paste0(rep("0",digit-nchar(el)),collapse = ""),el) })
}


#' Convert a named factor vector to a list
#'
#' @param factorValues A vector of factor. It has to be named if `factorNames=NULL`.
#' @param factorNames A character vector for providing the names separately.
#'
#' @return A list. Each element is named by a factor level of `factorValues`, and contains the provided names that had this level has a value.
#'
#' @examples
#' x<-factor(c("a","a","b","b","c","c","c"))
#' names(x)<-paste0("x",1:7)
#' factorToVectorList(x)
#'
#' @seealso VectorListToFactor
factorToVectorList<-function(factorValues,factorNames=NULL){
	if(is.null(factorNames)) factorNames<-names(factorValues)
	res<-lapply(levels(factorValues),function(x) factorNames[factorValues==x])
	names(res)<-levels(factorValues)
	res
}

#'  Convert a list to a named factor vector
#'
#' @param listOfVector A named list. Each element must contain a character vector.
#'
#' @return A named factor vector.
#' @export
#'
#' @examples
#' VectorListToFactor(list(a=c("x1","x2"),b=c("x3","x4"),c=c("x5","x6","x7")))
#'
#' @seealso factorToVectorList
VectorListToFactor<-function(listOfVector){
	res<-factor(unlist(lapply(seq_along(listOfVector),function(i) rep(names(listOfVector)[i],length(listOfVector[[i]])))),
							levels=names(listOfVector))
	names(res)<-unlist(listOfVector)
	res
}


#' Transform a range of value to another by a linear relationship.
#'
#' @description
#' Similar to the javascript function `d3.scaleLinear()`.
#'
#' @param vals A numeric vector. Values to be transposed in the new range.
#' @param newRange A vector of two numeric values: the new minimum and maximum.
#' @param returnFunction Logical. Return the linear scale as a function instead of the transposed values in a new scale. If set to `TRUE`, `vals` argument can be also a vector of 2 numeric corresponding to the minimum and maximum of the old range.
#' @return A vector of value or a function if `returnFunction=TRUE`.
#' @export
#'
#' @examples
#' oldValues<-1:10
#' linearScale(oldValues,c(0,1),returnFunction = FALSE)
#' scaleFun<-linearScale(c(1,10),c(0,1),returnFunction = TRUE)
#' scaleFun(oldValues)
linearScale <- function(vals,newRange,returnFunction = TRUE) {
	if(!is.numeric(vals)) stop("x should be a vector of numerics")
	if(length(newRange)!=2 | !is.numeric(newRange)) stop("newRange should be a vector of 2 numerics")

	oldMin<-min(vals)
	oldMax<-max(vals)
	newMin<-newRange[1]
	newMax<-newRange[2]

	mfac<-(newMax-newMin)/(oldMax-oldMin)
	scaleFun<-function(x) newMin+(x-oldMin)*mfac

	if(returnFunction){
		scaleFun
	}else{
		scaleFun(vals)
	}
}





