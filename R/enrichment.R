#' Compute p-value and other metric for enrichment of set intersection.
#'
#' @param intersectionSize Numeric. Number of common element between sets.
#' @param setSizes Vector of numeric. Number of element in each set.
#' @param universeSize Total number of element in the universe.
#'
#' @details
#' If two sets, perform an hypergeometric test, if more the p-value is computed from binomial distribution.
#'
#' @return A list of values:
#' - observed: observed number of elements in intersection (same value as `intersectionSize`)
#' - expected: expected number of elements in intersection.
#' - OEdeviation: Effect size of the association observed - expected / sqrt(universeSize)
#' - pval: pvalue (upper tail) where NULL hypotheses is that observed intersection is from the distribution of intersection size with p = nExpected / nUniverse.
#'
#' @export
#'
#' @examples
#' setSizes<-c(150,200,250)
#' universeSize<-10000
#' intersectionSize<-5
#' enrichSetIntersection(intersectionSize, setSizes, universeSize)
enrichSetIntersection<-function(intersectionSize, setSizes, universeSize){
    nSet<-length(setSizes)
    expected<-prod(setSizes) / (universeSize^(nSet-1))
    if(nSet==2){
        pval<-phyper(q = intersectionSize-0.5, m = setSizes[1],n = universeSize-setSizes[1], k = setSizes[2], lower.tail=FALSE)
    }else if(nSet>2){
        pval<-pbinom(intersectionSize-1,universeSize, expected/universeSize,  lower.tail = FALSE)
    }else{
        stop("Number of set should be superior to 1")
    }
    OEdeviation<- (intersectionSize - expected) / sqrt(universeSize)
    return(list("observed"=intersectionSize,"expected"=expected,"OEdeviation"=OEdeviation,"pval"=pval))
}


