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

#' Very simple over abundance analysis via logistic regression
#'
#' @param qualVar A vector or factor or character.
#' @param df A DataFrame, same number of rows as length of qualVar.
#' @param cols Vector of character. Column names of df to be used in the regression model.
#' @param removeIntercepts Should the intercept be removed from the model? Useful if odds ratio / pvalue of all levels of a factor imncluded in cols are wanted.
#'
#' @return A dataframe containing 2 column per level of qualVar (Odds Ratio and pvalues).
#' @export
#'
#'
#' @examples
#' data(sampleAnnot)
#' simpleAbundanceAnalysis(sampleAnnot$culture_media, sampleAnnot,  c("line","run"))
#'

simpleAbundanceAnalysis<-function (varToExplain, df, cols, removeIntercepts = TRUE)
{
    varToExplain <- as.factor(varToExplain)
    levels(varToExplain) <- make.names(levels(varToExplain))
    lvls <- levels(varToExplain)
    dfReg <- df[, cols]
    formula <- ""
    for (col in cols) {
        if (is.numeric(dfReg[, col])) {
            formula <- paste0(formula, "+", col)
        }
        else {
            dfReg[, col] <- as.factor(dfReg[, col])
            if (removeIntercepts) {
                formula <- paste0(formula, "+", col, "-1")
            }
            else {
                formula <- paste0(formula, "+", col)
            }
        }
    }
    formula <- substr(formula, 2, nchar(formula))
    res <- lapply(lvls, function(lvl) {
        dfReg[, lvl] <- varToExplain == lvl
        m <- glm(formula = formula(paste0(lvl, "~", formula)),
                         data = dfReg, family = binomial, )
        resLm <- coef(summary(m))[, c(1, 4)]
        colnames(resLm) <- paste0(lvl, "_", c("OD", "pvalue"))
        resLm
    })
    names(res) <- lvls
    rown<-lapply(res,rownames) |> unlist() |> unique() |> sort()
    dtres<-data.frame(row.names = rown)
    for(lvl in lvls){
        odCol<-cn(res[[lvl]])[1]
        pvalCOl<-cn(res[[lvl]])[2]
        dtres[,odCol]<-NA
        dtres[,pvalCOl]<-NA
        dtres[rn(res[[lvl]]),odCol]<-res[[lvl]][,odCol]
        dtres[rn(res[[lvl]]),pvalCOl]<-res[[lvl]][,pvalCOl]
    }
    dtres
}


