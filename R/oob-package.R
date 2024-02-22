#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import basilisk
#' @import BiocParallel
#' @import ComplexHeatmap
#' @import ggbeeswarm
#' @import ggplot2
#' @import ggrepel
#' @import grDevices
#' @import grid
#' @import Rcpp
#' @import RcppArmadillo
#' @import rgl
#' @import scattermore
#' @import stringr
#' @import stats
#' @importFrom BiocGenerics
#'	colnames
#'  rownames
#'  intersect
#'  union
#'  setdiff
#'  sizeFactors
#' @importFrom rlang
#'  abort
#'  as_label
#'  caller_arg
#'  caller_env
#'  enquos
#'  env_get_list
#'  eval_tidy
#'  new_quosure
#'  quo_is_missing
#'  quo_is_call
#'  quo_is_symbolic
#'  quo_get_expr
#'  is_list
#'  is_symbolic
#'  is_quosure
#'  is_missing
#'  is_na
#'  is_null
#'  is_character
#'  is_vector
#' @importFrom batchelor fastMNN
#' @importFrom circlize colorRamp2
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom graphics axis
#' @importFrom labeling extended
#' @importFrom lsa cosine
#' @importFrom MASS bandwidth.nrd
#' @importFrom Matrix image
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom reticulate import
#' @importFrom SingleCellExperiment colData
#' @importFrom SingleCellExperiment sizeFactors
#' @importFrom SummarizedExperiment assay
#' @importFrom utils write.table
#' @useDynLib oob, .registration=TRUE
## usethis namespace: end
NULL

#' sampleAnnot
#'
#' @name sampleAnnot
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description
#' Sample Annotation sheet (often referred as colData) for 77 bulk RNA-Seq samples used in
#' Kilens, Meistermann et al. 2018.
#'
#'
#' @format A data.frame containing 77 observations for 18 features.
NULL

#' bulkLogCounts
#'
#' @name bulkLogCounts
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description
#' Log-expression count table for 77 samples sequenced from a bulk 3'digital RNA-Seq
#' protocol used in Kilens, Meistermann et al. 2018.
#' @format A matrix containing 77 observations for 16959 features/genes
NULL


#' DEgenesPrime_Naive
#'
#' @name DEgenesPrime_Naive
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description
#' Results of the differential expression analysis between primed and naive
#' human pluripotent stem cells in Kilens, Meistermann et al. 2018.
#' @format A data.frame containing 16959 observations for 7 features
NULL

#' geneLengthGRCh38
#' @name geneLengthGRCh38
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description
#' Gene length table for 25017 genes from the GRCh38 annotation
#'  used in Kilens, Meistermann et al. 2018.
#' @format A vector of 25017 numeric named with genes.
NULL

#' humanGeneIDtable
#' @name humanGeneIDtable
#' @docType data
#' @author Dimitri Meistermann
#' @references \url{https://doi.org/10.1038/s41467-017-02107-w}
#' @keywords data
#' @description
#' Gene ID table corresponance between ENTREZ ID, gene SYMBOL and ENSEMBL ID
#'  from the GRCh38 annotation. Used in Kilens, Meistermann et al. 2018.
#' @format A data frame containing 199476 rows for 3 columns
NULL
